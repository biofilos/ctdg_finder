import json
from glob import glob

import numpy as np
import os
import shutil
import tarfile
from collections import Counter
from sklearn.cluster import MeanShift
from termcolor import colored

from cgpFinder import annotation as feat
from cgpFinder import blast
from cgpFinder import phylo


def sp_loop(in_table, columns):
    """
    Filters a table by a number of columns
    :param in_table: Input table
    :param columns: Column names to perform the filter
    :return: list of tables for each combination of column names
    """
    # Initialize the list of tables
    list_results = []
    for col_set in np.unique(in_table.set_index(columns).index):
        filtered_table = in_table
        cols_to_filter = {col_name: col_value for col_name, col_value in zip(columns, col_set)}
        for var in cols_to_filter.keys():
            filtered_table = filtered_table.loc[filtered_table[var] == cols_to_filter[var]]
        list_results.append(filtered_table)
    return list_results


def create_folders(root_dir, folders):
    """
    Create the directories where all the
    results will be saved
    """
    if not os.path.exists(root_dir):
        os.makedirs(root_dir)
        for sub_dir in folders:
            target = root_dir + '/' + sub_dir
            if not os.path.exists(target):
                os.makedirs(target)
    print("Directory tree for {} was created.".format(root_dir))


def meanshift_cluster(sp_table):
    """
    Extract clusters from annotation using the mean shift algorithm
    :param sp_table: hmmer hits table for each species/chromosome
    :return: annotated hmmer table
    """
    # Set variables
    sp_mean_shift = sp_table['species'].values[0]
    chrom_mean_shift = sp_table['chromosome'].values[0]
    proteome = feat.all_genes.loc[(feat.all_genes['species'] == sp_mean_shift) &
                             (feat.all_genes['chromosome'] == chrom_mean_shift)]
    # If there is only one gene, set the cluster to zero -> single-copy gene in the chromosome
    if len(sp_table) == 1:
        sp_table['cluster'] = 0
        sp_table['order'] = 0
    else:
        # Calculate the mean between neighboring genes in the chromosome
        gene_distances = proteome['end'].values[1:] - proteome['start'][:-1]
        mean_distance = np.mean(gene_distances)
        # Calculate standard deviation of the distance between neighboring genes in the chromosome
        sd_distance = np.std(gene_distances)
        # The band_width in MeanShift is set to the mean intergenic distance plus its standard deviation
        # This is more or less the minumum distance between "true" neighboring members of a cluster
        band_width = mean_distance + sd_distance
        # Since MeanShift works on 2D data, a 'y' axis was generated using zeroes
        gene_starts = np.array([(x, y) for x, y in zip(sp_table.start.values.astype(int),
                                                       np.zeros(len(sp_table)))])
        # Set the MeanShift
        ms = MeanShift(bandwidth=band_width)
        # Execute the MeanShift algorithm to the coordinates
        ms.fit(gene_starts)
        # The labels are the position of the members of each cluster in the coordinates array
        labels = ms.labels_
        labels_unique = np.unique(labels)
        # Number of clusters
        n_clusters = len(labels_unique)
        # Start from cluster 1
        cluster_number = 1
        # Go through each cluster label
        for k in range(n_clusters):
            # Set a masking array with the elements that belong to each cluster
            members = labels == k
            # Assign members of the cluster to new array
            proto_cluster = gene_starts[members, 0]
            cluster = np.sort(proto_cluster)
            # If there is more than one gene in the cluster
            if len(proto_cluster) > 1:
                # Annotate name (<name of the family>-<chromosome>_consecutive number
                sp_table.loc[sp_table['start'].isin(cluster),
                             'cluster'] = "{}-{}_{}".format(feat.name_family, chrom_mean_shift, cluster_number)
                # Annotate order of genes "from start of the chromosome, according to annotation
                sp_table.loc[sp_table['start'].isin(cluster),
                             'order'] = range(1, len(cluster) + 1)
                cluster_number += 1
            else:
                # If there is only one gene in the "cluster", assign name and order to zero
                sp_table.loc[sp_table['start'].isin(cluster),
                             'cluster'] = 'na'
                sp_table.loc[sp_table['start'].isin(cluster),
                             'order'] = 'na'
    # Return annotated table per species (useful for parallelization)
    return sp_table


def gene_numbers(sp_table):
    """
    Count the number of paralogs, total number of genes and proportion of paralogs
    for each cluster
    :param sp_table: hmmer table
    :return: annotated numbers table
    """
    sp_numbers = sp_table['species'].values[0]
    chromosome = sp_table['chromosome'].values[0]
    cluster = sp_table['cluster'].values[0]
    paralog_genes = len(sp_table)
    # If there is more than one paralog:
    if paralog_genes > 1:
        start = sp_table['start'].min()
        end = sp_table['end'].max()
        # Count the genes in that chromosomal region
        total_genes = len(feat.all_genes.loc[(feat.all_genes['species'] == sp_numbers) &
                                             (feat.all_genes['chromosome'] == chromosome) &
                                             (feat.all_genes['start'] >= start) &
                                             (feat.all_genes['end'] <= end)])
        # Calculate proportion of paralogs
        proportion_paralogs = paralog_genes / total_genes
        # Return the results per cluster (useful for parallelization
        gene_numbers_list = [sp_numbers, chromosome, cluster, start, end,
                             paralog_genes, total_genes, proportion_paralogs]
        return gene_numbers_list
    else:
        # If there are no genes in the cluster, return None
        print("No cluster to annotate for {}, chromosome {}".format(sp_numbers, chromosome))
        return None


def grab_paralogs(blast_dir, acc_list, sp, chrom, paralogs):
    """
    Given a list of accession numbers, get the maximum number of paralogs between them
    :param acc_list: accession numbers list
    :param paralogs: Dictionary with parsed blast results for the chromosome being sampled
    :return: max number of paralogs
    """
    paralogs_list = []
    # paralogs = json.loads(open('{}/{}_{}.json'.format(blast_dir, sp.replace(' ', '_'), chrom)).read())
    # If there are no genes in the sample interval, return 0 paralogs
    if len(acc_list) == 0:
        paralogs_max = 0
    # If there are more than one, return the maximum number of paralogs
    elif type(acc_list) == list:
        for acc in acc_list:
            # If the accession is not in the query, it had no blast hits with any sequence
            if acc in paralogs.keys():
                paralogs_list.extend(list(set(paralogs[acc]).intersection(acc_list)))
        # If more than zero paralogs were found
        if len(paralogs_list) > 0:
            paralogs_max = Counter(paralogs_list).most_common(1)[0][1]
        else:
            paralogs_max = 0
    else:
        paralogs_max = Counter(paralogs[acc_list]).most_common(1)[0][1]
    return paralogs_max


def blast_sampling(pre_cluster_table, gw, samples, db):
    sp = pre_cluster_table.loc[:, 'species'].values[0]
    # Load blast hits with queries from the speices
    sp_paralogs = json.loads(open("{}/blasts/{}.json".format(db, sp.replace(' ', '_'))).read())
    cluster = pre_cluster_table.loc[:, 'cluster'].values[0]
    cluster_length = abs(pre_cluster_table.loc[:, 'end'].values[0] - pre_cluster_table.loc[:, 'start'].values[0])
    # Set output directory depending of the type of sampling
    if gw:
        root_folder = '{}/blast_samples_gw'.format(feat.name_family)
        msg = 'G'
        col_name = 'perc95_gw'
    else:
        root_folder = '{}/blast_samples'.format(feat.name_family)
        msg = 'C'
        col_name = 'perc95_chrom'
    file_name = "{}/{}_{}".format(root_folder, sp.replace(' ', '_'), cluster)
    sample_coords_file = open("{}/{}_{}.coords".format(root_folder, sp.replace(' ', '_'), cluster), 'w')
    # with open(file_name, 'w') as seqO:
    sample_coords = []
    max_paras_list = []
    for i in range(1, samples + 1):
        if gw:
            chromosome = np.random.choice(feat.genomes.loc[feat.genomes.sp == sp, 'chromosome'].values)
        else:
            chromosome = pre_cluster_table.loc[:, 'chromosome'].values[0]

        chr_len = feat.genomes.loc[(feat.genomes.sp == sp) &
                              (feat.genomes.chromosome == chromosome), 'length'].values[0]
        # Set sample coordinates
        random_start = np.random.randint(1, chr_len)
        random_end = random_start + cluster_length

        # Save sample coordinates in log file
        sample_coords_file.write(','.join([str(i), str(sp), str(chromosome),
                                           str(random_start), str(random_end)]) + '\n')
        # Get genes in sample regions
        proteome_slice = feat.all_genes.loc[(feat.all_genes.species == sp) &
                                            (feat.all_genes.chromosome == chromosome) &
                                            (feat.all_genes.start >= random_start) &
                                            (feat.all_genes.end <= random_end)].index.values
        if len(proteome_slice) > 0:
            chrom_paralogs = sp_paralogs[chromosome]
            max_paralogs = grab_paralogs('blasts', list(proteome_slice), sp, chromosome, chrom_paralogs)
            current_sample = [i, sp, chromosome, random_start, random_end, max_paralogs]
            sample_coords.append(current_sample)
        else:
            max_paralogs = 0
        max_paras_list.append(max_paralogs)
    paralogs_log = open(file_name + ".samples", 'w')
    paralogs_data = ','.join([str(x) for x in max_paras_list])
    paralogs_log.write("{},{},{}\n".format(sp, cluster, paralogs_data))
    paralogs_log.close()
    try:
        original_paralogs = pre_cluster_table['paralogs'].values[0]
    except KeyError:
        print("Error in :{}, {}".format(sp, pre_cluster_table.cluster))
        return pre_cluster_table
    with_perc = pre_cluster_table.loc[:]
    with_perc[col_name] = np.percentile(max_paras_list, 95)
    query_perc = (sp, cluster, np.percentile(max_paras_list, 95))
    status_msg = "{:<25} {:<9} {:<25} {} ({})".format(sp, original_paralogs, cluster, msg, round(query_perc[2], 3))
    if gw and original_paralogs >= query_perc[2]:
        print(colored(status_msg, 'green'))
    elif gw and original_paralogs < query_perc[2]:
        print(colored(status_msg, 'red'))
    else:
        print(status_msg)
    return query_perc



def delete_intermediates(out_dir):
    """
    Remove intermediate files from hmmer sampling
    :return: None
    """
    sample_tar_o = tarfile.open('{}/report/blast_samples.tgz'.format(feat.name_family), 'w:gz')
    for sample_folder in ['blast_samples', 'blast_samples_gw']:
        # Process fastas
        for x in glob("{}/{}/*.fa".format(feat.name_family, sample_folder)):
            os.remove(x)
        # Process coordinates used for the sampling step
        for coord in glob("{}/{}/*.coords".format(feat.name_family, sample_folder)):
            sample_tar_o.add(coord, recursive=False)
            os.remove(coord)
        # Save blast outputs in tarfile, and remove the files
        for blast_file in glob("{}/{}/*.blast*".format(feat.name_family, sample_folder)):
            sample_tar_o.add(blast_file, recursive=False)
            os.remove(blast_file + "_filtered")
            os.remove(blast_file)
    sample_tar_o.close()
    for folder, _, __ in os.walk(feat.name_family):
        if len(os.listdir(folder)) == 0:
            os.removedirs(folder)
    # Move result to output directory
    shutil.move(feat.name_family, '{}/{}'.format(out_dir, feat.name_family))


def set_sample_tables(tabs, dataset):
    sample_cols = {-1: 'species', 0: 'cluster'}
    tabs = tabs.rename(columns=lambda x: x - 1 if type(x) != int else x)
    tabs = tabs.rename(columns=sample_cols)
    tabs['dataset'] = dataset
    return tabs

