import json
import os
import pybedtools
import sys
from collections import defaultdict, Counter
# import gffutils
from HTSeq import GFF_Reader, GenomicInterval
from sklearn.cluster import MeanShift
import numpy as np


# Load configuration


def load_config(config_json_dict):
    """
    Check that the config file has all the relevant parameters
    and that the files referenced in it exist, and return a dictionary
    
    Arguments:
        config_json_dict {str or dict} -- JSON file containing relevant parameters or dictionary of parsed JSON
    """
    files = ["gff", "chroms"]
    variables = ["feat_type", "top_level_feat",
                 "samples", "percentile_threshold",
                 "overwrite"]
    paths = ["out_dir"]
    parameters = files + variables

    # 2. Load config file
    if type(config_json_dict) == dict:
        config = config_json_dict
    else:
        # 1. Check that the config file exists
        assert os.path.exists(config_json_dict), config_json_dict + " does not exist"
        config = json.load(open(config_json_dict))
    # 3. Check that all the necessary parameters are specified 
    for parameter in parameters:
        assert parameter in config, parameter + " is missing in the config file"
    # 4. Check that all the required files exist
    for path in files:
        config_path = config[path]
        assert os.path.exists(config_path), config_path + " does not exist"
    # 5. Check that the paths to be created are specified
    for path in paths:
        assert path in config["paths"], path + " not specified"
    return config


def create_dirs(paths):
    """
    Create necesasary directories
    :param paths: list of directories to be created
    :return: None
    """
    for path in paths.values():
        os.makedirs(path, exist_ok=True)


def get_genes(gff, feat_type):
    """
    Extract features to be clustered
    :param gff: GFF containing the features to be clustered
    :param feat_type: feature to be clustered
    :return: generator of selected features
    """
    # Yield features of the correct type, and with family annotation
    for rec in GFF_Reader(gff):
        if rec.type == feat_type and "families" in rec.attr and "None" not in rec.attr["families"]:
            yield rec


def get_chrom_info(gff, top_level):
    """
    Parses bandwidth information from the GFF as a dictionary where the keys are top-level
    sequences, and the values are the bandwidth parameter calculated for that sequence
    :param gff: path name
    :param top_level: name of the top-level sequence (chromosome, scaffold, etc)
    :return: dict
    """
    # Parse features of type `top_level`
    chroms = [rec for rec in GFF_Reader(gff) if rec.type in top_level]
    # Check that all the top-level sequences have a bandwidth GFF attribute
    all_with_bw = all(["bandwidth" in rec.attr for rec in chroms])
    assert all_with_bw, "Some top-level features lack the bandwidth parameter"
    # Generate dictionary
    chrom_info = {rec.iv.chrom: float(rec.attr["bandwidth"]) for rec in chroms}
    return chrom_info


def group_chrom(gff, feature_type):
    """
    Group features of a specific type in chromosomes as a dictionary with
    keys = chromosomes and values a list of features belonging to that chromosome
    :param gff: path to GFF file
    :param feature_type: feature to be extracted
    :return: dictionary
    """
    # initialize dictionary
    chrom_dict = defaultdict(list)
    # Extract features of type `feature_type`
    for feat in get_genes(gff, feature_type):
        # Add the gene to its respective chromosome key
        feat_chrom = feat.iv.chrom
        chrom_dict[feat_chrom].append(feat)
    return chrom_dict


def get_chrom_homologs(gff, feat_type):
    """
    Group the features in each chromosome by the family they belong to
    :param gff: path to gff file
    :param feat_type: features to be extracted
    :return: dictionary: {chrom:{family1:[gene1,gene2],family2:[gene3,gene2]}}
    """
    # Get features grouped by chromosome
    genes_per_chrom = group_chrom(gff, feat_type)
    # Initialize dictionary
    families = {x: defaultdict(list) for x in genes_per_chrom.keys()}
    # Scan the genes in each chromosome
    for chrom, chrom_genes in genes_per_chrom.items():
        for gene in chrom_genes:
            # Add genes to the family (or families) they belong to
            for fam in gene.attr["families"].split(","):
                families[chrom][fam].append(gene)
        # Order the genes by start coordinate in case the GFF was not sorted
        for fam, genelist in families[chrom].items():
            families[chrom][fam] = sorted(genelist, key=lambda x: x.iv.start)
    return families


def meanshift(gene_list, chrom_info):
    """
    Cluster list of features using the MeanShift algorithm
    :param gene_list: features to be clustered (assumed to be in the same chromosome)
    :param chrom_info: dictionary {chromosome: bandwidth}
    :return: labels of the clusters
    """
    # Process only lists with more than one gene
    n_genes = len(gene_list)
    if n_genes == 1:
        return [0]
    else:
        # Here, I am going to subtract the length of n gene from
        # all the genes next to it. This to eliminate the
        # influence of gene lengths in the intergenic distance
        xy = [[0, 0] for x in gene_list]
        for ix in range(n_genes - 1):
            gene1 = gene_list[ix].iv
            for ix2 in range(ix + 1, n_genes):
                gene2 = gene_list[ix2].iv
                ig = gene2.start - gene1.end
                xy[ix + 1][0] = ig

        # Get bandwidth from one of the genes
        chromosome = gene_list[0].iv.chrom

        bandwidth = chrom_info[chromosome]
        # The MeanShift works on two-dimensional data, so a zero will be inputed as the "y coordinate"
        # xy = [(gene.iv.start, 0) for gene in gene_list]
        # Run the MeanShift algorithm
        ms = MeanShift(bandwidth=bandwidth).fit(xy)
        return ms.labels_


def genes_in_interval(record, feature, chrom, start, end):
    """
    Check if a feature is of a specific type and is in a specified iv
    (used for filtering)
    :param record: BedTool (feature)
    :param feature: feature type
    :param chrom: chromosome name
    :param start: start of the iv
    :param end: end of the iv
    :return: bool
    """
    # Check that the record is of the specified feature
    is_feat = record.fields[2] == feature
    # Check that the record is in the correct chromosome
    in_chrom = record.chrom == chrom
    # Check if the record overlaps the specified iv
    inside_interval = record.start < end and record.end > start
    # Return True if all conditions are True
    return is_feat and in_chrom and inside_interval


def cluster_stats(gff, feature, chrom_lengths, cluster_start,
                  cluster_end, samples, percentile):
    """
    Performs a statistical assesment of the clusteredness of a cluster:
    This takes a number of regions of the same length of the cluster
    and calculates the percentile X of the distribution of genes of the same
    family in each sample
    :param gff: path to gff file
    :param feature: feature type that was clustered
    :param chrom_lengths: genome file for bedtools
    :param cluster_start: start coordinate of the cluster
    :param cluster_end: end coordinate of the cluster
    :param samples: number of samples to be taken
    :param percentile: percentile threshold
    :return: percentile (int)
    """
    # Calculate the length of the cluster
    cluster_len = cluster_end - cluster_start
    # Generate `samples` number of random genomic region of length `cluster_len`
    random_bedtool = pybedtools.BedTool()
    random_intervals = random_bedtool.random(l=cluster_len, n=samples, g=chrom_lengths)
    # Initialize list of the distribution of features from the same family
    # all_max_dups = []
    # Scan random intervals
    # for rand_interval in random_intervals:
    #     # Get coordinates of the iv
    #     r_start = rand_interval.start
    #     r_end = rand_interval.end
    #     r_chrom = rand_interval.chrom
    #     # Use HTSeq data structures (maybe pybedtools can be faster
    #     sample_interval = GenomicInterval(r_chrom, r_start, r_end)
    #     # Extract the families of the features overlapping the random region
    #     potential_dups = [x.attr["families"].split(",") for x in GFF_Reader(gff) if
    #                       sample_interval.overlaps(x.iv) and x.type == feature and "None" not in x.attr["families"]]
    # f = open("debug.log", "a")
    # f.write("cl_chrom,cl_start,cl_end,sample,sample_chrom,sample_start,sample_end,feat_chrom,feat_start,feat_end\n")

    rand_ivs = [GenomicInterval(x.chrom, x.start, x.end) for x in random_intervals]
    potential_dups = [[] for x in range(samples)]
    for feat in GFF_Reader(gff):
        if feat.type == feature:
            for ix, sample_interval in enumerate(rand_ivs):
                if sample_interval.overlaps(feat.iv):
                    families = feat.attr["families"].split(",")
                    if "None" not in families:
                        potential_dups[ix] += families
                        # f.write("{},{},{},{},{},{},{},{},{},{}\n".format(cluster_chrom, cluster_start, cluster_end,
                        #                                         ix, sample_interval.chrom, sample_interval.start,
                        #                                         sample_interval.end, feat.iv.chrom,
                        #                                         feat.iv.start, feat.iv.end))
    # f.close()

    # Flatten the list of families
    # samples_log_file = "samples.log"
    # if os.path.exists(samples_log_file):
    #     samples_lst = json.load(open(samples_log_file))
    # else:
    #     samples_lst = []
    # ff = open(samples_log_file, "w")

    all_max_dups = [0 for x in potential_dups]
    for ix, sample in enumerate(potential_dups):
        if sample:
            n_dups = max(Counter(sample).values())
            if n_dups == 1:
                n_dups = 0
            all_max_dups[ix] = n_dups
    # potential_dups_flat = []
    # for i in potential_dups:
    #     potential_dups_flat += i
    # if potential_dups_flat:
    #     # Count how many featues have the same family
    #     family_counter = Counter(potential_dups_flat)
    #     # Add the maximum number of families represented in the sample region
    #     all_max_dups.append(max(family_counter.values()))
    # else:
    #     # If no features were included, use 0
    #     all_max_dups.append(0)
    # Calculate the X percentile
    perc_samples = np.percentile(all_max_dups, percentile)
    # log_data = {
    #     "{}_{}_{}".format(cluster_chrom,cluster_start,cluster_end): (all_max_dups, perc_samples)
    # }
    # samples_lst.append(log_data)
    # ff.write(json.dumps(samples_lst))
    # ff.close()
    # samples_lst.append({})
    return perc_samples


def sort_genes(gene_rec):
    """
    to be used by the `key` argument of the python builtin `sorted()` function
    to order features by their start coordinate
    :param gene_rec:
    :return: start coordinate
    """
    return gene_rec.iv.start


def cluster_with_stats(cluster_family, gff, feature, chrom_lenghts, samples, percentile, cluster_counter):
    """
    Performs statistical sampling and returns a GFF string of a cluster
    :param cluster_family: list of genes belonging to a cluster from MeanShift
    :param gff: path to the annotation file
    :param feature: feature to be clustered
    :param chrom_lenghts: genome file
    :param samples: number of samples for the statistical step
    :param percentile: percentile threshold to be used as inclusion criterion
    :param cluster_counter: index used for naming clusters
    :return:
    """
    # Initialize gff string
    gff_str = ""
    # Pick one gene to extract chromosome information
    gene = cluster_family[0]
    cluster_chrom = gene.iv.chrom
    # Check that all the genes are in the same chromosome
    chroms = {x.iv.chrom for x in cluster_family}
    in_the_chrom = len(chroms) == 1
    assert in_the_chrom, "Not all the genes are in the same chromosome ({})".format(chroms)
    # Calculate number of genes in the cluster
    cluster_n_genes = len(cluster_family)
    # Get list of starting and ending positions of the genes in the cluster
    gene_starts = [x.iv.start for x in cluster_family]
    gene_ends = [x.iv.end for x in cluster_family]
    # Calculate start, end and length of the cluster
    # IMPORTANT: HTSeq changes things to 0-based, so
    # I have to add one to all coordinates here
    cluster_start = min(gene_starts) + 1
    cluster_end = max(gene_ends) + 1
    cluster_len = cluster_end - cluster_start
    # Calculate percentile (statistical assesment)
    percentile_thresh = cluster_stats(gff, feature, chrom_lenghts,
                                      cluster_start, cluster_end,
                                      samples, percentile)
    # Check if the cluster has more duplicates than the sample
    # Take the first part of the GFF path as the species name
    sp = gff.split("/")[-1].replace("_genes.gff", "")
    cluster_name = "{}_{}".format(cluster_chrom, cluster_counter)
    good_cluster = cluster_n_genes >= percentile_thresh
    if good_cluster:
        msg = "*"
    else:
        msg = ""
    msg_cl_name = "...{}".format(cluster_name[-10:])
    chrom_info = "{}: {} - {}".format(msg_cl_name, cluster_start, cluster_end)
    print("{} {} D {} P95: {}{}".format(sp.ljust(30), chrom_info.ljust(35),
                                        str(cluster_n_genes).ljust(3),
                                        percentile_thresh.round(4), msg))
    if good_cluster:
        # Extract the names of the gene families represented in the cluster
        families = [x.attr["families"].split(",") for x in cluster_family]
        families_flat = []
        for f in families:
            families_flat += f
        families_flat = set(families_flat)
        # Encode the families in the cluster as a comma-separated string
        cl_families = ",".join(families_flat)
        # Format cluster name (chrom_index)

        # Encode the attributes in the GFF
        attrs = "ID={};length={};duplicates={};percentile={};families={}\n".format(cluster_name, cluster_len,
                                                                                   cluster_n_genes,
                                                                                   percentile_thresh, cl_families)
        # Format cluster information into a GFF string
        cluster_gff = [cluster_chrom, "CTDGFinder", "cluster",
                       str(cluster_start), str(cluster_end), ".", ".", ".",
                       attrs]
        cluster_gff_line = "\t".join(cluster_gff)
        # Add cluster to the GFF string
        gff_str += cluster_gff_line

        # Print progress

        # Process the genes in the cluster (ordered by start coordinate)
        for ix, gene in enumerate(sorted(cluster_family, key=sort_genes)):
            # Set the name of the cluster as the parent of the gene
            gene.attr["Parent"] = cluster_name
            # Set the order of the gene in the cluster
            gene.attr["order"] = ix + 1
            # Remove unused field
            del gene.attr["meanshift_cluster"]
            # Format GFF line
            gene_gff_pre_line = gene.get_gff_line()
            # Re-format GFF line to comply with GFF3
            # When returning the gff line, HTSeq works with 1-based coordinates, so changing them is not necessary
            gene_line = gene_gff_pre_line.replace("; ", ";").replace(" ", "=").replace('"', "")
            # Add gene to GFF string
            gff_str += gene_line
    return gff_str


def clustering(feature, chrom_homologs, chrom_info, gff, chrom_lenghts, samples, percentile, chrom_include):
    """
    Process a list of features grouped by chromosome and family and
    calls clusters
    :param chrom_include: Only analyze chromosomes in this list (mostly for debugging purposes)
    :param feature: feature to be clustered ('gene')
    :param chrom_homologs: dictionary of features from the same family {chrom:{fam1:[gene1,gene2]}}
    :param chrom_info: dictionary containing the bandwidth parameter calculated for each chromosome {chrom1:bandwidth}
    :param gff: path to gff file
    :param chrom_lenghts: genome file (tab-separated file of top-level sequences and their lengths
    :param samples: number of samples to be performed
    :param percentile: percentile threshold to call clusters
    :return: GFF3 string
    """

    # Initialize GFf string
    gff_str = ""
    # Parse out genes in each family in each chromosome
    for chrom, fam_dict in chrom_homologs.items():
        if len(chrom_include) == 0 or chrom in chrom_include:
            cluster_counter = 1
            fams_for_removal = []
            for family, genelist in fam_dict.items():
                # Perform MeanShift clustering
                meanshift_labels = meanshift(genelist, chrom_info)
                # Remove meanshift clusters if they contain only one gene
                cluster_n_genes = len(genelist)
                if cluster_n_genes == 1:
                    # This is not being used. Remove
                    fams_for_removal.append((chrom, family))
                else:
                    # gene_starts, gene_ends = [0] * cluster_n_genes, [0] * cluster_n_genes
                    # cluster_families = []
                    # Initialize groups of clusters
                    cluster_groups = defaultdict(list)
                    for ms_cluster, gene in zip(meanshift_labels, genelist):
                        # Assign a cluster ID (label from MeanShift) to the cluster candidate to its GF attributes field
                        ms_cluster_name = "{}_{}".format(chrom, ms_cluster)
                        gene.attr["meanshift_cluster"] = ms_cluster_name
                        # Add the gene to the cluster
                        cluster_groups[ms_cluster_name].append(gene)
                        # cluster_families += gene.attr["families"].split(',')
                    for ms_cl in cluster_groups:
                        cluster_gene_list = cluster_groups[ms_cl]
                        # Skip meanshift clusters with only one gene
                        if len(cluster_gene_list) > 1:
                            cluster_gff = cluster_with_stats(cluster_gene_list, gff, feature,
                                                             chrom_lenghts, samples,
                                                             percentile, cluster_counter)
                            cluster_counter += 1
                            gff_str += cluster_gff

    return gff_str


def merge_clusters(in_clusters, merged_clusters, feature_to_cluster):
    # Merge overlapping clusters
    # merged_clusters = in_clusters.replace("clusters", "merged_clusters")
    with open(merged_clusters, "w") as f:
        merged = pybedtools.BedTool(in_clusters).filter(lambda x: x.fields[2] == "cluster").sort().merge()
        chrom_cts = defaultdict(int)
        for new_cluster in merged:
            cl_chrom = new_cluster.chrom
            cl_start = new_cluster.start
            cl_end = new_cluster.end
            chrom_cts[cl_chrom] += 1
            new_name = "{}_{}".format(cl_chrom, chrom_cts[cl_chrom])
            clust_genes = sorted([x for x in GFF_Reader(in_clusters) if
                                  x.type == feature_to_cluster and x.iv.chrom == cl_chrom and
                                  x.iv.start < cl_end and x.iv.end > cl_start],
                                 key=lambda x: x.iv.start)
            # write cluster info

            cl_gff_str = ""
            cl_fams = []
            cl_start = []
            cl_end = []
            for ix, gene in enumerate(clust_genes):
                cl_fams += gene.attr["families"].split(",")
                cl_start.append(gene.iv.start)
                cl_end.append(gene.iv.end)
                ct = ix + 1
                gene.attr["order"] = ct
                gene.attr["Parent"] = new_name
                new_gene_gff = gene.get_gff_line().replace("; ", ";").replace(" ", "=")
                cl_gff_str += new_gene_gff
            cl_fams = ",".join(set(cl_fams))
            cl_start = min(cl_start)
            cl_end = max(cl_end)
            cl_length = cl_end - cl_start
            cl_attrs = "ID={};duplicates={};length={};families={}".format(new_name, len(clust_genes),
                                                                          cl_length, cl_fams)
            cl_summary_gff_str = "{}\tCTDGFinder\tcluster\t{}\t{}\t.\t.\t.\t{}\n".format(cl_chrom, cl_start, cl_end,
                                                                                         cl_attrs)
            f.write(cl_summary_gff_str)
            f.write(cl_gff_str)


def run(config_file_dict):
    # Check that arguments were passed

    config = load_config(config_file_dict)
    gff = config["gff"]
    chrom_lens = config["chroms"]
    samples = config["samples"]
    percentile_threshold = config["percentile_threshold"]
    feature_to_cluster = config["feat_type"]
    top_level_feat = config["top_level_feat"]
    create_dirs(config["paths"])
    out_suffix = "_clusters.gff"
    outfile = gff.replace("_genes.gff", out_suffix).replace("_genes.gff3", out_suffix).split("/")[-1]
    out_dir = config["paths"]["out_dir"]
    out_path = "{}/{}".format(out_dir, outfile)
    merged_clusters = out_path.replace("clusters", "merged_clusters")
    if os.path.exists(merged_clusters) and not config["overwrite"]:
        print("Not overwriting {}".format(out_path))
    else:
        print("Working on {}".format(gff))
        # Get groups of homologous genes per chromosome
        chrom_homologs = get_chrom_homologs(gff, feature_to_cluster)
        # Extract bandwidth parameter for MeanShift for all the chromosomes
        chrom_info = get_chrom_info(gff, top_level_feat)
        clusters_gff = clustering(feature_to_cluster, chrom_homologs, chrom_info,
                                  gff, chrom_lens, samples, percentile_threshold, chrom_include=[])

        with open(out_path, "w") as f:
            f.write(clusters_gff)
        merge_clusters(out_path, merged_clusters, feature_to_cluster)

    print("DONE")


if __name__ == "__main__":
    assert len(sys.argv) > 1, "No configuration file was provided"
    config_file = sys.argv[1]
    run(config_file)
