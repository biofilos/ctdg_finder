import argparse
from concurrent import futures
from glob import glob
from time import time

import os
import pandas as pd
import shutil
from functools import partial

####
import cgpFinder as cgp

# Define input data location
parser = argparse.ArgumentParser("Paralogs cluster annotation pipeline")
parser.add_argument("--name_family", "-n",
                    action='store',
                    help="Name of the gene family")

parser.add_argument("--ref_seq", "-r",
                    action='store',
                    help="reference gene family sequence")

parser.add_argument("--blast_samples", "-b",
                    action="store",
                    type=int,
                    help="Number of samples to build empirical distribution of paralogs")
parser.add_argument("--sp", "-s",
                    action="append",
                    help="Species to run the analysis on (must be written in quotes)")
parser.add_argument("--out_dir", "-o",
                    action="store",
                    default='cgp_out',
                    help="Output directory")
parser.add_argument("--cpu", "-c",
                    action="store",
                    type=int,
                    default=1,
                    help="CPU to be used")
parser.add_argument("--db", "-d",
                    action="store",
                    help="Directory with gene annotation and blast database")


# Check if blast is installed
blast_path = shutil.which('blastp')
assert bool(blast_path), "blastp is not installed or it is not in your PATH"

args = parser.parse_args()

# check if database exists
assert os.path.exists(args.db), "directory {} with database does not exist".format(args.db)
# Check if blast database exists
blast_files = sum([os.path.exists('{}/all_seqs.fa.'.format(args.db) + x) for x in ['phr', 'pin', 'psq']])
assert blast_files == 3, "Some of the blast database files with prefix all_seqs.fa do not exist"

out_dir = args.out_dir
name_family = args.name_family
# Remove unfinished analysis
if os.path.exists(name_family):
    os.removedirs(name_family)
    print("partial results for {} found, removing".format(name_family))
# If results directory already exists in output directory, finish
if os.path.exists("{}/{}".format(out_dir, name_family)):
    print("Results for {} are already saved in {}".format(name_family, out_dir))
else:

    genomes_file = '{}/chromosomes.csv'.format(args.db)
    assert os.path.exists(genomes_file), 'Chromosome information file does not exist'
    cgp.annotation.genomes = pd.DataFrame.from_csv(genomes_file)
    cgp.annotation.genomes = cgp.annotation.genomes.dropna()
    cgp.annotation.genomes['chromosome'] = cgp.annotation.genomes['chromosome'].astype(str)
    # Genes annotation
    genes_annotation_file = '{}/genes_parsed.csv'.format(args.db)
    assert os.path.exists(genomes_file), 'Genes annotation file does not exist'
    cgp.annotation.all_genes = pd.DataFrame.from_csv(genes_annotation_file)
    cgp.annotation.all_genes['chromosome'] = cgp.annotation.all_genes['chromosome'].astype(str)
    # Remove genes in non-assembled chromosomes
    cgp.annotation.all_genes = cgp.annotation.all_genes.loc[~(cgp.annotation.all_genes.chromosome.isnull()) &
                                                             (cgp.annotation.all_genes.chromosome != 'NA')]
    cgp.annotation.all_genes['chromosome'] = cgp.annotation.all_genes['chromosome'].astype(str)
    # Genes sequences
    cgp.annotation.all_genes_fasta = '{}/all_seqs.fa'.format(args.db)
    # Reference sequence for Blast
    ref_seq = args.ref_seq
    assert os.path.exists(ref_seq), 'Reference sequence does not exist'
    # Name of gene family
    cgp.annotation.name_family = name_family
    # Number of samples for the empirical distributions
    blast_samples_size = args.blast_samples
    # Check that the species is in the allowed species
    if args.sp:
        species = args.sp
        valid_species = set(cgp.annotation.genomes.sp.values)
        valid_sp_list = "\n".join(valid_species)
        for sp in species:
            assert sp in valid_species, "\n{} is not a valid species\n"\
            "Please select species from the following list{}\n".format(sp, valid_sp_list)
    else:
        species = []
    # Number of available CPUs
    cpus = args.cpu

    pd.options.mode.chained_assignment = None
    init_time = time()
    # Run analysis
    # Create directory structure
    cgp.create_folders(name_family, ['blast_samples', 'blast_samples_gw', 'report'])
    # Run Blast
    cgp.blast.exe(cpus, ref_seq, cgp.annotation.all_genes_fasta, "{0}/report/{0}.blast".format(name_family))
    if len(species) > 0:
        print("Running analysis with species:")
        print("\n".join(species))
    else:
        print("Running analysis with all species")
    family_blast = cgp.blast.parse("{0}/report/{0}.blast".format(name_family), 2, True, species)
    family_blast['species'] = family_blast['species'].map(lambda x: str(x.replace(' ', '_')))
    if len(family_blast) > 0:
        # Remove genes in non-assembled chromosomes

        family_blast = family_blast.loc[family_blast.chromosome != 'NA']

        # Set a list of rows, so that the MeanShift section can be run in parallel
        blast_rows = [sptable for sptable in cgp.sp_loop(family_blast, ['species', 'chromosome'])]

        with futures.ProcessPoolExecutor() as pool:
            # Run MeanShift cluster annotation
            family_blast_list = pool.map(cgp.meanshift_cluster, blast_rows)
        family_blast = pd.concat(family_blast_list)

        # Same as above, exclude non-clustered genes
        only_clusters = family_blast.loc[~(family_blast['cluster'].isin(['na', '0', 0]))]
        only_cluster_rows = [sptable for sptable in cgp.sp_loop(only_clusters, ['species', 'chromosome', 'cluster'])]

        with futures.ProcessPoolExecutor() as pool:
            # Get number of paralogs, total genes and proportions of paralogs per cluster
            family_numbers_list = pool.map(cgp.gene_numbers, only_cluster_rows)
        # Polish table
        # Remove empty entries from gene_numbers()
        family_numbers_list = [x for x in family_numbers_list if x is not None]
        # Only run the rest of the analysis if there are any annotated clusters so far

        if len(family_numbers_list) > 0:
            family_numbers = pd.DataFrame(family_numbers_list, columns=['species', 'chromosome', 'cluster',
                                                                        'start', 'end', 'paralogs',
                                                                        'total_genes', 'proportion_paralogs'])
            cluster_rows = [inTab for inTab in cgp.sp_loop(family_numbers, ['species', 'chromosome', 'cluster']) if
                            inTab is not None]

            # Partly fill required arguments on these functions, so that map() can be run only
            # with the table list
            one_arg_blast_samples = partial(cgp.blast_sampling, gw=False, samples=blast_samples_size, db=args.db)
            one_arg_blast_samples_gw = partial(cgp.blast_sampling, gw=True, samples=blast_samples_size, db=args.db)
            # Run the sampling algorithm
            print("Analyzing {} proto-cluster(s)".format(len(cluster_rows)))
            print("{:<25} {:<9} {:<25} {}".format('species','paralogs','proto-cluster','sample (95P)'))
            with futures.ProcessPoolExecutor(cpus) as pool:
                # Only picking regions from the same chromosome
                chrom_wide = pool.map(one_arg_blast_samples, cluster_rows)
                genome_wide = pool.map(one_arg_blast_samples_gw, cluster_rows)
            # chrom_wide = map(one_arg_blast_samples, cluster_rows)
            # genome_wide = map(one_arg_blast_samples_gw, cluster_rows)
            for chrom_data, genome_data in zip(chrom_wide, genome_wide):
                family_numbers.loc[(family_numbers.species == chrom_data[0]) &
                                   (family_numbers.cluster == chrom_data[1]), 'perc95_chrom'] = chrom_data[2]

                family_numbers.loc[(family_numbers.species == genome_data[0]) &
                                   (family_numbers.cluster == genome_data[1]), 'perc95_gw'] = genome_data[2]
            # Identify clusters that didn't make it through the threshold
            for_deletion = family_numbers.loc[family_numbers['paralogs'] < family_numbers['perc95_gw'],
                                              ['species', 'chromosome', 'cluster']]
            # Get accession numbers for the genes that are going to be reaanotated
            accs_for_deletion = pd.merge(family_blast, for_deletion, how='inner',
                                         on=['species', 'chromosome', 'cluster'])['prot_acc']
            # Remove annotation for cluster name and cluster order
            family_blast.loc[family_blast['prot_acc'].isin(accs_for_deletion), ['cluster', 'order']] = 'na'
            family_numbers.loc[family_numbers['paralogs'] < family_numbers['perc95_gw'], 'cluster'] = 'na'

            # # # Save cluster tables
            family_numbers.to_csv("{0}/report/{0}_numbers.csv".format(name_family))
            family_numbers.loc[~(family_numbers['cluster'].isin(['na', '0', 0]))].to_csv(
                "{0}/report/{0}_numbers_clean.csv".format(name_family))
        else:
            print("NO CLUSTERS in the {} family".format(name_family))
        # # Save hmmer hits table with annotations of their clusters
        family_blast.to_csv("{0}/report/{0}_genes.csv".format(name_family))
        # Concatenate sample files in one table per set of samples

        chrom_files = glob("{}/blast_samples/*.samples".format(name_family))
        if len(chrom_files) > 0:
            chrom_tabs = pd.concat([pd.read_table(x, header=None, sep=',') for x in chrom_files])
            gw_files = glob("{}/blast_samples_gw/*.samples".format(name_family))
            gw_tabs = pd.concat([pd.read_table(x, header=None, sep=',') for x in gw_files])
            samples_summary = pd.concat([cgp.set_sample_tables(chrom_tabs, 'chromosome'),
                                         cgp.set_sample_tables(gw_tabs, 'genome_wide')])
            # save samples to report dir
            samples_summary.to_csv("{0}/report/{0}.samples".format(name_family))
            # Remove samples per sp
            _ = [os.remove(x) for x in chrom_files + gw_files]
        else:
            print("No clusters for {}".format(name_family))
        # Sequence analysis
#        if args.seq_analysis:
#            print("Downloading sequences")
#            cgp.phylo.download_genes(family_blast)
#            print("Aligning")
#            cgp.phylo.align()
    else:
        print("Not enough blast hits({}) to perform the analysis".format(len(family_blast)))
    # End
    # Delete intermediate files
    cgp.delete_intermediates(out_dir)
    print("DONE")
    run_time = time() - init_time
    seconds = int(run_time % 60)
    minutes = int(run_time / 60)
    hours = int(minutes / 60)
    minutes = minutes - hours * 60
    print("Results for {} were saved in {}".format(name_family, out_dir))
    print("Run time: {}:{}:{}\n".format(str(hours).zfill(2), str(minutes).zfill(2), str(seconds).zfill(2)))
