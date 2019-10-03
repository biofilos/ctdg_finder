"""
Calculate basic statistics about the clusters in a genome
Usage:
python cluster_stats.py genome_file.txt annotation.gff clusters.gff
"""

from HTSeq import GFF_Reader as GFF
from sys import argv
from collections import defaultdict
import pandas as pd


FEATURE = "gene"


def count_genes(gff, feature):
    """
    Count the number of features per chromosome
    :param gff: path to GFF file
    :param feature: feature to be counted
    :return: dictionary {chrom: Number_of_features
    """
    done = [] # this might not be necessary,
    # but just in case there are genes represented twice
    gene_d = defaultdict(int)
    for rec in GFF(gff):
        if rec.type == feature:
            chrom = rec.iv.chrom
            gene = rec.name
            if gene not in done:
                gene_d[chrom] += 1
                done.append(gene)
    return gene_d


def measure_clusters(genome, clusters):
    """
    Calculate the total length of each top-level sequence
    that is clustered, and compare it with the length of
    their top-level sequences
    :param genome:
    :param clusters:
    :return: pd.DataFrame with total lengths of clusters and
    their top-level sequences, and the proportion of the
    chromosome length spanned by clusters
    """
    clusters_len = defaultdict(int)
    chromosomes_len = defaultdict(int)
    # Parse cluster lengths
    for feat in GFF(clusters):
        if feat.type == "cluster":
            clusters_len[feat.iv.chrom] += int(feat.attr["length"])
    # Parse chromosome lengths
    for line in open(genome):
        chrom, str_len = line.strip().split("\t")
        length = int(str_len)
        chromosomes_len[chrom] += length
    # Assemble measurements in dataframe
    cluster_tab = pd.concat([pd.Series(x) for x in [chromosomes_len, clusters_len]], 1)
    col_names = ("chromosome_length", "clustered_regions_length")
    cluster_tab.columns = col_names
    cluster_tab.loc[:, "prop_cluster_length"] = cluster_tab[col_names[1]] / cluster_tab[col_names[0]]
    return cluster_tab


def cluster_stats(genome, annotation, clusters, feature):
    """
    Assemble basic stats on numbers of clustered genes and the
    length of the clusters across the genome
    :param genome: path to genome file
    :param annotation: path to annotation gff file
    :param clusters: path to merged clusters gff file
    :param feature: feature to be clustered
    :return: pd.DataFrame
    """
    # Count genes
    all_genes = count_genes(annotation, feature)
    # Count clustered genes
    clu_genes = count_genes(clusters, feature)
    # Assemble table
    gene_table = pd.concat([pd.Series(x) for x in [all_genes, clu_genes]], 1).dropna()
    gene_table.columns = ("genes", "clustered_genes")
    gene_table.loc[:, "prop_clustered"] = gene_table["clustered_genes"] / gene_table["genes"]
    gene_table.sort_values("prop_clustered", inplace=True)
    # Include cluster/top-level sequence information
    cluster_tab = measure_clusters(genome, clusters)
    big_table = pd.concat([gene_table, cluster_tab], 1)
    return big_table.dropna()


if __name__ == "__main__":
    GENOME = argv[1]
    ANNOTATION = argv[2]
    CLUSTERS = argv[3]
    stats_table = cluster_stats(GENOME, ANNOTATION, CLUSTERS, FEATURE)
    #comment

