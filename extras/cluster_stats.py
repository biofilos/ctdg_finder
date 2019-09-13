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

# def measure_clusters(genome, clusters):


def cluster_stats(genome, annotation, clusters, feature):
    all_genes = count_genes(annotation, feature)
    clu_genes = count_genes(clusters, feature)

    gene_table = pd.concat([pd.Series(x) for x in [all_genes, clu_genes]], 1).dropna()
    gene_table.columns = ("genes", "clustered_genes")
    gene_table.loc[:, "prop_clustered"] = gene_table["clustered_genes"] / gene_table["genes"]
    gene_table.sort_values("prop_clustered", inplace=True)
    return gene_table


if __name__ == "__main__":
    GENOME = argv[1]
    ANNOTATION = argv[2]
    CLUSTERS = argv[3]
    stats_table = cluster_stats(GENOME, ANNOTATION, CLUSTERS, FEATURE)
