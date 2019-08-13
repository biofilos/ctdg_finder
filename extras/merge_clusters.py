import os
import sys
from concurrent import futures
from glob import glob
from collections import namedtuple
import pandas as pd

"""
Argument 1: all genes annotation (genes_parsed.csv
Argumnet 2: CTDG_out directory
Argument 3: gene table output
Argument 4: numbers table output
Argument 5: CPUS
"""
####
def collect_results(initial_folder):
    # Initialize list of gene and cluster summary tables list
    m_genes, m_numbers = [], []
    # Go through all the output folders for each run
    for folder in glob("{}/*".format(initial_folder)):
        # Extract family name
        pfam = folder.split("/")[-1]
        # Format gene and table file names
        g_file = "{}/report/{}_genes.csv".format(folder, pfam)
        n_file = g_file.replace("genes", "numbers")
        # If the files exist and they are valid tables, read them
        try:
            g_df = pd.read_csv(g_file)
            #g_df.sort_values(["species", "chromsome", "start"], inplace=True)
            n_df = pd.read_csv(n_file)
            # If they contain data, include them: The gene table should
            # contain at least two genes, and the cluster summary table
            # at least one cluster
            if len(g_df.index) > 1 and len(n_df.index) > 0:
                m_genes.append(g_df)
                m_numbers.append(n_df)
        except:
            pass
    # Concatenate all genes and cluster summary tables
    all_genes = pd.concat(m_genes)
    all_numbers = pd.concat(m_numbers)
    # return genes and cluster summary tables
    return all_genes, all_numbers

genes, numbers = collect_results(sys.argv[2])
genes.loc[:, "chromosome"] = genes["chromosome"].astype(str)
numbers.loc[:, "chromosome"] = numbers["chromosome"].astype(str)
genes.sort_values(["species", "chromosome", "start"], inplace=True)
numbers.dropna(subset=["species"], inplace=True)
numbers = numbers.loc[~(numbers["cluster"].isin(["na_95", "na_ms", "0", 0,
                                                 "0.0", 0.0]))]
numbers.to_csv("all_numbers.csv")
genes.dropna(subset=["species"], inplace=True)
genes.set_index("acc").to_csv("all_genes.csv")
####


annotation = pd.read_csv(sys.argv[1])
annotation.loc[:, "chromosome"] = annotation["chromosome"].astype(str)
annotation.set_index("acc", inplace=True)
annotation.sort_values(["species", "chromosome", "start"], inplace=True)
# genes = pd.read_csv(sys.argv[2])
genes.set_index("acc", inplace=True)
# numbers = pd.read_csv(sys.argv[3])

# get list of clustered genes
cl_accs = list(genes.loc[~(genes["cluster"].isin(["na_95", "na_ms", "0", 0,
                                                  "0.0", 0.0]))].index.unique())

def merge_ranges(intervals):
    """
    Merge overlapping intervals
    :param intervals: list
    :return: merged list
    """
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged


# for sp, chrom in pd.np.unique(genes.set_index(["species", "chromosome"]).index):


def update_cluster(sp_chrom):
    sp, chrom = sp_chrom
    #sp = sp_chrom.species
    #chrom = sp_chrom.chromosome
    print("{}: {}".format(sp, chrom))
    updated_clusters, updated_numbers = [], []
    numbers_cl = numbers.loc[(numbers["species"] == sp) & (numbers["chromosome"] == chrom)]
    intervals = merge_ranges(numbers_cl[["start", "end"]].values)
    cluster_counter = 0
    for interval in intervals:
        interval_df = annotation.loc[(annotation["species"] == sp) &
                                     (annotation["chromosome"] == chrom) &
                                     (annotation["start"] >= interval[0]) &
                                     (annotation["end"] <= interval[1])]
        interval_df.sort_values("start", inplace=True)
        interval_df["cluster"] = "0"
        interval_df["order"] = "0"
        intervals_accs_clustered = set(interval_df.index).intersection(set(cl_accs))
        interval_clustered = interval_df.loc[intervals_accs_clustered]
        interval_not_clustered = interval_df.loc[~(interval_df.index.isin(intervals_accs_clustered))]
        cluster_size = interval_clustered.shape[0]
        if cluster_size > 1:
            cluster_counter += 1
            cluster_name = "{}_{}".format(chrom, cluster_counter)
            interval_clustered.loc[:, "cluster"] = cluster_name
            interval_clustered.sort_values("start", inplace=True)
            interval_clustered.loc[:, "order"] = range(1, cluster_size + 1)
            # Update numbers table
            updated_numbers.append([sp, chrom, cluster_name,
                                    interval[0], interval[1],
                                    cluster_size, interval_df.shape[0]])
            updated_clusters.append(interval_clustered)
            interval_not_clustered = interval_not_clustered.loc[~(interval_not_clustered.index.isin(interval_clustered.index))]
        updated_clusters.append(interval_not_clustered)
    if len(updated_clusters) > 0:
        updated_chrom = pd.concat(updated_clusters)
    else:
        updated_chrom = pd.DataFrame()

    numbers_df = pd.DataFrame(updated_numbers, columns=["species", "chromosome", "cluster", "start", "end",
                                                        "duplicates", "total_genes"])
    numbers_df.loc[:, "proportion_duplicates"] = numbers_df["duplicates"] / numbers_df["total_genes"]
    return updated_chrom, numbers_df

# Get unique list of species-chromosome pairs
sp_chroms = []
field = namedtuple("sp_chrom",["species","chromosome"])
for sp_chrom in genes[["species", "chromosome"]].values:
    sp_chrom_l = "{}XXXX{}".format(*list(sp_chrom))
    sp_chroms.append(sp_chrom_l)
    #sp_chroms.append(field(species=sp_chrom[0], chromosome=sp_chrom[1]))

non_overlap = [x.split("XXXX") for x in set(sp_chroms)]
with futures.ThreadPoolExecutor(int(sys.argv[5])) as pool:
    updated_chroms = pool.map(update_cluster, non_overlap)


updated_results = list(updated_chroms)

updated_sp_list, updated_numbers_list = [], []
for cl, num in updated_results:
    updated_sp_list.append(cl)
    updated_numbers_list.append(num)
#
updated_sp = pd.concat(updated_sp_list)
updated_sp.to_csv(sys.argv[3])
updated_sp.loc[~(updated_sp["order"].isin(["0", "0.0", 0, 0.0]))].to_csv(sys.argv[3].replace(".csv", "_clean.csv"))
#
# # Update numbers table
updated_numbers = pd.concat(updated_numbers_list)
updated_numbers.to_csv(sys.argv[4])
updated_numbers.to_csv(sys.argv[4].replace(".csv", "_clean.csv"))
