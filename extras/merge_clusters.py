import os
import sys
from concurrent import futures
from glob import glob

import pandas as pd

"""
Argument 1: all genes annotation (genes_parsed.csv
Argumnet 2: genes table
Argument 3: numbers table
Argument 4: gene table output
Argument 5: numbers table output
"""

annotation = pd.read_csv(sys.argv[1], index_col=0)
genes = pd.read_csv(sys.argv[2], index_col=0)
numbers = pd.read_csv(sys.argv[3])

# get list of clustered genes
cl_accs = list(genes.loc[~(genes["order"].isin(["na_95", "na_ms", "0", 0, "0.0", 0.0]))].index)


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
    print("{}: {}".format(sp, chrom))
    updated_clusters, updated_numbers = [], []
    numbers_cl = numbers.loc[(numbers["species"] == sp) & (numbers["chromosome"] == chrom)]
    intervals = merge_ranges(numbers_cl[["start", "_end"]].values)
    cluster_counter = 0
    for interval in intervals:
        interval_df = annotation.loc[(annotation["species"] == sp) & (annotation["chromosome"] == chrom) &
                                     (annotation["start"] >= interval[0]) & (annotation["end"] <= interval[1])]
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
            interval_clustered["cluster"] = cluster_name
            interval_clustered.sort_values("start", inplace=True)
            interval_clustered["order"] = range(1, cluster_size + 1)
            # Update numbers table
            updated_numbers.append([sp, chrom, cluster_name, interval[0], interval[1],
                                    cluster_size, interval_df.shape[0]])
            updated_clusters.append(interval_clustered)
            interval_not_clustered = interval_not_clustered.loc[~(interval_not_clustered.index.isin(interval_clustered.index))]
        updated_clusters.append(interval_not_clustered)
    if len(updated_clusters) > 0:
        updated_chrom = pd.concat(updated_clusters)
    else:
        updated_chrom = pd.DataFrame()

    numbers_df = pd.DataFrame(updated_numbers, columns=["species", "chromosome", "cluster", "start", "_end",
                                                        "duplicates", "total_genes"])
    numbers_df.loc[:, "proportion_duplicates"] = numbers_df["duplicates"] / numbers_df["total_genes"]
    return (updated_chrom, numbers_df)


with futures.ProcessPoolExecutor() as pool:
    updated_chroms = pool.map(update_cluster, pd.np.unique(genes.set_index(["species", "chromosome"]).index))


updated_results = list(updated_chroms)

updated_sp_list, updated_numbers_list = [], []
for cl, num in updated_results:
    updated_sp_list.append(cl)
    updated_numbers_list.append(num)
#
updated_sp = pd.concat(updated_sp_list)
updated_sp.to_csv(sys.argv[4])
updated_sp.loc[~(updated_sp["order"].isin(["0", "0.0", 0, 0.0]))].to_csv(sys.argv[4].replace(".csv", "_clean.csv"))
#
# # Update numbers table
updated_numbers = pd.concat(updated_numbers_list)
updated_numbers.to_csv(sys.argv[5])
updated_numbers.to_csv(sys.argv[5].replace(".csv", "_clean.csv"))
