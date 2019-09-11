import json
import os
import pybedtools
import sys
from collections import defaultdict, Counter
import gffutils
from HTSeq import GFF_Reader, GenomicInterval, GenomicArrayOfSets
from sklearn.cluster import MeanShift
import numpy as np
from ipdb import set_trace
# Load configuration


def load_config(config_file):
    """
    Check that the config file has all the relevant parameters
    and that the files referenced in it exist, and return a dictionary
    
    Arguments:
        config_file {json file} -- JSON file containing relevant parameters
    """
    files = ["gff","chroms"]
    variables = ["feat_type", "top_level_feat",
                 "samples","percentile_threshold", "overwrite"]
    paths = ["out_dir"]
    parameters = files + variables

    # 2. Load config file
    if type(config_file) == dict:
        config = config_file
    else:
        # 1. Check that the config file exists
        assert os.path.exists(config_file), config_file + " does not exist"
        config = json.load(open(config_file))
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
    for path in paths.values():
        os.makedirs(path, exist_ok=True)
    
def check_gff(config):
    """
    Quality check of GFF file
    
    Arguments:
        config {json file} -- json file containing relevant parameters
    """
    gff = config["gff"]
    records = [x.strip().split("\t") for x in open(gff)]
    

# def filter_field(interval, kind, attr):
#     """
#     Return intervals of a specific type that
#     has non-empy (with 'None') attribute
    
#     Arguments:
#         interval {BedTool} -- Interval
#         interval {str} -- type (third column of gff)
#         attr {str} -- attribute to be tested (from 9th col)
#     """
#     desired_type = kind == interval.fields[2]
#     with_attr = attr in interval.attrs
#     return desired_type and with_attr and "None" not in interval[attr]

def get_genes(gff, feat_type):
    """Extract the features to be clustered
    
    Arguments:
        gff {file_path} -- GFF file containing features to be clustered
        feat_type {str} -- feature to be clustered (genes, etc)
    
    Returns:
        generator of valid genes
    """
    
    # Yield features of the correct type, and with family annotation
    for rec in GFF_Reader(gff):
        if rec.type == feat_type and \
            "families" in rec.attr and \
            not "None" in rec.attr["families"]:
            yield rec

def get_chrom_info(gff, top_level):
    chroms = [rec for rec in GFF_Reader(gff) if rec.type in top_level]
    all_with_bw = all(["bandwidth" in rec.attr for rec in chroms])
    assert all_with_bw, "Some top-level features lack the bandwidth parameter"
    chrom_info = {rec.iv.chrom: float(rec.attr["bandwidth"]) for rec in chroms}
    return chrom_info
        
    
def group_chrom(gff, feature_type):
    chrom_dict = defaultdict(list)
    for feat in get_genes(gff, feature_type):
        feat_chrom = feat.iv.chrom
        chrom_dict[feat_chrom].append(feat)
    return chrom_dict
            

def get_chrom_homologs(gff,feat_type):
    genes_per_chrom = group_chrom(gff, feat_type)
    families = {x: defaultdict(list) for x in genes_per_chrom.keys()}
    for chrom, chrom_genes in genes_per_chrom.items():
        for gene in chrom_genes:
            for fam in gene.attr["families"].split(","):
                families[chrom][fam].append(gene)
    return families


def meanshift(gene_list, chrom_info):
    # Process only lists with more than one gene
    if len(gene_list) == 1:
        return [0]
    else:
        # Get bandwidth from one of the genes
        chromosome = gene_list[0].iv.chrom
        bandwidth = chrom_info[chromosome]
        xy = [(gene.iv.start, 0) for gene in gene_list]
        MS = MeanShift(bandwidth=bandwidth).fit(xy)
        return MS.labels_


def genes_in_interval(record, feature, chrom, start, end):
    is_feat = record.fields[2] == feature
    in_chrom = record.chrom == chrom
    inside_interval = record.start < end and record.end > start
    return is_feat and in_chrom and inside_interval

def cluster_stats(gff, feature, chrom_lengths, cluster_start,
                  cluster_end, samples, percentile):
    cluster_len = cluster_end - cluster_start            
    # Genome intervals
    
    # Random intervals
    random_bedtool = pybedtools.BedTool()
    random_intervals = random_bedtool.random(l=cluster_len,
                                                n=samples,
                                                g=chrom_lengths)
    all_max_dups = []
    for rand_interval in random_intervals:
        
        n_dups = 0
        r_start = rand_interval.start
        r_end = rand_interval.end
        r_chrom = rand_interval.chrom
        sample_interval = GenomicInterval(r_chrom, r_start, r_end)
        potential_dups = [x.attr["families"] for x in GFF_Reader(gff) if sample_interval.overlaps(x.iv) and x.type == feature and not "None" in x.attr["families"]]
        
        # for x in pybedtools.BedTool(gff):
        #     if genes_in_interval(x, feature, r_chrom,r_start, r_end) and x.attrs["families"] != "None":
        #         potential_dups.append(x.attrs["families"])
        if potential_dups:
            family_counter = Counter(potential_dups)
            all_max_dups.append(max(family_counter.values()))
            perc_samples = np.percentile(all_max_dups, percentile)
            return perc_samples
        else:
            return 0
def sort_genes(gene_rec):
    return gene_rec.iv.start


def cluster_with_stats(cluster_family, chrom_info, gff, feature,chrom_lenghts, samples,percentile, cluster_counter):
    gff_str = ""
    gene = cluster_family[0]
    cluster_n_genes = len(cluster_family)
    cluster_chrom = gene.iv.chrom
    gene_starts = [x.iv.start for x in cluster_family]
    gene_ends = [x.iv.end for x in cluster_family]
    cluster_start = min(gene_starts)
    cluster_end = max(gene_ends)
    cluster_len = cluster_end - cluster_start
    percentile_thresh = cluster_stats(gff, feature, chrom_lenghts,
                                        cluster_start, cluster_end,
                                        samples, percentile)
    # Check if the cluster has more duplicates than the sample
    if cluster_n_genes >= percentile_thresh:
        families = {x.attr["families"] for x in cluster_family}
        cl_families = ",".join(families)
        cluster_name = "{}_{}".format(cluster_chrom, cluster_counter)
        attrs = "ID={};length={};duplicates={};percentile={};families={}\n".format(cluster_name, cluster_len,
                                                                                    cluster_n_genes,
                                                                                    percentile_thresh,cl_families)
        cluster_gff = [cluster_chrom, "CTDGFinder", "cluster",
                        str(cluster_start), str(cluster_end), ".", ".", ".",
                        attrs]
        cluster_gff_line = "\t".join(cluster_gff)
        gff_str += cluster_gff_line
        sp = gff.split("/")[-1].replace(".gff", "")
        print("{}\t{}: {} duplicates (P95: {})".format(sp, cluster_name, cluster_n_genes, percentile_thresh))
        for ix, gene in enumerate(sorted(cluster_family, key=sort_genes)):
            gene.attr["Parent"] = cluster_name
            gene.attr["order"] = ix + 1
            del gene.attr["meanshift_cluster"]
            gene_gff_pre_line = gene.get_gff_line()
            gene_line = gene_gff_pre_line.replace("; ", ";").replace(" ", "=").replace('"', "")
            gff_str += gene_line
    return gff_str


def clustering(feature, chrom_homologs, chrom_info, gff, chrom_lenghts, samples, percentile):
    gff_str = ""
    for chrom, fam_dict in chrom_homologs.items():
        cluster_counter = 1
        fams_for_removal = []
        for family, genelist in fam_dict.items():
            meanshift_labels = meanshift(genelist, chrom_info)
            # Remove meanshift clusters if they contain only one gene
            cluster_n_genes = len(genelist)
            if cluster_n_genes == 1:
                fams_for_removal.append((chrom,family))
            else:
                gene_starts, gene_ends = [0] * cluster_n_genes, [0] * cluster_n_genes
                # cluster_families = []
                cluster_groups = defaultdict(list)
                for ms_cluster, gene in zip(meanshift_labels, genelist):
                    ms_cluster_name = "{}_{}".format(chrom,ms_cluster)
                    gene.attr["meanshift_cluster"] = ms_cluster_name
                    cluster_groups[ms_cluster_name].append(gene)
                    # cluster_families += gene.attr["families"].split(',')
                for ms_cl in cluster_groups:
                    cluster_gene_list = cluster_groups[ms_cl]
                    # Skip meanshift clusters with only one gene
                    if len(cluster_gene_list) > 1:
                        cluster_gff = cluster_with_stats(cluster_gene_list, chrom_info,
                                                        gff, feature,
                                                        chrom_lenghts, samples,
                                                        percentile, cluster_counter)
                        cluster_counter += 1
                        gff_str += cluster_gff
                    

    return gff_str
def merge_clusters(out_path, feature_to_cluster):
    # Merge overlapping clusters
    merged_path = out_path.replace("clusters", "merged_clusters")
    with open(merged_path, "w") as f:
        merged = pybedtools.BedTool(out_path).filter(lambda x: x.fields[2]=="cluster").sort().merge()
        chrom_cts = defaultdict(int)
        for new_cluster in merged:
            cl_chrom = new_cluster.chrom
            cl_start = new_cluster.start
            cl_end = new_cluster.end
            chrom_cts[cl_chrom] += 1
            new_name = "{}_{}".format(cl_chrom, chrom_cts[cl_chrom])
            clust_genes = sorted([x for x in GFF_Reader(out_path) if x.type == feature_to_cluster and x.iv.chrom == cl_chrom and x.iv.start < cl_end and x.iv.end > cl_start], key=lambda x: x.iv.start)
            # write cluster info
            f.write("{}\tCTDGFinder\tcluster\t{}\t{}\t.\t.\t.\tID={};duplicates={}\n".format(cl_chrom, cl_start, cl_end, new_name, len(clust_genes)))
            for ix, gene in enumerate(clust_genes):
                ct = ix + 1
                gene.attr["order"] = ct
                gene.attr["Parent"] = new_name
                new_gene_gff = gene.get_gff_line().replace("; ", ";").replace(" ", "=")
                f.write(new_gene_gff)
# def group_ms_clusters(meanshift_clusters):
#     for 
    
def run(config_file):
    # Check that arguments were passed
    
    config = load_config(config_file)
    gff = config["gff"]
    chrom_lens = config["chroms"]
    samples= config["samples"]
    percentile_threshold = config["percentile_threshold"]
    feature_to_cluster = config["feat_type"]
    top_level_feat = config["top_level_feat"]
    create_dirs(config["paths"])
    out_suffix = "_clusters.gff"
    outfile = gff.replace("_genes.gff", out_suffix).replace("_genes.gff3", out_suffix).split("/")[-1]
    out_dir = config["paths"]["out_dir"]
    out_path = "{}/{}".format(out_dir, outfile)
    if os.path.exists(out_path) and not config["overwrite"]:
        print("Not overwriting {}".format(out_path))
    else:
        print("Working on {}".format(gff))
        # Get groups of homologous genes per chromosome
        chrom_homologs = get_chrom_homologs(gff, feature_to_cluster)
        # Extract bandwidth parameter for MeanShift for all the chromosomes
        chrom_info = get_chrom_info(gff, top_level_feat)
        clusters_gff = clustering(feature_to_cluster, chrom_homologs, chrom_info,
                                        gff, chrom_lens, samples, percentile_threshold)

        with open(out_path, "w") as f:
            f.write(clusters_gff)
        merge_clusters(out_path, feature_to_cluster)
                    
            
    print("DONE")

if __name__ == "__main__":
    assert len(sys.argv) > 1, "No configuration file was provided"
    config_file = sys.argv[1]
    run(config_file)