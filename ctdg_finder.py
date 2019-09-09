import json
import os
import pybedtools
import sys
from collections import defaultdict
import gffutils
from HTSeq import GFF_Reader
from sklearn.cluster import MeanShift

# Load configuration


def load_config(config_file):
    """
    Check that the config file has all the relevant parameters
    and that the files referenced in it exist, and return a dictionary
    
    Arguments:
        config_file {json file} -- JSON file containing relevant parameters
    """
    files = ["gff","chroms"]
    variables = ["feat_type", "top_level_feat"]
    paths = ["out_dir"]
    parameters = files + variables
    # 1. Check that the config file exists
    assert os.path.exists(config_file), config_file + " does not exist"
    # 2. Load config file
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
    chroms = [rec for rec in GFF_Reader(gff) if rec.type == top_level]
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
        return 0
    else:
        # Get bandwidth from one of the genes
        chromosome = gene_list[0].iv.chrom
        bandwidth = chrom_info[chromosome]
        xy = [(gene.iv.start, 0) for gene in gene_list]
        MS = MeanShift(bandwidth=bandwidth).fit(xy)
        return MS.labels_

if __name__ == "__main__":
    # Check that arguments were passed
    assert len(sys.argv) > 1, "No configuration file was provided"
    config_file = sys.argv[1]
    config = load_config(config_file)
    gff = config["gff"]
    feature_to_cluster = config["feat_type"]
    top_level_feat = config["top_level_feat"]
    create_dirs(config["paths"])
    # Get groups of homologous genes per chromosome
    chrom_homologs = get_chrom_homologs(gff, feature_to_cluster)
    # Extract bandwidth parameter for MeanShift for all the chromosomes
    chrom_info = get_chrom_info(gff, top_level_feat)
    
    print("DONE")