import json
import os
import pybedtools
import sys
from collections import defaultdict
import gffutils

# Load configuration


def load_config(config_file):
    """
    Check that the config file has all the relevant parameters
    and that the files referenced in it exist, and return a dictionary
    
    Arguments:
        config_file {json file} -- JSON file containing relevant parameters
    """
    files = ["gff","chroms"]
    variables = ["feat_type"]
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
    

def filter_field(interval, kind, attr):
    """
    Return intervals of a specific type that
    has non-empy (with 'None') attribute
    
    Arguments:
        interval {BedTool} -- Interval
        interval {str} -- type (third column of gff)
        attr {str} -- attribute to be tested (from 9th col)
    """
    desired_type = kind == interval.fields[2]
    with_attr = attr in interval.attrs
    return desired_type and with_attr and "None" not in interval[attr]

def get_genes(config):
    """Extract the features to be clustered
    
    Arguments:
        config {json file} -- Json file containing parameters
    
    Returns:
        dict -- Dictionry of Bed tools for features with annotated families grouped by chromosome
    """
    # Load config
    gff = config["gff"]
    feat_type = config["feat_type"]
    out_dir = config["paths"]["out_dir"]
    species = gff.split(".gff")[0].split("/")[-1]    
    
    
    # Parse the gff
    feat_db_fn = "{}/{}_feats.db".format(out_dir, species)
    if os.path.exists(feat_db_fn):
        features = gffutils.FeatureDB(feat_db_fn)
    else:
        gff_str = ""
        for line in open(gff):
            data = line.split("\t")
            if line.startswith("##") or (data[2] == feat_type and not "None" in data[-1]):
                gff_str += line
        features = gffutils.create_db(gff_str, feat_db_fn, force=False, from_string=True)
    return features


def group_chrom(features):
    """
    Group features by chromosome
    
    Arguments:
        features {BedTool} -- Group of features
    Return:
        feat_dict {dict} -- Dictionary of features with chromosomes as keys
    """
    feat_dict = defaultdict(list)
    for feat in features:
        chrom = feat.chrom
        feat_dict[chrom].append(feat)
    return feat_dict

if __name__ == "__main__":
    # Check that arguments were passed
    assert len(sys.argv) > 1, "No configuration file was provided"
    config_file = sys.argv[1]
    config = load_config(config_file)    
    create_dirs(config["paths"])
    genes = get_genes(config)
    print("DONE")