"""
Helper script to run several analyses
Usage:
python run_all config.json path_to_genome_files CPUS
"""
from glob import glob
from sys import argv
from concurrent import futures
import ctdg_finder
import json

d_base = json.load(open(argv[1]))

def run_all(chrom):
    """
    This runs an analysis only changing the GFF and genome files
    It assumes that both files are in the same directory, and that
    the GFF ends in 'genes.gff' while the genome file ends in 'chroms.txt'
    :param chrom: path to the genome file
    :return:
    """
    gff = chrom.replace("chroms.txt", "genes.gff")
    # out = gff.replace("genes.gff", "clusters.gff")
    d_base["chroms"] = chrom
    d_base["gff"] = gff

    ctdg_finder.run(d_base)


all_chroms = glob(argv[2] + "/*.txt")
cpu = int(argv[3])
with futures.ProcessPoolExecutor(cpu) as pool:
    res = pool.map(run_all, all_chroms)
_ = list(res)