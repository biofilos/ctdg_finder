import pandas as pd
import json
from sys import argv
import os
from concurrent import futures
"""
Build JSON dictionaries, for each species/chromosome, where the keys are each
gene in the chromosome, and the values is a list containing all the genes in
that chromosome that share at least one domain as the gene in the key
"""
db = argv[1]
cpus = int(argv[2])
# Load dictionary of Pfam hits (via HMMer)
pfams = json.load(open(db + "/hmmer/pfam.json"))
families = pfams
# Load gene annotation
genes = pd.read_csv(db + "/genes_parsed.csv")


def verbose(flag=1):
    if flag > 0:
        print("Number of gene families/domains: {}:".format(len(pfams)))
        print("Number of genes: {}".format(len(genes["acc"].unique())))
        accs_csv = set(genes["acc"].values)
        accs_json = set()
        for pfam in pfams:
            pfam_accs = pfams[pfam]
            if flag > 1:
                print("{} genes in gene family/domain {}".format(len(pfam_accs),pfam))
            accs_json.update(pfam_accs)
        print("Total number of genes in hmmer result: {}".format(len(accs_json)))
        in_both = accs_json.intersection(accs_csv)
        len_both = len(in_both)
        print("Number of genes in hmmer result that are also in genes_parsed: {}".format(len_both))
        if len_both != len(accs_json):
            print("""The following genes are in the hmmer output but are not
                  in the csv:""")
            for not_in_csv in accs_json - accs_csv:
                print(not_in_csv)
            print("""The following genes are in the csv but not in the hmmer output""")
            for not_in_json in accs_csv - accs_json:
                print(not_in_json)
        else:
            print("All genes are in both data structures")

def rec2dict_sp(sp):
    """
    sp_chrom: list, where the first element is species
    """
# for sp, chrom in set(genes.set_index(["species", "chromosome"]).index.values):
    sp_json_file = db + "/sp_hmmers/{}.json".format(sp)
    if not os.path.exists(sp_json_file):
        # Initialize dictionary of Pfam hits
        pfam_in_chrom = {}
        for chrom in genes.loc[genes["species"] == sp, "chromosome"].unique():
            print("Processing {} ({})".format(sp, chrom))
            # Set of genes in the chromosome
            chrom_genes = set(genes.loc[(genes["species"] == sp) &
                                        (genes["chromosome"] == chrom), "acc"].values)
            # Go through each gene
            for gene in chrom_genes:
                same_pfams = []
                for pfam in pfams:
                    # Set of genes for each pfam accession
                    pfam_set = set(pfams[pfam])
                    # If the reference gene has that domain
                    if gene in pfam_set:
                        # Add all the genes in the chromosome that have that domain
                        same_pfams += list(chrom_genes.intersection(pfam_set))
                # Add all those genes as a non-redundant list
                pfam_in_chrom[gene] = list(set(same_pfams))
        # Save to file
        chrom_out = open(sp_json_file, "w")
        chrom_out.write(json.dumps(pfam_in_chrom))
        chrom_out.close()
    print("*{} done".format(sp))

verbose(2)
with futures.ProcessPoolExecutor(cpus) as pool:
    res = pool.map(rec2dict_sp, list(set(genes.set_index("species").index.values)))
_ = list(res)
