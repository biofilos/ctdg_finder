import pandas as pd
import json
from sys import argv
import os
from concurrent import futures

"""
Build JSON dictionaries, for each species/chromosome, where the keys are each
gene in the chromosome, and the values is a list containing all the genes in
that chromosome that share the same Pfam domain as the gene in the key
"""
db = argv[1]
cpus = int(argv[2])
# Load dictionary of Pfam hits (via HMMer)
pfams = json.load(open(db + "/hmmer/pfam.json"))
families = pfams
# Load gene annotation
genes = pd.read_csv(db + "/genes_parsed.csv")

# Go through each chromosome
# def rec2dict(sp_chrom):
#     """
#     sp_chrom: list, where the first element is species, and the second
#     a chromosome name
#     """
# # for sp, chrom in set(genes.set_index(["species", "chromosome"]).index.values):
#     sp, chrom = sp_chrom
#     if not os.path.exists(db + "/hmmers/{}_{}.json".format(sp, chrom)):
#         print(sp, chrom)
#         # Initialize dictionary of Pfam hits
#         pfam_in_chrom = {}
#         # Set of genes in the chromosome
#         chrom_genes = set(genes.loc[(genes["species"] == sp) &
#                                     (genes["chromosome"] == chrom),
#                                     "acc"].values)
#         # Go through each gene
#         for gene in chrom_genes:
#             same_pfams = []
#             for pfam in pfams:
#                 # Set of genes for each pfam accession
#                 pfam_set = set(pfams[pfam])
#                 # If the reference gene has that domain
#                 if gene in pfam_set:
#                     # Add all the genes in the chromosome that have that domain
#                     same_pfams += list(chrom_genes.intersection(pfam_set))
#             # Add all those genes as a non-redundant list
#             pfam_in_chrom[gene] = list(set(same_pfams))
#         # Save to file
#         chrom_out = open(db + "/hmmers/{}_{}.json".format(sp, chrom), "w")
#         chrom_out.write(json.dumps(pfam_in_chrom))
#         chrom_out.close()
#     else:
#         print("{} {} done".format(sp, chrom))

# with futures.ProcessPoolExecutor(cpus) as pool:
#     res = pool.map(rec2dict, list(set(genes.set_index(["species", "chromosome"]).index.values)))
# _ = list(res)


# def homologs_chrom(sp):
#     json_file = "{}/sp_hmmers/{}.json".format(db, sp)
#     print("Processing {}".format(sp))
#     if not os.path.exists(json_file+"d"):
#         sp_df = genes.loc[genes["species"] == sp]
#         families_dict = {gene: set() for gene in sp_df["acc"].values}
#         for family in families:
#             # print("Procesing {}".format(family))
#             genes_in_family = set(families[family])
#             for chrom in sp_df["chromosome"].unique():
#                 genes_in_chrom = set(sp_df.loc[sp_df["chromosome"] == chrom, "acc"].values)
#                 genes_in_family_chrom = genes_in_chrom.intersection(genes_in_family)
#                 for gene in genes_in_family_chrom:
#                     families_dict[gene].update(genes_in_family_chrom)
#         json_out = open(json_file, "w")
#         json_out.write(json.dumps(families_dict))
#         json_out.close()
#         print("* {} Done".format(sp))


# with futures.ProcessPoolExecutor(cpus) as pool:
#     res = pool.map(homologs_chrom, genes["species"].unique())
# _ = list(res)

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


with futures.ProcessPoolExecutor(cpus) as pool:
    res = pool.map(rec2dict_sp, list(set(genes.set_index("species").index.values)))
_ = list(res)
