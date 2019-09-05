"""
Converts csv to gff and include annotation information 
in teh gff attributes
Usage
python csv2gff.csv chromosomes.json genes.csv {gene:fam1,fam2} dir_out
"""
import json
import os
import sys
from collections import namedtuple

import pandas as pd

#1. Parse arguments
_, CHROMS, GENES, FAMILIES, DIR_OUT = sys.argv

chroms = json.load(open(CHROMS))
families = json.load(open(FAMILIES))
genes = pd.read_csv(GENES)
# Correct strand notation
genes.loc[genes["strand"]== 1, "strand"] = "+"
genes.loc[genes["strand"]== -1, "strand"] = "-"
# add columns needed in the gff
genes.loc[:, "source"] = "Annotation"
genes.loc[:, "kind"] = "gene"
genes.loc[:, "dummy1"] = "."
genes.loc[:, "dummy2"] = "."
genes = genes.astype(str)
genes.loc[:, "chromosome"] = "chr" + genes["chromosome"].astype(str)
# Create directories if necessary
os.makedirs(DIR_OUT, exist_ok=True)
#2. Parse gff header and chromosome lines

def gene_attrs(gene_id):
    if gene_id in families:
        fams = ",".join(families[gene_id])
    else:
        fams = "None"
    attrs = "ID={};families={};biotype=protein_coding".format(gene_id, fams)
    return attrs

genes.loc[:, "attrs"] = genes["acc"].apply(gene_attrs)
# Sort genes in order of the gff
# genes = genes[["chromosome","source","kind","start","end","dummy1","strand","dummy2","attrs"]]
sp_lst = []
for sp, chrom_data in chroms.items():
    sp_lst.append(sp)
    gff_head = ["##gff-version 3"]
    gff_lst = []
    chrom_file = "{}/{}_chroms.txt".format(DIR_OUT, sp)
    gff_file = "{}/{}_genes.gff".format(DIR_OUT, sp)
    with open(chrom_file, "w") as chrom_f, open(gff_file, "w") as gff_f:
        for chrom, (length, bw) in chrom_data.items():
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            # Parse the header of the GFF
            head_line = "##sequence-region\t{} 1 {}".format(chrom, length)
            # Write on chromosome size files
            gff_head.append(head_line)
            chrom_f.write("{}\t{}\n".format(chrom, length))
            attrs = "ID={};bandwidth={}".format(chrom, bw)
            gff_item = [chrom, "Annotation", "chromosome", str(1), str(length),
                        ".", "+", ".",attrs]
            # Add chromosome to gff list
            gff_lst.append("\t".join(gff_item))
            # Parse genes of each chromosome
            chrom_genes = genes.query("chromosome == @chrom and species == @sp")
            col_order = ["chromosome","source","kind","start","end",
                         "dummy1","strand","dummy2","attrs"]
            chrom_genes = chrom_genes[col_order].set_index("chromosome")
            for gene in chrom_genes.to_records():
                gff_gene = "\t".join(gene)
                gff_lst.append(gff_gene)
        sp_gff = gff_head + gff_lst
        gff_f.write("\n".join(sp_gff))
