"""
Calculates the bandwidth parameter in a GFF and include that
information as attributes in the chromosome tracks
"""
import pybedtools as bt
from sys import argv
from collections import defaultdict
from HTSeq import GFF_Reader
import numpy as np


TOP_LEVEL = "chromosome"
gff = argv[1]
genome = argv[2]
gff_out = argv[3]

chroms = []
genes = defaultdict(list)
# Parse chromosomes and genes
for feat in GFF_Reader(gff):
    if feat.type == TOP_LEVEL:
        chroms.append(feat)
    else:
        genes[feat.iv.chrom].append(feat)

ig_d = {}
for chrom, c_genes in genes.items():
    # If there is only one gene, the bandwidth is zero
    n_genes = len(c_genes)
    if n_genes == 1:
        bw = 1
    # If there are two genes, the bandwidth is
    # the distance between them
    elif n_genes == 2:
        bw = c_genes[1].iv.start - c_genes[0].iv.end
        print(1)
    else:
        intergenics = [0] * (n_genes - 1)
        for ix1 in range(n_genes - 1):
            gene1 = c_genes[ix1].iv
            for ix2 in range(ix1 + 1, n_genes):
                gene2 = c_genes[ix2].iv
                intergenic = gene2.start - gene1.end
                intergenics[ix1] = intergenic
        bw = int(np.round(np.mean(intergenics) + np.std(intergenics), 0))

    ig_d[chrom] = bw


# Include the bandwidth parameter in the GFF and save it
with open(gff_out, "w") as f:
    for feat in GFF_Reader(gff):
        if feat.type == "chromosome":
            feat.attr["bandwidth"] = ig_d[feat.iv.chrom]
        f.write(feat.get_gff_line().replace('"', "").replace("; ", ";").replace(" ", "="))
print(1)