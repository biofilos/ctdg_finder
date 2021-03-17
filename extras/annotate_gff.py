"""
Calculates the bandwidth parameter in a GFF and include that
information as attributes in the chromosome tracks.
Also, annotate gene families in gff
"""
import pybedtools as bt
from sys import argv
from collections import defaultdict
from HTSeq import GFF_Reader
import numpy as np


def parse_hmmer(pfam_in):
    # Initialize pfam dictionary
    domains_db = {}

    for line in open(pfam_in):
        # Skip comment lines
        if not line.startswith("#"):
            # Process each line
            data_line = [x for x in line.split(' ') if x != ""]
            acc = data_line[0]
            seq = acc
            pfam = data_line[3]
            # Include entries with a minimum coverage
            # in both the hmm and the query of 30%
            threshold = 0.3
            hmm_start = int(data_line[15])
            hmm_end = int(data_line[16])
            hmm_aln_len = hmm_end - hmm_start
            hmm_len = int(data_line[5])
            hmm_coverage = hmm_aln_len / hmm_len

            q_len = int(data_line[2])
            q_end = int(data_line[20])
            q_start = int(data_line[19])
            q_aln_len = q_end - q_start
            query_coverage = q_aln_len / q_len
            if query_coverage >= threshold and hmm_coverage >= threshold:
                if seq not in domains_db:
                    domains_db[seq] = [pfam]
                else:
                    domains_db[seq].append(pfam)

    return domains_db


TOP_LEVEL = argv[1]
gff = argv[2]
hmmer_dom_in = argv[3]
gff_out = argv[4]
chroms_out = argv[5]

domains_db = parse_hmmer(hmmer_dom_in)

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
    c_genes = sorted(c_genes, key=lambda x: x.iv.start)
    # If there is only one gene, the bandwidth is zero
    n_genes = len(c_genes)
    if n_genes == 1:
        bw = 1
    # If there are two genes, the bandwidth is
    # the distance between them
    elif n_genes == 2:
        bw = c_genes[1].iv.start - c_genes[0].iv.end
        if bw <= 0:
            bw = 1
    else:

        intergenics = [0] * (n_genes - 1)
        for ix1 in range(n_genes - 1):
            gene1 = c_genes[ix1].iv
            ix2 = ix1 + 1
            gene2 = c_genes[ix2].iv
            # In case of overlapping genes, take the absolute value
            intergenic = abs(gene2.start - gene1.end)
            intergenics[ix1] = intergenic
        bw = int(np.round(np.mean(intergenics) + np.std(intergenics), 0))

    ig_d[chrom] = bw


# Include the bandwidth parameter in the GFF and save it
with open(gff_out, "w") as f, open(chroms_out, "w") as c_out:
    for feat in GFF_Reader(gff):
        if feat.type == "chromosome":
            feat.attr["bandwidth"] = ig_d[feat.iv.chrom]
            c_out.write(f"{feat.iv.chrom}\t{feat.iv.end}\n")
        elif feat.name in domains_db:
            families = set(domains_db[feat.name])
            feat.attr["families"] = ",".join(families)
        else:
            feat.attr["families"] = "None"
        f.write(feat.get_gff_line().replace('"', "").replace("; ", ";").replace(" ", "="))
print(gff_out)
