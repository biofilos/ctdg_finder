import pandas as pd
from sys import argv
from concurrent import futures

"""
Given a database directory, calculate the bandwidth parameter for
MeanShift and include it in the chromosomes.csv file
"""

db_dir = argv[1]


def intergenic(c_table):
    # If there is only one gene in the chromosome
    # the intergenic distance will be the end
    # coordinate of the gene
    c_table.sort_values("start", inplace=True)
    if c_table.shape[0] == 0:
        bandwidth = 0
    elif c_table.shape[0] == 1:
        bandwidth = c_table["end"]
    # If there are two, calculate it
    elif c_table.shape[0] == 2:
        bandwidth = abs(c_table.iloc[1]["start"] - c_table.iloc[0]["end"])
    else:
        pre_distances = c_table.iloc[1:]["start"].values - c_table.iloc[:-1]["end"].values
        # Set intergenic distances of overlapping genes to 0
        pre_distances[pre_distances < 0] = 0
        gene_distances = pre_distances
        if len(gene_distances) == 0:
            bandwidth = c_table.sort_values("end")["end"]
        elif len(gene_distances) == 1:
            bandwidth = gene_distances[0]
        else:
            bandwidth = gene_distances.mean() + gene_distances.std()
    return pd.np.round(float(bandwidth), 4)


def get_bw(ch_ix):
    chr_df = chromosomes.loc[ch_ix]
    sp = chr_df["species"]
    chrom = chr_df["chromosome"]
    length = int(chr_df["length"])
    sp_ch_genes = genes.loc[(genes["species"] == sp) &
                            (genes["chromosome"] == chrom)]
    bw = intergenic(sp_ch_genes)
    str_out = "{},{},{},{}\n".format(sp, chrom, length, bw)
    return str_out


genes = pd.read_csv(db_dir + "/genes_parsed.csv")
chromosomes = pd.read_csv(db_dir + "/pre_chromosomes.csv")
genes.sort_values(["species", "chromosome", "start"], inplace=True)

# Calculate bandwidth for each chromosome
with futures.ProcessPoolExecutor(8) as pool:
    bw_map = pool.map(get_bw,  chromosomes.index)

# Set file
chrom_file = open(db_dir + "/chromosomes.csv", "w")
chrom_file.write("species,chromosome,length,bandwidth\n")
for chrom in bw_map:
    chrom_file.write(chrom)
chrom_file.close()
