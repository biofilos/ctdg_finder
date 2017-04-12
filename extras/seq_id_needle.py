import os
import sys
from concurrent import futures
from glob import glob

import pandas as pd
from Bio import SeqIO
import shutil

# Load fasta file
in_file = sys.argv[1]
# Out file
out_file = sys.argv[2]
# annotation file
gene_file = sys.argv[3]

# annotation file output
gene_selected_out = sys.argv[4]

# Percentage similarity threshold
perc_id = float(sys.argv[5])

# CPUS
cpus = int(sys.argv[6])

def clean_seqname(seq_path):
    return seq_path.replace(".fa", "").replace("seqs/", "")


def get_id(hit):
    """
    Calculate identity and similarity between a pair of sequences
    :param hit: list of 2 sequence filenames
    :return:list: [sequence1, sequence2, identity, similarity]
    """
    seq1, seq2 = hit
    seq1_name = clean_seqname(seq1)
    seq2_name = clean_seqname(seq2)
    out = "{}_{}.needle".format(seq1.replace(".fa", ""), seq2.replace(".fa", "").replace("seqs/", ""))
    os.system(
        "needle -asequence {} -bsequence {} -gapopen 10 -gapextend 0.5 -outfile {}".format(seq1,
                                                                                           seq2,
                                                                                           out))
    for line in open(out):
        if "Identity" in line:
            iden = float(line.split("(")[1].split("%")[0].strip(" "))
        elif "Similarity" in line:
            simil = float(line.split("(")[1].split("%")[0].strip(" "))
            break
    os.remove(out)
    return [seq1_name, seq2_name, iden, simil]

# def hit_name(row):


# create sequence directory
if os.path.exists("seqs"):
    shutil.rmtree("seqs")
os.makedirs("seqs")
for seq in SeqIO.parse(in_file, "fasta"):
    #if not seq.name.split("|")[-1] in ["0", "na_95", "na_ms"]:
    seq_name = "seqs/" + seq.name.split("|")[2] + ".fa"
    SeqIO.write(seq, seq_name, "fasta")

# Get all pairs of sequences
seq_pairs = []
for seq1 in glob("seqs/*.fa"):
    for seq2 in glob("seqs/*.fa"):
        if seq1 != seq2:
            seq_pairs.append([seq1, seq2])

with futures.ProcessPoolExecutor(cpus) as pool:
    id_calculated = pool.map(get_id, seq_pairs)
# Load similarity data as a dataframe
id_data = []
for datum in id_calculated:
    if datum[3] > perc_id:
        id_data.append(datum)
id_hits = pd.DataFrame(id_data, columns=["seq1", "seq2", "identity", "similarity"])

# Read gene annotation
genes = pd.read_csv(gene_file)
# Extract gene annotation for sequences in question
genes_selected = genes.loc[genes["acc"].isin(set(id_hits["seq1"]))]
genes_selected.to_csv(gene_selected_out)

# Save hits (directed graph)
id_hits.set_index("seq1", inplace=True)
id_hits = id_hits.loc[id_hits["similarity"] >= perc_id]
id_hits.to_csv(out_file + "_{}directed".format(perc_id))

# Use only one identity-similarity (undirected graph)
id_dict = {}
records = id_hits.to_records()
for rec in records:
    hit_name = "XXX".join(sorted([rec[0], rec[1]]))
    if hit_name not in id_dict:
        id_dict[hit_name] = [rec[2], rec[3]]
    else:
        id_dict[hit_name][0] = pd.np.mean([id_dict[hit_name][0], rec[2]])
        id_dict[hit_name][1] = pd.np.mean([id_dict[hit_name][1], rec[3]])

# Format undirected graph
undirected_list = []
for hit in id_dict:
    seq1, seq2 = hit.split("XXX")
    identity, similarity = id_dict[hit]
    undirected_list.append([seq1, seq2, identity, similarity])
# Save undirected graph data
undirected_df = pd.DataFrame(undirected_list, columns=["seq1", "seq2", "identity", "similarity"])
undirected_df.set_index("seq1", inplace=True)
undirected_df.to_csv(out_file + "_{}undirected".format(perc_id))

###
# Export a fasta file for each component
###
annotation = genes
annotation.set_index("acc", inplace=True)
net_df = undirected_df
net_df.set_index(["seq1", "seq2"], inplace=True)

# Initialize graph
G = nx.Graph()
# Load edges in graph
G.add_edges_from(net_df.index)
# Annotate edges
for edge in G.edges_iter():
    for edge_info in net_df.columns:
        # In case the nodes of an edge are in different order, flip them
        try:
            G.edge[edge[0]][edge[1]][edge_info] = net_df.loc[edge, edge_info]
        except KeyError:
            G.edge[edge[0]][edge[1]][edge_info] = net_df.loc[(edge[1], edge[0]), edge_info]

cc_ix = 1
ccs = {}
un_cl_set, cl_set = set(), set()
cc_df = pd.DataFrame(columns=["name", "nodes", "clustered", "unclustered"])
cc_df.set_index("name", inplace=True)

rule = {str:lambda cluster: "na" in cluster or cluster == "0",
       pd.Series: lambda cluster:("na" in cluster or cluster == "0").all()}
for cc in nx.connected_component_subgraphs(G):
    cc_name = "cc_" + str(cc_ix)
    cc_ix += 1
    nodes = len(cc)
    cl, un_cl = 0, 0
    for cc_node in cc.nodes():
        cluster = cc.node[cc_node]["cluster"]
        if rule[type(cluster)](cluster):
            un_cl += 1
            un_cl_set.add(cc_node)
        else:
            cl += 1
            cl_set.add(cc_node)
    cc_df.loc[cc_name] = [nodes, cl, un_cl]
    ccs[cc_name] = cc

cc_df.loc[:, "prop_cl"] = cc_df["clustered"] / cc_df["nodes"]
seq_dict = SeqIO.to_dict(SeqIO.parse(in_file, "fasta"), key_function=lambda x: x.name.split("|")[2])
# Save sequences from each component in a different file
basename = in_file.split(".")[0]
for cc in ccs:
    fileO = open("{}_{}.fa".format(basename, cc), "w")
    for node in ccs[cc].nodes():
        SeqIO.write(seq_dict[node], fileO, "fasta")
    fileO.close()
