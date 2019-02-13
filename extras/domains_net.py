import networkx as nx
import pandas as pd
import json
import sys
from concurrent import futures
import os


"""
Transform a json dictionary of the form {dom1:[gene1,gene2]} into
a network where an edge between two genes exist if they share one domain (or
gene family)
"""


domain_dict = json.load(open(sys.argv[1]))
genes_df = pd.read_csv(sys.argv[2])
genes_df.loc[:, "chromosome"] = genes_df["chromosome"].astype(str)
out_dir = sys.argv[3]
cpus = int(sys.argv[4])
def species_graph(sp_chrom):
    sp, chrom = sp_chrom
    if not os.path.exists("{}/{}_{}.graphml".format(out_dir, sp, chrom)):
        print("{} ({})".format(sp, chrom))
        dom_G = nx.Graph()
        sp_genes_df = genes_df.loc[(genes_df["species"] == sp) &
                                (genes_df["chromosome"] == chrom)]
        sp_genes_df.sort_values("start", inplace=True)
        sp_genes = sp_genes_df["acc"].values
        for domain in domain_dict:
            genes_with_dom = list({x for x in domain_dict[domain] if x in sp_genes})
            for ix1, gene1 in enumerate(genes_with_dom):
                # dom_G.add_edge(gene1, gene1)
                for gene2 in genes_with_dom[ix1 + 1:]:
                    dom_G.add_edge(gene1, gene2)
        nx.write_graphml(dom_G, "{}/{}_{}.graphml".format(out_dir, sp, chrom))

genes_df.reset_index(inplace=True)
sp_chroms = genes_df.set_index(["species", "chromosome"]).index.unique()
with futures.ProcessPoolExecutor(cpus) as pool:
    res = pool.map(species_graph, sp_chroms)
_ = list(res)
