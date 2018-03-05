import networkx as nx
import sys
from glob import glob
from concurrent import futures
import os

"""
python parse_matrices.py graphs_dir
"""
graphs = glob(sys.argv[1] + "/*.graphml")
out_matrices = (sys.argv[2])


def parse_matrix(graph):
    filename = "{}/{}".format(out_matrices, graph.split("/")[-1].replace(".graphml", ".csv"))
    if not os.path.exists(filename):
        G = nx.read_graphml(graph)
        df = nx.to_pandas_adjacency(G)
        df.to_csv(filename)
        print(filename)


with futures.ProcessPoolExecutor(4) as pool:
    _ = pool.map(parse_matrix, graphs)

d = list(_)
