from glob import glob
from Bio import SeqIO, pairwise2
import pandas as pd
import numpy as np
from concurrent import futures


# seq_raw = SeqIO.to_dict(SeqIO.parse("../cgpf_ncbi/all_seqs.fa", "fasta"))
# seq_dict = {x.split('|')[2]: seq_raw[x].seq for x in seq_raw if x.split("|")[0]=='Homo sapiens'}


def cluster_pid(folder):
    result = []
    f_name = folder.split("/")[-1]
    try:
        genes = pd.read_csv(folder + "/report/" + f_name + "_genes.csv")
        genes = genes.loc[~(genes['cluster'].isin(['na', '0', 0])) & (genes['species'] == 'Homo sapiens')]
        if genes.shape[0] > 0:
            for cluster in set(genes['cluster']):
                pids = []
                accs = genes.loc[genes['cluster'] == cluster, 'prot_acc'].values
                for seq1 in accs:
                    for seq2 in accs:
                        seq_1 = [x.seq for x in SeqIO.parse("../cgpf_ncbi/all_seqs.fa", 'fasta') if
                                 x.name.split("|")[2] == seq1]
                        seq_2 = [x.seq for x in SeqIO.parse("../cgpf_ncbi/all_seqs.fa", 'fasta') if
                                 x.name.split("|")[2] == seq2]
                        aln = pairwise2.align.globalxx(seq_1[0], seq_2[0])[0]
                        mean_len = (len(aln[0]) + len(aln[1])) / 2
                        pids.append(aln[2] / mean_len)

                n_genes = len(pids)
                mean_pid = np.mean(pids)
                sd_id = np.std(pids)
                result.append(cluster, n_genes, mean_pid, sd_id)
                print(cluster)
        return result
    except OSError:
        return None


with futures.ProcessPoolExecutor() as pool:
    pid_results = pool.map(cluster_pid, glob("../cgpf_ncbi/p_atlas/finish/*"))

with open("../cgpf_ncbi/p_atlas_pids.csv", "w") as file_o:
    file_o.write("cluster,n_genes,mean_pid,sd_pid\n")
    for line in pid_results:
        if line:
            file_o.write(line+'\n')
