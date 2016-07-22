import sys
from glob import glob

import os
from Bio import SeqIO

sys.path.append("..")
from cgp_exec import CGP

'''
This script runs the blast step of CGPFinder.
Importantly, the sequence must be named Species_name.fa
In that way, the species name can be matched to the table.
The reason for this filter (the [sp] argument in cgp.blast.parse)
is that for the sampling process, only the blast hits for the same species
are needed
'''
# Specify directory where the blast database is located
all_seqs = sys.argv[1]
#set a directory "blasts" as output
out_dir = all_seqs.replace(all_seqs.split('/')[-1], 'blasts')
# Set number of cpus to be used
cpus = sys.argv[2]

# Create directory if necessary
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
# Generate species-specific proteomes
for gene in SeqIO.parse(all_seqs, 'fasta'):
    sp_name = gene.name.split("|")[0]
    fasta_name = "{}/{}.fasta".format(out_dir, sp_name)
    fileO = open(fasta_name, 'a')
    SeqIO.write(gene, fileO, 'fasta')
    fileO.close()

# Run Blast step for each proteome
for input_file in glob("{}/*.fasta".format(out_dir)):
    cgp = CGP(name_family="all_blast", evalue=1, out_dir=out_dir, db=all_seqs,
              blast_samples=0, sp=[], ref_sequence=input_file)
    output_file = input_file.replace('fasta', 'blast')
    sp = input_file.replace('.fasta', '').split('/')[-1]

    cgp.blast_exe(cpus, cgp.ref_sequence, all_seqs, output_file)
    
    sub_table = cgp.blast_parse(output_file, acc_col=2, tab=True, sp_list=cgp.sp, for_dict=True)

