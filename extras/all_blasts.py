import sys
from glob import glob

import os
from Bio import SeqIO
ctdg_path = sys.argv[0].replace("extras/all_blasts.py", "")

sys.path.append(ctdg_path)
from ctdg_finder import CtdgRun, CtdgConfig

'''
This script runs the blast step of CTDGFinder.
Importantly, the sequence must be named Species_name.fa
In that way, the species name can be matched to the table.
The reason for this filter (the [sp] argument in ctdg.blast.parse)
is that for the sampling process, only the blast hits for the same species
are needed
'''
# Specify directory where the blast database is located
all_seqs = sys.argv[1]
out_dir = sys.argv[2]
db = sys.argv[3]
cpus = sys.argv[4]

# Create directory if necessary
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
if not os.path.exists("{}/done".format(out_dir)):
    os.makedirs("{}/done".format(out_dir))

# Generate species-specific proteomes
for gene in SeqIO.parse(all_seqs, 'fasta'):
    sp_name = gene.name.split("|")[0]
    if not (os.path.exists("{}/{}.blast_out".format(out_dir, sp_name)) or
                os.path.exists("{}/{}.blast".format(out_dir, sp_name))):
        fasta_name = "{}/{}.fasta".format(out_dir, sp_name)
        fileO = open(fasta_name, 'a')
        SeqIO.write(gene, fileO, 'fasta')
        fileO.close()

# Run Blast step for each proteome
ctdg_config = CtdgConfig(out_dir=out_dir, db=db, blast_samples=1000, evalue=1, sp=[])
ctdg_config.check_arguments()
ctdg_config.init_tables()

for input_file in glob("{}/*.fasta".format(out_dir)):
    ctdg = CtdgRun(ctdg_config, name_family="all_blast", ref_sequence=input_file)
    output_file = input_file.replace('fasta', 'blast')
    sp = input_file.replace('.fasta', '').split('/')[-1]
    if not os.path.exists(output_file):
        ctdg.blast_exe(cpus, ctdg.ref_sequence, all_seqs, output_file)

    sub_table = ctdg.blast_parse(output_file, 2, [sp], True, True)
    for out in [input_file, output_file, output_file + "_filtered"]:
        os.rename(out, out.replace(out.split('/')[-1], "done/{}".format(out.split('/')[-1])))
