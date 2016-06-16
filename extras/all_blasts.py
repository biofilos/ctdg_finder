import cgpFinder as cgp
import sys
from Bio import SeqIO
from glob import glob
'''
This script runs the blast step of CGPFinder.
Importantly, the sequence must be named Species_name.fa
In that way, the species name can be matched to the table.
The reason for this filter (the [sp] argument in cgp.blast.parse)
is that for the sampling process, only the blast hits for the same species
are needed
'''
all_seqs = sys.argv[1]
out_dir = all_seqs.replace(all_seqs.split('/')[-1], 'blasts')
cpus = sys.argv[2]

# Generate species-specific proteomes
#for gene in SeqIO.parse(all_seqs,'fasta'):
#    sp_name = gene.name.split("|")[0]
#    fileO = open("{}/{}.fasta".format(out_dir,sp_name),'a')
#    SeqIO.write(gene,fileO,'fasta')
#    fileO.close()

# Run Blast step for each proteome
for input_file in glob("{}/*.fasta".format(out_dir)):
    output_file = input_file.replace('fasta', 'blast')
    sp = input_file.replace('.fasta', '').split('/')[-1]
    
    cgp.blast.exe(cpus, input_file, all_seqs, output_file)
    
    sub_table = cgp.blast.parse(output_file, 2, True, [sp], True)

