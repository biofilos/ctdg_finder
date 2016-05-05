import cgpFinder as cgp
import sys
'''
This script runs the blast step of CGPFinder.
Importantly, the sequence must be named Species_name.fa
In that way, the species name can be matched to the table.
The reason for this filter (the [sp] argument in cgp.blast.parse)
is that for the sampling process, only the blast hits for the same species
are needed
'''
input_file = sys.argv[1]
all_seqs = sys.argv[2]
cpus = sys.argv[3]
output_file = input_file.replace('fasta', 'blast')
sp = input_file.replace('.fasta', '').split('/')[-1]

cgp.blast.exe(cpus, input_file, all_seqs, output_file)

sub_table = cgp.blast.parse(output_file, 2, True, [sp], True)

