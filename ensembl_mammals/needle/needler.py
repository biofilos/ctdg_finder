import os
import sys

from Bio import SeqIO


"""
this script runs the Needleman-Wunsch algorithm for each of the sequence in a fasta file
against all the other sequences in the file. The only parameter is the input fasta file
"""
seqs_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'fasta'))

for seq1 in seqs_dict:
    subject_file = "{}_subject.fas".format(seq1.split("|")[2])
    needle_out = "_".join(seq1.split("|")[:3]) + '.needle'
    query_file = '{}_query.fas'.format(seq1.split("|")[2])
    query = open(query_file, 'w')
    seq1_record = seqs_dict[seq1]
    # seq1_record.name = seq1_record.name.replace('|')
    SeqIO.write(seqs_dict[seq1], query, 'fasta')
    query.close()
    if os.path.exists(subject_file):
        os.remove(subject_file)

    subject = open(subject_file, 'a')
    for seq2 in seqs_dict:
        if seq1 != seq2:
            SeqIO.write(seqs_dict[seq2], subject, 'fasta')
    subject.close()
    # Set command line for needleall
    cmd = 'needleall -sformat1 pearson -asequence {} '.format(query_file) + \
          '-sformat2 pearson -bsequence {} -gapopen 10 '.format(subject_file) + \
          '-gapextend 0.5 -aformat srspair -outfile {}'.format(needle_out)
    os.system(cmd)
    os.remove(query_file)
    os.remove(subject_file)
os.rename(sys.argv[1], 'done/{}'.format(sys.argv[1]))
