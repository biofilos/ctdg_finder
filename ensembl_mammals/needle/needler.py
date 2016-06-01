import os
import sys

from Bio import SeqIO

seqs_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'fasta'))

for seq1 in seqs_dict:
    query = open('query.fa', 'w')
    seq1_record = seqs_dict[seq1]
    # seq1_record.name = seq1_record.name.replace('|')
    SeqIO.write(seqs_dict[seq1], query, 'fasta')
    query.close()
    if os.path.exists('subject.fa'):
        os.remove('subject.fa')

    subject = open('subject.fa', 'a')
    for seq2 in seqs_dict:
        if seq1 != seq2:
            SeqIO.write(seqs_dict[seq2], subject, 'fasta')
    subject.close()
    # Set command line for needleall
    cmd = 'needleall -sformat1 pearson -asequence query.fa -sformat2 pearson -bsequence subject.fa -gapopen 10 ' \
          '-gapextend 0.5 -aformat srspair -outfile {}.needle'.format(seq1.split("|")[2])
    os.system(cmd)
