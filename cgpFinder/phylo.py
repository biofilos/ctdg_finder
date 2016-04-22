from Bio import Entrez, SeqIO
import os
import numpy as np
from cgpFinder.annotation import name_family
Entrez.email = 'juan.f.ortiz@vanderbilt.edu'

def download_genes(ordered):
    """Download protein and dna sequences
    """
    # If sequence file is empty or does not exist
    if not os.path.exists('{0}/genes/genes.fa'.format(name_family)) or os.stat(
            '{0}/genes/genes.fa'.format(name_family)).st_size == 0:
        # Initialize a list of already downloaded sequences
        dones = []
    else:
        # If some of the sequences were already downloaded, extract their names, so that they are not
        # downloaded again
        dones = [secs.id.split("|")[1] for secs in SeqIO.parse('{0}/genes/genes.fa'.format(name_family), 'fasta')]
        print("{} genes already downloaded".format(len(dones)))

    # Set file names
    dna_file = open('{0}/genes/genes_rna.fa'.format(name_family), 'a')
    aa_file = open('{0}/genes/genes.fa'.format(name_family), 'a')
    # If the hmmer table is not indexed by the accession number, set it
    if not ordered.index.name == 'prot_acc':
        ordered.set_index('prot_acc', inplace=True)
    # Download all the sequences that are not in the "dones" list
    for aa_acc in set(ordered.index) - set(dones):
        # Set fields
        parsing = ordered.loc[aa_acc]
        species = parsing['species'].replace(' ', '_')
        symbol = parsing['symbol']
        cluster = parsing['cluster']
        order = parsing['order']
        # Download protein record
        prot_gi = Entrez.read(Entrez.esearch(db='protein', term=aa_acc))['IdList'][0]
        record = SeqIO.read(Entrez.efetch(db='protein', id=str(prot_gi), rettype='gb'), 'gb')
        # Extract sequence from record
        prot = str(record.seq)
        dna = ''
        for feats in record.features:
            # Extract CDS annotation
            if feats.type == 'CDS' and dna == '':
                # Extract accession number of nucleotide record that codes for the CDS
                dna_id_coord = feats.qualifiers['coded_by'][0].split(':')
                dna_id = dna_id_coord[0]
                dna_coord = dna_id_coord[1].split('..')
                # Download mRNA record
                dna_gi = Entrez.read(Entrez.esearch(db='nucleotide', term=dna_id))['IdList'][0]
                record_dna = SeqIO.read(Entrez.efetch(db='nucleotide', id=dna_gi, rettype='gb'), 'gb')
                try:
                    # Extract mRNA sequence from record
                    dna = record_dna.seq[int(dna_coord[0]) - 1:int(dna_coord[1])]
                    # Make the coordinates uniform. For some reason, sometimes they are parsed as real
                    # and some others as integers
                    if type(order) == float:
                        # Set sequence names
                        aa_file.write(">{}|{}|{}|{}|{:.0f}\n{}\n".format(species, aa_acc, symbol, cluster, order, prot))
                        aa_file.flush()
                        dna_file.write(">{}|{}|{}|{}|{:.0f}\n{}\n".format(species, aa_acc, symbol, cluster, order, dna))
                        dna_file.flush()
                    else:
                        aa_file.write(">{}|{}|{}|{}|{}\n{}\n".format(species, aa_acc, symbol, cluster, order, prot))
                        aa_file.flush()
                        dna_file.write(">{}|{}|{}|{}|{}\n{}\n".format(species, aa_acc, symbol, cluster, order, dna))
                        dna_file.flush()
                except ValueError:
                    # ValueError in the coordinates indicate broken annotation, these sequence will be omitted
                    print("Incomplete annotation: {0} -> {1}".format(dna_id, dna_coord))
                    # print(e)
                # When the CDS is processed, don't keep reading annotation, go to the next gene
                break
    # Reset index of hmmer table
    ordered.reset_index(inplace=True)


def align():
    """
    Align the sequences and write convenient Slurm files to run raxml and modeltest
    :return:
    """
    # Align protein sequences using MAFFT
    run_mafft = 'mafft --anysymbol {0}/genes/genes.fa > {0}/genes/genes.aln'.format(name_family)
    print("Protein alignment for {}, calculated".format(name_family))
    os.system(run_mafft)
    # Back-translate protein alignment to codons using Pal2Nal and save the results in fasta and paml formats
    run_pal2nal_paml = 'pal2nal.pl -output paml \
    {0}/genes/genes.aln {0}/genes/genes_rna.fa > {0}/paml/genes_codon.aln'.format(name_family)
    os.system(run_pal2nal_paml)
    run_pal2nal_fasta = 'pal2nal.pl -output fasta \
    {0}/genes/genes.aln {0}/genes/genes_rna.fa > {0}/raxml/genes_codon.aln'.format(
        name_family)
    os.system(run_pal2nal_fasta)
    # Write slurm script for modeltest
    modeltest_cmd = open('{0}/raxml/{0}_modeltest.slurm'.format(name_family), 'w')
    modeltest_text = '''#!/bin/bash

#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=7-00:00:00
#SBATCH --output=raxml_{0}.out
#SBATCH --mail-user={1}
#SBATCH --mail-type=ALL
#SBATCH --job-name=modeltest_{0}

export I_MPI_PMI_LIBRARY=/usr/scheduler/slurm/lib/libmpi.so

cd /data/ortizjf/placenta/{0}/modeltest
date
echo $SLURM_JOB_NODELIST

java -jar /home/ortizjf/bin/jmodeltest-2.1.5/jModelTest.jar -d genes_codon.aln -o {0}.modeltest -f -i -g 4 -AIC -tr 8
date'''.format(name_family, Entrez.email)
    modeltest_cmd.write(modeltest_text)
    modeltest_cmd.close()
    # Generate random numbers for seeds in RAxML
    raxml_random1 = np.random.randint(10000)
    raxml_random2 = np.random.randint(10000)
    # Write slurm script for RAxML
    raxml = open('{0}/raxml/{0}_raxml.slurm'.format(name_family), 'w')
    raxml_text = '''#!/bin/bash

#SBATCH --nodes=10
#SBATCH --tasks-per-node=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=7-00:00:00
#SBATCH --output=raxml_{0}.out
#SBATCH --mail-user={1}
#SBATCH --mail-type=ALL
#SBATCH --job-name=RAxML_{0}

setpkgs -a mpiexec
#export I_MPI_PMI_LIBRARY=/usr/scheduler/slurm/lib/libmpi.so


cd /data/ortizjf/placenta/{0}/raxml

date
echo $SLURM_JOB_NODELIST
mpirun -n 80 raxmlHPC-MPI-SSE3 -f a -s {0}_codon.phy -n {0} -m GTRGAMMA -p {2} -x {3} -N 1000
date'''.format(name_family, Entrez.email, raxml_random1, raxml_random2)
    raxml.write(raxml_text)
    raxml.close()
    print('Codon alignment for {}, generated'.format(name_family))

