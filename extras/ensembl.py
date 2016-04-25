import sys
import os

import requests
import numpy as np
import pandas as pd

"""
This is the first step in the processing of the pipeline. In this step,
assembled chromosomes names and lengths
will be gathered from the Ensembl rest API
Run it as ./assembly_info.py species_list.txt
"""

species_list = [sp.lower().strip('\n') for sp in open(sys.argv[1])]

# Filter out species without annotations


def check_available(sp_list):
    """
    Check if there is assembly information for
    a species
    :param sp_list: list of species in the format species_name
    :return:
    """
    good_sp = []
    for sp in sp_list:
        server = "http://rest.ensembl.org"
        ext = "/info/assembly/{}?".format(sp)
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
        if r.status_code != 200:
            print("{} is not in the database".format(sp))
        else:
            good_sp.append(sp)
    return good_sp


def get_assembly(sp):
    """
    This function downloads assembly information form Ensembl, and filter out scaffolds and
    non assembled chromosomes
    :param sp: species list in a text file
    :return: Dataframe
    """
    # Set server location
    server = "http://rest.ensembl.org"
    ext = "/info/assembly/{}?".format(sp)
    # Query information
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    try:
        # Check for transaction integrity (?)
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        # Convert JSON format to dictionary
        decoded = r.json()
        print(sp)
        # Filter out non-assembled chromosomes
        assembled_chroms = {}
        for chrom in decoded['top_level_region']:
            # Some unassembled chromosomes in horse are named UnXXX, so an explicit statement is set
            # to filter those out. Also, skip mitochondrial chromosomes
            if chrom['coord_system'] == 'chromosome' and not (chrom['name'].startswith('Un') or
                                                              chrom['name'].startswith('un') or
                                                              chrom['name'] in ['MT', 'cutchr'] or
                                                              'random' in chrom['name'] or
                                                              'hap' in chrom['name']):
                assembled_chroms[chrom['name']] = [sp, chrom['length']]
        # Organize results in data frame
        # Process assembly if there were assemblied chromosomes
        if len(assembled_chroms) == 0:
            return None
        else:
            assemblies_table = pd.DataFrame.from_dict(assembled_chroms, orient='index')
            assemblies_table.reset_index(inplace=True)
            assemblies_table.columns = ['chromosome', 'species', 'length']
        return assemblies_table
    except requests.HTTPError:
        print()


def get_annotation(assemblies_df):
    """
    Extract annotations for all genes in the species under study
    :param assemblies_df: dataframe of species and chromosomes
    :return: dataframe
    """
    # Get species information
    assemblies_df.chromosome = assemblies_df.chromosome.astype(str)
    assemblies_df.length = assemblies_df.length.astype(int)
    species_table = assemblies_df.set_index(['species', 'chromosome'])
    # Set server address
    server = "http://rest.ensembl.org"
    annotation_list = []
    for sp, chrom in species_table.index:
        print(sp, chrom)
        length = species_table.loc[(sp, chrom), 'length']
        for start in np.arange(1, length, 4e6):
            end = start + 4e6
            if end > length:
                end = length
            ext = "/overlap/region/{}/{}:{:.0f}-{:.0f}?feature=gene".format(sp, chrom, start, end)
            r = requests.get(server + ext, headers={"Content-Type": "application/json"})
            if not r.ok:
                r.raise_for_status()
                sys.exit()
            region_table = pd.DataFrame.from_dict(r.json())
            region_table['species'] = sp
            annotation_list.append(region_table)
    annotation = pd.concat(annotation_list)
    return annotation


def get_transcript_id(gene_id):
    """
    Extract transcript ID from gene ID
    :param gene_id: gene ID
    :return: transcript ID
    """
    transcript = ''
    # Open log file
    out_file = 'ensembl/geneid_transcript.csv'
    file_o = open(out_file, 'a')
    #Check if out file exists, and load it if necessary
    if not os.stat(out_file).st_size == 0:
        transcript_csv = pd.read_csv(out_file, header=None)
        transcript_csv.set_index(0, inplace=True)
    else:
        transcript_csv = pd.DataFrame(columns=[0, 1])
    if gene_id in transcript_csv.index:
        transcript = transcript_csv.loc[gene_id, 1]
    else:
        # Set server
        server = "http://rest.ensembl.org"
        ext = "/overlap/id/{}?feature=transcript".format(gene_id)
        # wait a little bit
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
        # Check transaction
        if not r.ok:
            print("wrong")
            r.raise_for_status()
            sys.exit()
        # Filter out transcripts that are not directly associated with the gene accession
        isoforms = [x for x in r.json() if x['Parent'] == gene_id]
        # Sort transcript by length
        isoforms_lengths = {i['id']: i['end'] - i['start'] for i in isoforms}
        # Extract the longest isoform
        transcript = sorted(isoforms_lengths, key=isoforms_lengths.get)[-1]
        file_o.write("{},{}\n".format(gene_id, transcript))
        file_o.close()
    return transcript


def download_seq(transcript):
    """
    Downloads transcript sequences given a transcript ID
    :param transcript: python list of transcripts ID
    :return:
    """
    server = "http://rest.ensemblgenomes.org"
    ext = "/sequence/id/{}?object_type=transcript;type=cds".format(transcript)
    print(ext)
    r = requests.get(server + ext, headers={"Content-Type": "text/x-fasta"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    print("done")
    return decoded

print("Checking species annotation")
good_sp = check_available(species_list)

print("Getting assembly information")
if not os.path.exists('ensembl/chromosomes.csv'):
    assemblies = pd.concat(map(get_assembly, good_sp))
    # Save to file
    assemblies.to_csv('ensembl/chromosomes.csv')
else:
    assemblies = pd.read_csv('ensembl/chromosomes.csv')

print("Downloading annotations")
if not os.path.exists('ensembl/annotation.csv'):
    annotation_table = get_annotation(assemblies)
    annotation_table.to_csv('ensembl/annotation.csv')
else:
    annotation_table = pd.read_csv('ensembl/annotation.csv')

print("Downloading transcript IDs")
if 'transcript' in annotation_table.columns:
    pass
else:
    protein_coding = annotation_table['biotype'] == 'protein_coding'
    annotation_table.loc[protein_coding,
                         'transcript'] = annotation_table.loc[protein_coding, 'gene_id'].map(get_transcript_id)
    annotation_table.to_csv("ensembl/annotation.csv")

coding_annotation = annotation_table.loc[annotation_table.biotype == 'protein_coding']
coding_annotation.to_csv("ensembl/coding_annotation.csv")

_ = download_seq(coding_annotation.transcript.values[0])
