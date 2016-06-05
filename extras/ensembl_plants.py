import sys

import numpy as np
import os
import pandas as pd
import requests
from Bio import SeqIO

"""
This is the first step in the processing of the pipeline. In this step,
assembled chromosomes names and lengths
will be gathered from the Ensembl rest API
Run it as ./assembly_info.py species_list.txt
"""

species_list = []
server = "http://rest.ensemblgenomes.org"
for taxa in ['rosids', 'Chlamydomonas reinhardtii', 'Oryza sativa', 'Zea mays', 'Triticum aestivum']:
    ext = '/info/genomes/taxonomy/{}'.format(taxa)
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    for spp in r.json():
        species_list.append(spp['species'])


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
                                                              chrom['name'] in ['MT', 'cutchr', 'Mt', 'Pt'] or
                                                              'random' in chrom['name'] or 'hap' in chrom['name']):
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
    annotation_list = []
    # Parse a table for each chromosome
    for sp, chrom in species_table.index:
        print(sp, chrom)
        length = species_table.loc[(sp, chrom), 'length']
        # Process chromosomes by regions of 4e6 nt in length
        for start in np.arange(0, length, 4e6):
            end = start + 4e6
            if end > length:
                end = length
            # Get annotations for region
            ext = "/overlap/region/{}/{}:{:.0f}-{:.0f}?feature=gene".format(sp, chrom, start, end)
            r = requests.get(server + ext, headers={"Content-Type": "application/json"})
            # Fail gracefully if necessary
            if not r.ok:
                r.raise_for_status()
                sys.exit()
            # Format annotation as dataframe
            region_table = pd.DataFrame.from_dict(r.json())
            region_table['species'] = sp
            region_table['chromosome'] = chrom
            annotation_list.append(region_table)
    annotation = pd.concat(annotation_list)
    return annotation


def remove_overlaps(chrom_only):
    """
    Removes overlapping genes according to these rules:
    If the second gene (sorted by start coordinate) starts before the gene before ends
    AND both genes are in the same strand
    :param chrom_only: dataframe for each chromosome on a species
    :return: Dataframe
    """
    # Remove duplicated accession numbers
    chrom_only.reset_index(inplace=True)
    chrom_only.drop_duplicates('gene_id', inplace=True)
    # Set gne ids as indices
    chrom_only.set_index('gene_id', inplace=True)
    # Sort table by start coordinates
    chrom_only.sort_values('start', inplace=True)
    # Table from the first record to the une before the last record
    first_to_one_to_last = chrom_only.index.values[:-1]
    # Table from the second record to the last record
    second_to_last = chrom_only.index.values[1:]
    # If a record x has a start coordinate smaller than a record x - 1, both are overlapping
    overlapping = chrom_only.loc[second_to_last,
                                 'start'] < chrom_only.loc[first_to_one_to_last, 'end']
    # Consecutive genes in the same strand
    same_strand = chrom_only.loc[second_to_last, 'strand'] == chrom_only.loc[first_to_one_to_last, 'strand']
    # Flag for removal overlapping genes in the same strand
    for_removal = same_strand & overlapping
    for_removal_indices = for_removal.loc[for_removal].index
    # Remove overlapping genes in the same strand
    non_overlapping_table = chrom_only.drop(for_removal_indices)
    return non_overlapping_table


def download_seq(ensembls):
    """
    Downloads transcript sequences given a transcript ID
    :param ensembls: pandas dataframe with genes annotation
    :return:
    """
    plant_fasta = "plants/all_seqs.fa"
    all_seqs = open(plant_fasta, "a")
    #timer_to_flush = 0
    for gene in ensembls.index:
        seq_record = ensembls.loc[gene]
        sp = seq_record['species']
        chrom = str(seq_record['chromosome'])
        # seq = gene
        start = int(seq_record['start'])
        end = int(seq_record['end'])
        strand = int(seq_record['strand'])
        symbol = seq_record['external_name'].replace(' ','--')
        if symbol is np.nan:
            symbol = gene
        # Check if the gene has already been downloaded
        seq_name = "{}|{}|{}|{}|{}|{}|{}".format(sp, chrom, gene, symbol,
                                                 start, end, strand)
        seq_name = seq_name.replace(' ', '--')
        #timer_to_flush += 1
        ext = "/sequence/id/{}?multiple_sequences=1;type=protein".format(gene)
        headers = {"Content-Type": "application/json"}
        r = requests.get(server + ext, headers=headers)

        if not r.ok:
            r.raise_for_status()
            sys.exit()
        decoded = r.json()
        longest = {'seq': ''}
        if len(decoded) > 1:
            for isoform in decoded:
                if len(isoform['seq']) > len(longest['seq']):
                    longest = isoform
        else:
            longest = decoded[0]

        all_seqs.write(">{}\n{}\n".format(seq_name, longest['seq']))
        # Write in file every 10 genes
        #if timer_to_flush == 200:l_scaffold_0006_259
        #    all_seqs.flush()
        #    timer_to_flush = 0
        print(sp, symbol)
    all_seqs.close()


print("Checking species annotation")
good_sp = check_available(species_list)

print("Getting assembly information")
if not os.path.exists('plants/chromosomes.csv'):
    assemblies = pd.concat(map(get_assembly, good_sp))
    # Save to file
    assemblies.to_csv('plants/chromosomes.csv')
else:
    assemblies = pd.read_csv('plants/chromosomes.csv')

print("Downloading annotations")
if not os.path.exists('plants/annotation.csv'):
    annotation_table = get_annotation(assemblies)
    annotation_table.to_csv('plants/annotation.csv')
else:
    annotation_table = pd.read_csv('plants/annotation.csv')

protein_coding = annotation_table['biotype'] == 'protein_coding'
coding_annotation = annotation_table.loc[protein_coding]

# Remove overlapping genes
non_overlap_list = []
reindex = coding_annotation.set_index(['species', 'chromosome'])
for sp, chromosome in set(reindex.index):
    chrom_tab = coding_annotation.loc[(coding_annotation['species'] == sp) &
                                      (coding_annotation['chromosome'] == chromosome)]
    non_overlap_list.append(remove_overlaps(chrom_tab))
coding_no_overlap = pd.concat(non_overlap_list)
coding_no_overlap.to_csv("plants/coding_annotation.csv")
if os.path.exists("plants/all_seqs.fa"):
    downloaded = [gene.id.split("|")[2] for gene in SeqIO.parse("plants/all_seqs.fa", "fasta")]
else:
    downloaded = []
to_download = coding_no_overlap.loc[~(coding_no_overlap.index.isin(downloaded))]
download_seq(to_download)
