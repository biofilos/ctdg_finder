import sys
import json
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



# Check in which database each species is
file_in = sys.argv[2]
# Load assemblies dictionary
assemblies = sys.argv[1]
out_dir = file_in + '_db'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
ensembls = json.loads(open(assemblies).read())

## This part of the code was implemented as a separate script
#for server in ["http://rest.ensemblgenomes.org", "http://rest.ensembl.org"]:
#    for ext in ["/info/species?", "/info/species?division=EnsemblPlants"]:
#        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
#        if not r.ok:
#            r.raise_for_status()
#            sys.exit()
#        else:
#            decoded = r.json()
#            for sp in decoded['species']:
#                species = sp['name']
#                ensembls[species] = server

# Get the server for each of the species of interest
pre_list = [line.strip('\n') for line in open(file_in) if line != '\n']

sp_server = {}
for sp in pre_list:
    sp_server[sp] = ensembls[sp]

# Filter out species without annotations


def get_assembly(sp_item):
    """
    This function downloads assembly information form Ensembl, and filter out scaffolds and
    non assembled chromosomes
    :param sp_item: species list in a text file
    :return: Dataframe
    """
    # Set server location
    server = ensembls[sp_item]
    ext = "/info/assembly/{}?".format(sp_item)
    # Query information
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    try:
        # Check for transaction integrity (?)
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        # Convert JSON format to dictionary
        decoded = r.json()
        # Filter out non-assembled chromosomes
        assembled_chroms = {}
        for chrom in decoded['top_level_region']:
            # Some unassembled chromosomes in horse are named UnXXX, so an explicit statement is set
            # to filter those out. Also, skip mitochondrial chromosomes
            if chrom['coord_system'] == 'chromosome' and not (chrom['name'].startswith('Un') or
                                                              chrom['name'].startswith('un') or
                                                              'random' in chrom['name'] or
                                                              'hap' in chrom['name'] or
                                                              chrom['name'] in ['MT', 'cutchr', 'Mt', 'Pt',
                                                                                'random', 'hap',
                                                                                'dmel_mitochondrion_genome',
                                                                                'Mito', 'MtDNA']):
                assembled_chroms[chrom['name']] = [sp_item, chrom['length']]
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
        server = sp_server[sp]
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


print("Getting assembly information")
if not os.path.exists('{}/chromosomes.csv'.format(out_dir)):
    assemblies = pd.concat(map(get_assembly, sp_server))
    # Save to file
    assemblies.to_csv('{}/chromosomes.csv'.format(out_dir))
else:
    assemblies = pd.read_csv('{}/chromosomes.csv'.format(out_dir))

print("Downloading annotations")
if not os.path.exists('{}/annotation.csv'.format(out_dir)):
    annotation_table = get_annotation(assemblies)
    annotation_table.to_csv('{}/annotation.csv'.format(out_dir))
else:
    annotation_table = pd.read_csv('{}/annotation.csv'.format(out_dir))

protein_coding = annotation_table['biotype'] == 'protein_coding'
coding_annotation = annotation_table.loc[protein_coding]

# Remove overlapping genes
non_overlap_list = []
reindex = coding_annotation.set_index(['species', 'chromosome'])
if not os.path.exists("{}/coding_no_overlap.csv".format(out_dir)):
    for sp, chromosome in set(reindex.index):
        chrom_tab = coding_annotation.loc[(coding_annotation['species'] == sp) &
                                          (coding_annotation['chromosome'] == chromosome)]
        non_overlap_list.append(remove_overlaps(chrom_tab))
    coding_no_overlap = pd.concat(non_overlap_list)
    coding_no_overlap.to_csv('{}/coding_no_overlap.csv'.format(out_dir))
else:
    coding_no_overlap = pd.read_csv("{}/coding_no_overlap.csv".format(out_dir))
coding_no_overlap['length'] = coding_no_overlap['end'] - coding_no_overlap['start']
coding_no_overlap['length'] = coding_no_overlap['length'].astype(int)
coding_no_overlap.reset_index(inplace=True)
coding_new_names = coding_no_overlap.loc[:, ['gene_id', 'start', 'end', 'length',
                                              'strand', 'chromosome', 'species', 'external_name']]
coding_new_names.rename(columns={'gene_id':'acc','external_name':'symbol'},inplace=True)
null_symbols = coding_new_names['symbol'].isnull()
coding_new_names.loc[null_symbols, 'symbol'] = coding_new_names.loc[null_symbols, 'acc']
coding_new_names.set_index('acc', inplace=True)
coding_new_names.to_csv("{}/genes_parsed.csv".format(out_dir))

