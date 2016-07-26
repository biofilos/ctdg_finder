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



# Check in which database each species is
print("Checking databases")
ensembls = {}
for server in ["http://rest.ensemblgenomes.org", "http://rest.ensembl.org"]:
    for ext in ["/info/species?", "/info/species?division=EnsemblPlants"]:
        r = requests.get(server+ext, headers={"Content-Type": "application/json"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        else:
            decoded = r.json()
            for sp in decoded['species']:
                species = sp['name']
                ensembls[species] = server
                print(sp)

# Get the server for each of the species of interest
pre_list = ["homo_sapiens", "anolis_carolinensis", "ciona_intestinalis", "caenorhabditis_elegans",
            "gallus_gallus", "latimeria_chalumnae", "drosophila_melanogaster", "petromyzon_marinus",
            "monodelphis_domestica", "xenopus_tropicalis", "danio_rerio"]

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
            # if chrom['coord_system'] == 'chromosome' and not (chrom['name'].startswith('Un') or
            #                                                   chrom['name'].startswith('un') or
            #                                                   chrom['name'] in ['MT', 'cutchr', 'Mt', 'Pt'] or
            #                                                   'random' in chrom['name'] or 'hap' in chrom['name']):
            if chrom['coord_system'] == 'chromosome' and not (chrom['name'].startswith('Un') or
                                                              chrom['name'].startswith('un') or
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


def download_seq(ensembls):
    """
    Downloads transcript sequences given a transcript ID
    :param ensembls: pandas dataframe with genes annotation
    :return:
    """
    plant_fasta = "all_seqs.fa"
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
        symbol = seq_record['symbol']
        if symbol is np.nan:
            symbol = gene
        # Check if the gene has already been downloaded
        seq_name = "{}|{}|{}|{}|{}|{}|{}".format(sp, chrom, gene, symbol,
                                                 start, end, strand)
        seq_name = seq_name.replace(' ', '--')
        #timer_to_flush += 1
        server = sp_server[sp]
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
        #if timer_to_flush == 200:
        #    all_seqs.flush()
        #    timer_to_flush = 0
        print(sp, symbol)
    all_seqs.close()


print("Getting assembly information")
if not os.path.exists('chromosomes.csv'):
    assemblies = pd.concat(map(get_assembly, sp_server))
    # Save to file
    assemblies.to_csv('chromosomes.csv')
else:
    assemblies = pd.read_csv('chromosomes.csv')

print("Downloading annotations")
if not os.path.exists('annotation.csv'):
    annotation_table = get_annotation(assemblies)
    annotation_table.to_csv('annotation.csv')
else:
    annotation_table = pd.read_csv('annotation.csv')

protein_coding = annotation_table['biotype'] == 'protein_coding'
coding_annotation = annotation_table.loc[protein_coding]

# Remove overlapping genes
non_overlap_list = []
reindex = coding_annotation.set_index(['species', 'chromosome'])
if not os.path.exists("coding_no_overlap.csv"):
    for sp, chromosome in set(reindex.index):
        chrom_tab = coding_annotation.loc[(coding_annotation['species'] == sp) &
                                          (coding_annotation['chromosome'] == chromosome)]
        non_overlap_list.append(remove_overlaps(chrom_tab))
    coding_no_overlap = pd.concat(non_overlap_list)
    coding_no_overlap.to_csv('coding_no_overlap.csv')
else:
    coding_no_overlap = pd.read_csv("coding_no_overlap.csv", index_col=0)
coding_no_overlap.loc[:, 'length'] = coding_no_overlap['end'] - coding_no_overlap['start']
coding_no_overlap.loc[:, 'length'] = coding_no_overlap['length'].astype(int)
coding_no_overlap = coding_no_overlap.loc[:, ['start', 'end', 'length', 'strand',
                                              'chromosome', 'species', 'external_name']]
coding_no_overlap.loc[:, 'chromosome'] = coding_no_overlap['chromosome'].astype(str)

coding_no_overlap.reset_index(inplace=True)
coding_no_overlap.rename(columns={'gene_id': 'acc', 'external_name': 'symbol'}, inplace=True)
null_symbols = coding_no_overlap['symbol'].isnull()
coding_no_overlap.loc[null_symbols, 'symbol'] = coding_no_overlap.loc[null_symbols, 'acc']
coding_no_overlap.set_index('acc', inplace=True)
coding_no_overlap.to_csv("genes_parsed.csv")

if os.path.exists("all_seqs.fa"):
    downloaded = [gene.id.split("|")[2] for gene in SeqIO.parse("all_seqs.fa", "fasta")]
else:
    downloaded = []
to_download = coding_no_overlap.loc[~(coding_no_overlap.index.isin(downloaded))]
download_seq(to_download)
