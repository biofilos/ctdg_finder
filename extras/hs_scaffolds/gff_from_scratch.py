import glob
from concurrent import futures
import numpy as np
import os
import pandas as pd
import time
from Bio import SeqIO, Entrez

CPUS = 6
def parse_gff(in_file, sp):
    """
    Parses a GFF3 file, and returns a chromosome and CDS table
    :param in_file: GFF# file name
    :param sp: Species name
    :return: two pd.Dataframes
    """
    regions, cds = [], []
    # sp = 'hs_GRCh38.p7'
    with open(in_file) as fileIn:
        for line in fileIn:
            if not line.startswith("#"):
                record = line.rstrip().split('\t')
                chromosome = record[0]
                type_feat = record[2]
                attr_dict = {}
                for attr in record[8].split(';'):
                    fields = attr.split("=")
                    attr_dict[fields[0]] = fields[1]
                if type_feat == 'CDS' and 'protein_id' in attr_dict:
                    start, end = record[3], record[4]
                    strand = record[6]
                    accession = attr_dict['protein_id']
                    if 'gene' in attr_dict:
                        symbol = attr_dict['gene']
                    else:
                        symbol = accession

                    if 'product' in attr_dict:
                        product = attr_dict['product']
                    else:
                        product = attr_dict['gene']
                    data = [sp, chromosome, start, end, strand, symbol, accession, product]
                    cds.append(data)
                elif type_feat == 'region' and 'Name' in attr_dict:
                    region_chrom = attr_dict['Name']
                    region_id = record[0]
                    region_length = record[4]
                    regions.append([sp, region_id, region_chrom, region_length])
    # Genes table
    feat_table = pd.DataFrame(cds, columns=['species', 'chromosome', 'start', 'end', 'strand',
                                            'symbol', 'acc', 'name'])
    feat_table.loc[:,'strand'] = feat_table['strand'].map(lambda x: 1 if x == '+' else -1)
    feat_table.loc[:,'start'] = feat_table['start'].astype(int)
    feat_table.loc[:,'end'] = feat_table['end'].astype(int)
    feat_table.loc[:,'length'] = abs(feat_table['end'] - feat_table['start'])

    # Chromosomes table
    chrom_table = pd.DataFrame(regions, columns=['sp', 'acc', 'chromosome', 'length'])
    chrom_table = chrom_table.loc[chrom_table['acc'].isin(set(feat_table['chromosome']))]

    return feat_table, chrom_table

# feat, chrom = parse_gff("hs-sc_GRCh38.p2.gff3", "hs_38-sc")

feats, chroms = [], []

for genome in glob.glob("*.gff3"):
    sp = genome.split("_")[0]
    feat, chrom = parse_gff(genome, sp)
    feats.append(feat)
    chroms.append(chrom)

cds_table = pd.concat(feats).set_index('acc')
chromosomes = pd.concat(chroms)

####


def sp_loop(in_table, columns):
    """
    Filters a table by a number of columns
    :param in_table: Input table
    :param columns: Column names to perform the filter
    :return: list of tables for each combination of column names
    """
    # Initialize the list of tables
    list_results = []
    for col_set in np.unique(in_table.set_index(columns).index):
        filtered_table = in_table
        cols_to_filter = {col_name: col_value for col_name, col_value in zip(columns, col_set)}
        for var in cols_to_filter.keys():
            filtered_table = filtered_table.loc[filtered_table[var] == cols_to_filter[var]]
        list_results.append(filtered_table)
    return list_results


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
    chrom_only.drop_duplicates('acc', inplace=True)
    chrom_only.set_index('acc', inplace=True)

    first_to_one_to_last = chrom_only.index.values[:-1]
    second_to_last = chrom_only.index.values[1:]

    overlapping = chrom_only.loc[second_to_last,
                                 'start'] < chrom_only.loc[first_to_one_to_last, 'end']

    same_strand = chrom_only.loc[second_to_last, 'strand'] == chrom_only.loc[first_to_one_to_last, 'strand']

    for_removal = same_strand & overlapping
    for_removal_indices = for_removal.loc[for_removal].index
    non_overlapping_table = chrom_only.drop(for_removal_indices)
    return non_overlapping_table


def remove_overlaps_sp(sp):
    """
    Loads a csv named sp.csv in the current directory.
    Adds a length column
    Removes genes with duplicated start coordinates.
    :param sp:
    :return:
    """
    sp_table_list = []
    # Load species table
    sp_table = out_data.loc[out_data.species == sp]

    # Set indices so that a data frame per chromosome with ordered genes by position and length can be easily selected
    sp_table.sort_values(['start', 'length'], inplace=True)

    # Remove entries that have the same start position and leave the longest one
    sp_table.drop_duplicates('start', inplace=True)

    # Remove overlapping genes
    # Get list of tables per chromosome
    chrom_tables = sp_loop(sp_table, ['species', 'chromosome'])
    # Get list of tables (one per chromosome) without overlapping genes
    non_overlapping = []
    for chrom_table in chrom_tables:
        before = chrom_table.shape[0]
        no_overlaps = remove_overlaps(chrom_table)
        after = no_overlaps.shape[0]
        while before != after:
            before = no_overlaps.shape[0]
            no_overlaps = remove_overlaps(no_overlaps)
            after = no_overlaps.shape[0]
        non_overlapping.append(no_overlaps)
    parsed_sp_table = pd.concat(non_overlapping)

    return parsed_sp_table


def dict_seq(fasta):
    seq_dict = {}
    for line in open(fasta):
        line = line.strip('\n')
        if line.startswith('>'):
            name = line.strip('>').split('|')[3]
            seq_dict[name] = ''
        else:
            seq_dict[name] += line
    return seq_dict


def seq_column(accession):
    spp = all_sp_table.loc[accession, 'species']
    chromosome = all_sp_table.loc[accession, 'chromosome']
    try:
        seq = seq_dic[accession]
    except KeyError:
        try:
            print("Downloading {0} from the NCBI ({1} ch: {2})".format(accession, spp, chromosome))
            seq = str(SeqIO.read(Entrez.efetch(db='protein', id=accession, rettype='fasta'), 'fasta').seq)
        except:
            time.sleep(1)
            print("Downloading {} from the NCBI again".format(accession))
            seq = str(SeqIO.read(Entrez.efetch(db='protein', id=accession, rettype='fasta'), 'fasta').seq)
        record = pd.Series({accession: seq})
    return record


# Parse Genbank
out_data = cds_table
out_data =  out_data.loc[~(out_data['name'].map(lambda x: 'isoform X' in x))]
# Save full annotation
out_data.loc[:, 'species'] = out_data['species'].map(lambda x: x.replace(' ','_'))
out_data.to_csv('genes_cds.csv')
# Remove some columns
mini = out_data.loc[:, ['species', 'chromosome', 'symbol', 'start', 'end', 'length', 'strand']]
mini.to_csv('mini_table.csv')
del mini



if os.path.exists('genes_parsed.csv'):
    all_sp_table = pd.read_csv('genes_parsed.csv')
else:
    with futures.ProcessPoolExecutor(CPUS) as pool:
        sp_tables = pool.map(remove_overlaps_sp, list(set(out_data.species.values)))
    all_sp_table = pd.concat(sp_tables)
    all_sp_table.to_csv('genes_parsed.csv')
    del out_data

print("GenBank annotation parsed")

# Generate fasta file
genes_selected = all_sp_table
genes_selected.reset_index(inplace=True)


seq_dic = dict_seq('proteomes.fa')
print("Sequence dictionary, generated")

#gene_ids = list(seq_dic.keys())
#for gene in gene_ids:
#    if gene not in all_sp_table['acc'].values:
#        del seq_dic[gene]
#    else:
#        print(gene)
print("Sequence file, parsed")
if not os.path.exists("all_seqs.fa"):
    with open('all_seqs.fa', 'w') as seq_db:
        for index in genes_selected.index:
            seq_table = genes_selected.loc[index]
            spp = seq_table['species'].replace(' ', '_')
            chromosome = seq_table['chromosome']
            accession = seq_table['acc']
            symbol = seq_table['symbol']
            start = seq_table['start']
            end = seq_table['end']
            strand = seq_table['strand']
            try:
                seq = seq_dic[accession]
                del seq_dic[accession]
            except:
                try:
                    print("Downloading {0} from the NCBI ({1} ch: {2})".format(accession, spp, chromosome))
                    seq = str(SeqIO.read(Entrez.efetch(db='protein', id=accession, rettype='fasta'), 'fasta').seq)
                except:
                    time.sleep(1)
                    print("Downloading {} from the NCBI again".format(accession))
                    seq = str(SeqIO.read(Entrez.efetch(db='protein', id=accession, rettype='fasta'), 'fasta').seq)
            structure = '>{0}|{1}|{2}|{3}|{4}|{5}|{6}\n{7}\n'.format(spp, chromosome, accession,
                                                                     symbol, start, end, strand, seq)
            seq_db.write(structure)
            seq_db.flush()
    seq_db.close()

print("Sequence file, generated")
del seq_dic
del all_sp_table
del genes_selected

assembly_table = chromosomes
assembly_table['sp'] = assembly_table['sp'].map(lambda x: x.replace(' ','_'))

assembly_table.to_csv('chromosomes.csv')
print("Assembly file, generated")


print("DONE")
