import json
from glob import glob
import sys
import pandas as pd

# Change directory to location of blast files
# os.chdir(sys.argv[1])
blast_folder = sys.argv[1]
# Load proteome annotation
all_genes = pd.read_csv(sys.argv[2])
# Make sure the chromosome is a string type
all_genes['chromosome'] = all_genes['chromosome'].astype(str)

# Parse parsed blasts
for table in glob('{}/*.blast_out'.format(blast_folder)):
    # Get species from file name
    sp = table.split('.')[0].split('/')[-1]
    print(sp)
    # Extract genes from that species
    sp_genes = all_genes.loc[all_genes['species'] == sp]
    # Load blast for that species
    sp_table = pd.read_csv(table)
    sp_table['chromosome'] = sp_table['chromosome'].astype(str)
    chromosomes = set(sp_genes['chromosome'].values)

    # Initialize dictionary with chromosome names
    table_dict = {chrom: {} for chrom in chromosomes}
    # Go through the chromosomes
    for chrom in chromosomes:
        print(chrom)
        # Get annotation from all the genes in the chromosome
        queries = sp_genes.loc[sp_genes['chromosome'] == chrom, 'acc'].values
        # Subset the blast with queries and subjects from that chromosome
        chrom_table = sp_table.loc[(sp_table['query'].isin(queries)) &
                                   (sp_table['chromosome'] == chrom)]

        for query in queries:
            value_array = list(chrom_table.loc[(chrom_table['query'] == query) &
                                               (chrom_table['species'] == sp), ['prot_acc', 'evalue']].values)
            value_list = [list(x) for x in value_array]
            table_dict[chrom][query] = value_list

#     # Save result to JSON file
    with open(table.replace('.blast_out', '.json'), 'w') as json_out:
        json_str = json.dumps(table_dict)
        json_out.write(json_str)
