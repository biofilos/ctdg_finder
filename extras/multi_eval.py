import os
import pandas as pd
import sys
from glob import glob
# This script generates multiple JSON dictionary files for paralogs databse for different E-value thresholds

# It will use a Blast run with an E-value = 1

# Path to the *.blast_filtered files

blast_dir = sys.argv[1]


for evalue in [1]:#, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5]:
    for blast_file in glob("{}/*.blast_filtered".format(blast_dir)):
        print("{} ({})".format(blast_file.split('/')[-1], evalue))
        filtered_blast = pd.read_csv(blast_file, index_col=0)
        filtered_blast = filtered_blast.loc[filtered_blast['eval'] <= evalue]
        filtered_blast['query_acc'] = filtered_blast['query'].map(lambda x: x.split("|")[2])

        # Attach the query name to the subject, to include it as a column
        filtered_blast.loc[:, 'subject'] = filtered_blast.query_acc + '|' + filtered_blast.subject + '|' +\
                                           filtered_blast['eval'].astype(str)

        # Extract table from subject names
        fields = map(lambda x: x.split('|'), filtered_blast.subject.values)


        # Set table
        sub_table = pd.DataFrame(list(fields), columns=['query', 'species', 'chromosome', 'prot_acc',
                                                        'symbol', 'start', 'end', 'strand', 'evalue'])

        # sub_table['species'] = sub_table['species'].apply(lambda x: str(x).replace('_', ' '))
        # Only consider hits from selected species
        sub_table['start'] = sub_table['start'].astype(int)
        sub_table['end'] = sub_table['end'].astype(int)
        sub_table['strand'] = sub_table['strand'].astype(int)
        sub_table['Length'] = abs(sub_table['end'] - sub_table['start'])

        sub_table.sort_values(by=['species', 'chromosome', 'start'], inplace=True)

        base_path = blast_file.split('/')[:-1]
        out_path = "/".join(base_path) + "/" + str(evalue).replace(',', '_')
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        out_filename = blast_file.split('/')[-1].replace('filtered', 'out')
        sub_table.to_csv(out_path + "/" + out_filename)

