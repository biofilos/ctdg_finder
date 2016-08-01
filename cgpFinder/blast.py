import numpy as np
import os
import pandas as pd
from pandas.computation.eval import eval


def blast_filter(series):
    """
    This function calculates the length ratio of subject and query in a blast hits table output
    :param series: row of blast output (pd.Series)
    :return: length ratio
    """
    cols = ['s_len', 'q_len']
    ratio = min(series[cols]) / max(series[cols])

    return ratio


def exe(cpus, ref, subject, out_file, evalue=0.0001):
    """
    Blast a reference sequence againsts all the sequences in a blast database
    :param evalue: E-Value threshold in Blast
    :param subject: Blast database
    :param cpus: Number of CPUs to be used by Blast
    :param ref: Query
    :param out_file: Output for blast
    :return: pd.DataFrame with filtered sequences
    """
    # Set blast command
    cmd = "blastp -query {0} " \
          "-db {1} " \
          "-out {2} " \
          "-evalue {3} " \
          "-outfmt '6 qseqid sseqid qlen slen qstart qend sstart send length gaps gapopen evalue bitscore' " \
          "-num_threads {4}".format(ref, subject, out_file, evalue, cpus)
    # Run blast
    os.system(cmd)


def parse(out_file, acc_col, tab=True, sp_list=[], for_dict=False):
    """
    Parse Blast output, and get annotation from the sequence names
    :param sp_list: optional species list to restrict the analysis
    :param out_file: output
    :param acc_col: sequence name field (separated by |) that contains the accession number
    :param tab: parse annotation from sequence name
    :param for_dict: Analysis for building the json files?
    :return: data frame
    """
    blast_size = os.stat(out_file).st_size / 1000
    print("Blast results size: {} kb".format(blast_size))
    # Read output
    if os.stat(out_file).st_size == 0:
        print("Blast output is empty")
        no_records = pd.DataFrame()
        return no_records
    else:
        result_tab = pd.read_table(out_file, header=None)
        result_tab.columns = ['query', 'subject', 'q_len', 's_len',
                              'q_start', 'q_end', 's_start', 's_end',
                              'aln_len', 'gaps', 'gap_open', 'eval', 'bitscore']
        # Sort the results by query, subject and alignment length, so that they can be filtered
        if tab and not for_dict:
            if len(np.unique(result_tab['query'].values)) > 1:
                new_query_name = result_tab['query'].values[0]
                result_tab['query'] = new_query_name

        # Filter out hits if they are less than one third (or greater
        #  than three times) the query length
        result_tab['len_ratio'] = result_tab.apply(blast_filter, axis=1)
        result_tab.sort_values(by=['query', 'subject', 'len_ratio'], inplace=True, ascending=True)
        result_tab.drop_duplicates(subset=['query', 'subject'], keep='last', inplace=True)

        filtered_blast = result_tab.loc[result_tab.len_ratio >= 0.3]
        # Save filtered output if run is for family definition
        filtered_blast.to_csv("{}_filtered".format(out_file))
        # If more than one query was used, rename all the queries with the same name
        if tab:
            try:
                filtered_blast['query_acc'] = filtered_blast['query'].map(lambda x: x.split("|")[acc_col])
            except IndexError:
                filtered_blast['query_acc'] = filtered_blast['query'].map(lambda x: x.split("|")[0])
            # Attach the query name to the subject, to include it as a column
            filtered_blast.loc[:, 'subject'] = filtered_blast.query_acc + '|' + filtered_blast.subject
            if for_dict:
                filtered_blast.loc[:, 'subject'] = filtered_blast.query_acc + '|' + filtered_blast['eval']
            # Extract table from subject names
            fields = map(lambda x: x.split('|'), filtered_blast.subject.values)

            # Set table
            sub_table = pd.DataFrame(list(fields), columns=['query', 'species', 'chromosome', 'prot_acc',
                                                            'symbol', 'start', 'end', 'strand', 'evlaue'])

            # sub_table['species'] = sub_table['species'].apply(lambda x: str(x).replace('_', ' '))
            # Only consider hits from selected species
            if len(sp_list) > 0:
                sub_table = sub_table.loc[sub_table['species'].isin(sp_list)]
            sub_table['start'] = sub_table['start'].astype(int)
            sub_table['end'] = sub_table['end'].astype(int)
            sub_table['strand'] = sub_table['strand'].astype(int)
            sub_table['Length'] = abs(sub_table['end'] - sub_table['start'])
            sub_table.sort_values(by=['species', 'chromosome', 'start'], inplace=True)
            sub_table.to_csv("{}_out".format(out_file))
            # Return filtered blast results
            return sub_table
        else:
            return filtered_blast

