from bioservices import ensembl, UniProt, biomart
import pandas as pd
import io
import requests
import json
import numpy as np
import time
import sys

#Load gene table
gene_table = sys.argv[1]
b = biomart.BioMart()
b.host = 'www.ensembl.org'

datasets = b.datasets("ENSEMBL_MART_ENSEMBL")
# Load annotation
big_annotation = pd.read_csv(gene_table)
big_annotation.set_index(['species','chromosome'], inplace=True)
table_list = []
chunk_size = 100000
for sp, chrom in set(big_annotation.index):
    print("{}: {}".format(sp, chrom))
    sp_code = sp[0]+sp.split('_')[-1]
    for sp_server in datasets:
        if sp_code in sp_server.lower():
            break
    ids = big_annotation.loc[(sp,chrom), 'acc'].values
    id_chunks = [ids[i:i+chunk_size] for i in range(0,len(ids),chunk_size)]
    # Configure search
    #filters = b.filters(sp_server)
    #attributes = b.attributes(sp_server)
    for chunk in id_chunks:
        # set query
        b.new_query()
        b.add_dataset_to_xml(sp_server)
        # These are the columns I want
        attributes = ['ensembl_gene_id','uniprot_genename',
                      'interpro','go_id','pfam','family']
        for attr in attributes:
            b.add_attribute_to_xml(attr)
        id_str = ','.join(chunk)
        b.add_filter_to_xml('ensembl_gene_id',id_str)
        xml = b.get_xml()
        result = requests.get('http://useast.ensembl.org/biomart/martservice?query='+xml)
        tab_result = pd.read_csv(io.StringIO(result.text),sep='\t',header=None)
        table_list.append(tab_result)
    time.sleep(3)

tab_result = pd.concat(table_list)
tab_result.columns = ['ensembl','uniprot','interpro','go','pfam','panther']
tab_result.set_index('ensembl',inplace=True)

annotation_dict = {}
for gene in set(tab_result.index):
    annotation_dict[gene] = {}
    for column in tab_result:
        values = tab_result.loc[gene,column]
        if type(values) in [np.float64,float,str]:
            items = values
        else:
            items = values.values
            annotation_dict[gene][column] = items




