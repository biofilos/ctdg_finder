from bioservices import import biomart
import sys

# Set biomart
b = biomart.BioMart()
b.host = 'biomart.vectorbase.org'
in_file = sys.argv[1]
species = [x.strip('\n') for x in open(in_file).readlines()]

for sp in species:
    dataset = sp[0]+sp.split("_")[1]+'eg_gene'

