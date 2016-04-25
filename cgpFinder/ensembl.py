import requests, sys
import biomart as bm

server = "http://rest.ensembl.org"
chromF = open("ensembl/chromosomes.csv", 'w')
chromF.write("species,chromosome,length\n")

assemblyF = open("ensembl/assemblies.csv", "w")
assemblyF.write("species,accession,date,name\n")

for sp in open("extras/species.txt"):
    spp = sp.lower().strip("\n")
    ext = "/info/assembly/{}?".format(spp)
    r = requests.get(server+ext, headers={"Content-Type": "application/json"})
    try:
        if not r.ok:
            r.raise_for_status()
            sys.exit()

        decoded = r.json()
        sp_dict = {}
        assembly = {}
        assembly_str = "{},{},{},{}\n".format(spp,
                                              decoded['assembly_accession'],
                                              decoded['assembly_date'],
                                              decoded['assembly_name'])
        assemblyF.write(assembly_str)

        for chromosome in decoded['top_level_region']:
            if chromosome['coord_system'] == 'chromosome' and not \
                    chromosome['name'] == 'MT':
                chrm_line = "{},{},{}\n".format(spp, chromosome['name'],
                                                chromosome['length'])
                chromF.write(chrm_line)
    except requests.HTTPError:
        print("{} is not available".format(spp))
assemblyF.close()
chromF.close()

