from ftplib import FTP
import os

species = [x.strip('\n') for x in open('species.txt').readlines()]

for sp in species:
    print(sp)
    # Connect to the NCBI
    ncbi = FTP("ftp.ncbi.nih.gov")
    ncbi.login()
    # Get into each species' directory
    ncbi.cwd('genomes/{}/Assembled_chromosomes/gbs/'.format(sp))
    # Go through the chromosomes
    for chrom in ncbi.nlst():
        # Only download annotations from reference assembly
        if 'ref' in chrom:
            chrom_name = chrom.split('_')[-1].split('.')[0]
            # And only for assembled chromosomes
            if chrom_name not in ['unplaced', 'unlocalized', ]:
                # Download chromosome annotation
                if not os.path.exists(chrom.replace('.gz', '')):
                    ret_cmd = "RETR {}".format(chrom)
                    ncbi.retrbinary(ret_cmd, open(chrom, 'wb').write)
                    # Extract annotation
                    extract = 'gunzip -d {}'.format(chrom)
                    os.system(extract)
    print("Downloading proteome")
    ncbi = FTP("ftp.ncbi.nih.gov")
    ncbi.login()
    # Go to protein directory for each species
    ncbi.cwd("genomes/{}/protein/".format(sp))
    if not os.path.exists("{}.fa".format(sp)):
        # Download proteome
        ncbi.retrbinary("RETR protein.fa.gz", open('protein.fa.gz', 'wb').write)
        # Extract proteome
        os.system("gunzip protein.fa.gz")
        os.rename("protein.fa", "{}.fa".format(sp))
