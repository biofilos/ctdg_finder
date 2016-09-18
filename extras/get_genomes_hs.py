from ftplib import FTP
import os

builds = [x.strip('\n') for x in open('builds.txt').readlines()]

for build in [x for x in builds if not x.startswith("#")]:
    print(build)
    # Connect to the NCBI

    # Get into each species' directory
    ncbi = FTP("ftp.ncbi.nih.gov")
    ncbi.login()
    ncbi.cwd('genomes/Homo_sapiens/ARCHIVE/{}'.format(build))

    chromosomes = [x for x in ncbi.nlst() if x.startswith("CHR")]
    for chrom in chromosomes:
        ncbi.cwd(chrom)
        # Only download annotations from reference assembly
        for file in [x for x in ncbi.nlst() if 'alt' not in x and 'gbs' in x]:
            chrom_name = chrom.split('_')[-1].split('.')[0]

            # Download chromosome annotation
            out_file = file.replace('.gz', '')
            if not os.path.exists(out_file):
                ret_cmd = "RETR {}".format(file)
                ncbi.retrbinary(ret_cmd, open(file, 'wb').write)
                # Extract annotation
                extract = 'gunzip -d {}'.format(file)
                os.system(extract)
                os.rename(out_file, "{}_{}.gbs".format(build.lower(), chrom_name))
                break
        ncbi.cwd("..")

    print("Downloading proteome")
    ncbi = FTP("ftp.ncbi.nih.gov")
    ncbi.login()
    # Go to protein directory for each species
    ncbi.cwd("genomes/Homo_sapiens/ARCHIVE/{}/protein/".format(build))
    if not os.path.exists("{}.fa".format(build)):
        # Download proteome
        ncbi.retrbinary("RETR protein.fa.gz", open('protein.fa.gz', 'wb').write)
        # Extract proteome
        os.system("gunzip protein.fa.gz")
        os.rename("protein.fa", "{}.fa".format(build.lower()))
