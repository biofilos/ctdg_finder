## Preprocessing
### From NCBI genomes
This section considers the case where the user wants to select a set of species from the [NCBI genomes](http://ftp.ncbi.nih.gov/genomes) database.  
In order to use CTDGFinder with custom databases, a set of files must me generated. The input files are:
* `proteomes.fa`: Fasta file containing most of the protein sequences of the species to be considered (the preprocessing script will attempt to download sequences not found in `proteomes.fa`, but it will severely affect performance, and generate unnecessary traffic on the NCBI server). Importantly, the accession number must be in the field 4 of the sequence name (fields being separated by "|"; this is the default format of protein sequences downloaded from the NCBI). 
* `genomes.gb`: GenBank file with genome annotations. Must contain CDS features, which will be considered for genes extraction. 
With both files in the same directory, run `preprocessing.py` (in the extras directory). It receives one argument, which is the number of cores required for the parsing process.