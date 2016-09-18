## Preprocessing
### From NCBI genomes
This section considers the case where the user wants to select a set of species from the [NCBI genomes](http://ftp.ncbi.nih.gov/genomes) database.  
In order to use CTDGFinder with custom databases, a set of files must be generated. The following two files are required:  

* `proteomes.fa`: FASTA file containing most of the protein sequences of the species to be considered (the preprocessing script will attempt to download sequences not found in `proteomes.fa`, but it will severely affect performance, and generate unnecessary traffic on the NCBI server). Importantly, the accession number must be in the field 4 of the sequence name (fields being separated by "|"; this is the default format of protein sequences downloaded from the NCBI). 
* `genomes.gb`: GenBank file with genome annotations. Must contain CDS features, which will be considered for genes extraction.  

With both files in the current directory, run `preprocessing.py` (from the extras directory). It receives one argument, which is the number of cores required for the parsing process.  
The script will generate three files that will be used in CTDGFinder (genes_cds.csv and mini_table.csv are intermediate files, useful only for debugging purposes):  

* `chromosomes.csv`: Contains chromosome information (species, chromosome name, taxonomy ID, accession number, assembly code, and length).
* `genes_parsed.csv`: Gene annotation. Contains accession number, genomic coordinates, species name, and symbol and gene names for non-overlapping genes.
* `all_seqs.fa`: FASTA file containing protein sequences, sequence names, and associated information for all the genes in the `genes_parsed.csv` file.

Now, a BLAST database has to be generated based on the `all_seqs.fa` file. In order to do that, run `makeblastdb -dbtype prot -in all_seqs.fa`.  

The next step consists on running an all vs. all BLASTP using all the sequences in `all_seqs.fa`. In order to automate the process, run `all_blasts.py` (from the extras directory). A sample run would look like: `python extras/all_blasts.py all_seqs.fa 8` where the first argument points to `all_seqs.fa` (and its associated BLAST database files), and the second argument is the number of cores used in the computation. It will generate a `blasts` directory. The files with a `.blast_out` extension will be used for the dictionary generation (the `done` directory contains raw BLAST output files and intermediate parsing files)

The last step consists on generating JSON files containing the BLAST hits for each gene in the same chromosome and their associated e-value. This task is accomplished by the `hits_dictionary.py` script (in the `extras` directory). A sample run will look like `python extras/hits_dictionary.py blasts genes_parsed.csv`, where the firs argument is the directory containing the `.blast_out` files, and the second argument should point to the `genes_parsed.csv` file.

### From Scratch
CTDGFinder can be run using custom databases. The format and content of a minimal example of the necessary files follows (Note: space characters should be avoided in all files):  

* `chromosomes.csv`: comma-delimited file. The index (first column) can be an arbitrary unique string, which will not be used for the analysis. The columns *sp* (species name), *chromosome* (chromosome name), and *length* (chromosome length) must be included.
* `genes_parsed.csv`: comma-delimited file. The index should have the name *acc*. It will be the accession number (or other unique identifier, which has to be consistent in this file and the sequence file). The rest of the columns are: *start*, *end* (start and end coordinates), *strand* (strand using the notation 1, -1), *chromosome* (chromosome name), *species* (species name), *symbol* (short gene name), *length* (length of the gene)
* `all_seqs.fa`: FASTA file. The name of each sequence must comply with the following format: \>species_name|chromosome|accession|symbol|start|end|strand  

Importantly, the information in each sequence name should match with the fields in `genes_parsed.csv` and `chromosomes.csv`. With these three files, you can continue using the scripts `all_blasts.py` and `hits_dictionary.py` as described above. 