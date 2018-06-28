# CTDGFinder. A comparative genomics approach to gene clusters

 
The software package is a python3 script (`ctdg_finder.py`), which runs the homology search (via HMMeR), process coordinates information, and finds clusters of gene duplicates (ctdg).  
CTDGFinder expects a directory structure and files. A description of how to build them can be found in the [preprocessing manual](manual/preprocessing.mg)

### Installation
Before installing ctdgFinder, the following python modules should be installed:  

* numpy 
* pandas 
* numexpr
* scipy
* scikit-learn
* biopython
* Aditionally, HMMeR must be installed.  

### Running instructions
ctdgFinder is run as a python script (`python ctdg_finder.py`)  
The following are the minimal options that should be specified  

* \-(-n)ame_family: name of the gene family  
* \-(-r)ef_seq: path to the fasta file used as query  
* \-(-b)last_samples: number of samples to build the empirical distribution of paralogs.  
* \-(-d)b: database: directory whith the required blast database and annotation tables (all_seqs.fa, genes_parsed.csv, etc)  

A minimal run would look like this: `$ python ctdg_finder.py --ref-seq path/to/reference.fa --name_family gene_family_name --blast_samples 1000 --db mini`  

Optional parameters  

* \-(-s)p: If specified, the analysis will be performed only with the specified species. Note: the species name should be typed in qotes ("Homo sapiens"). If more than one species is going to be specified, an independent --sp option should be used (e.g. \--sp "Homo sapiens" \--sp "Mus musculus"). If not specified, the analysis wll be run with all the available 23 species.  
* \-(-o)ut_dir: directory where the results of the run are going to be stored (default: ctdg_out/).  
* \-(-c)pu: Number of cores used (default: 1 core).  
* \-(-i)terative: Run in iterative mode. In iterative mode, an extra BLASTP search will be performed at the beginning using the user-provided query sequence(s), and the complete sequences of all the BLASTP hits will be used as a compound query for CTDGFinder.
* \-(-e)value: Evalue threshold for BLASTP. If none is specified, 0.001 is used by default.
* \-(-p): Accession number to use as query, instead of fasta query. Useful for large scale analyses 

### Output

Under the ctdg_out directory (or the directory specified by \--out_dir), a directory for each gene family (specified by \--name_family) (name_family) will be created. This is useful in case several analyses are going to be performed. In case an analysis is going to be run for a gene family that already exists in the output direcrory, ctdgFinder will terminate.  
Under the gene family directory, the results from the run are found in the "report" directory. The "report" directory contain the following files:  

* `name_family_genes.csv`: genes found in the analysis. Most of the columns are self-explanatory. The "cluster" column contains the name of the cluster where the gene was found, and the "order" column describes the order (according to the genomic coordinates) of the gene in the ctdg. If the "cluster" and "order" columns have 0 (zero), the gene was found as single copy in that chromosome. If the value of "cluster" and "order" is "na", the gene was found in a chromosome with cluster candidates, but were discarded by either the meanshift or the statistical sampling steps.
* `name_family_numbers.csv`: Coordinates of the cluster candidates. The values of the "cluster" column follow the convention for the name_family_genes.csv file. The column "perc95_chrom" and "perc95_gw" columns describe the 95th percentile of the paralogs empirical distribution that was sampled in the analysis.
* `name_family_numbers_clean.csv`: A version of the `name_family_numbers.csv` file witout rejected clusters ("na" value in the "cluster" column).
* `name_family.samples`: Table that describes the maximum number of paralogs in each cluster candidate. These are the distributions used to calculate the above mentioned percentiles.
