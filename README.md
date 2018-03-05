# CTDGFinder
## Database generation using HMMer

This version of CTDGFinder was modified to use HMMer instead of Blast. More recent advances in orthology calling are using HMM profiles to define orthologous groups and gene families (see [Panther](http://pantherdb.org/) and the [Quest of orthologs DB repository](https://questfororthologs.org/orthology_databases), and of course, [Pfam](https://pfam.xfam.org/)). For that reason, CTDGFinder wants to make it easier to use those resources on its pipeline.  
A substantial HMMsearch run has to be performed only once to annotate the proteomes with their domains, gene families, orthology groups, etc present in the HMM database. However, this potential downside can be leveraged for large-scale analyses. If the user wants to identify all the CTDGs present in a set of proteomes, they won't have to run HMMscan for each protein, but used the precomputed HMMscan search by using the new parameter `-p`, which instead of using a sequence as a query (do not use the `-f` parameter in this case), it will use the HMM accession specified by `-p` as a query. Given that a large-scale analysis can include several thousand HMMscan searches, this significantly improves the runtime of CTDGFinder on those cases. A test case example will be provided.  

### Preprocessing steps
All this steps are included in the script `extras/run_hmmer.sh`, but will be explained here in case they have to be tweaked by the used.  

#### Note on intergenic distances and overlapping genes.  
Since no `preprocessing.py` is going to be run, it is assumed that the genes are properly annotated. In the Blast version of `preprocessing.py`, overlapping genes were removed, to make sure no negative intergenic distances were called, and to reduce the chances of calling isoforms as genes. Although it is not very common, overlapping genes do exist in many species, and in this version, *no check for overlapping genes will be performed*. In case of encountering overlapping genes, the script `bandwith.py` will set its intergenic distance to zero.  

#### Directory structure
assuming that `ctdg_db` is the base directory where the annotation files reside (genes_parsed.csv, pre_chromosomes.csv, etc), the following directories are going to be created
* ctdg_db/hmmer
* ctdg_db/hmmers
* ctdg_db/graphs
* ctdg_db/matrices  

NOTE: Observe that chromosomes.csv has been changed to pre_chromosomes.csv. The reason is that in the next step, the bandwidth parameter for MeanShift is going to be calculated, and that step will generate chromosomes.csv which will contain a `bandwidth` column.  

#### Calculate bandwidth parameter for MeanShift  
In order to avoid calculating it on each run, it can be calculated for all chromosomes in one go, improving the performance of CTDGFinder

#### Run HMMsearch  
The gene family/orthologs source selected by the user (and properly formatted by using `hmmpress`) will be used as database, and a sequence file named seqs_hmmer.fa will be used as query. By default, the hmm file will be assumed to be stored as `ctdg_db/panther/panther.hmm` and `ctdg_db/seqs_hmmer.fa` where `ctdg_db` is the directory where all the annotation files are stored (chromosomes.csv, genes_parsed.csv, etc). The E-value used for the search is 0.0001. If the user wants to use a different E-value, it will have to be changed both in this step, as well as in the source code of CTDGFinder.  
Notes on annotation files structure:
* seqs_hmmer.fa has to have a header with the following format: >accession|species|chromosome|start|end|strand
* genes_parsed.csv requires the following columns: acc, species, chromosome, start, end, strand. acc is the accession of each gene, and should match the one in seqs_hmmer.fa
#### Parse HMMsearch output  
In this step, the HMMscan output will be parsed to generate a JSON dictionary containing each gene family as keys, and a list of all genes part of that gene family. **NOTE**: For my own research, I am using PantherDB, which is a gene family database. For that reason, I put an extra filter to the HMMscan results, so that only sequence hits covering 30% of both query and subject will be kept. If you use a domain-based database (e.g. Pfam), you might want to remove this constraint. Contact me if you want to do it in that way, because so far, this feature is hard-coded in CTDGFinder.

#### Generate chromosome-specific JSON files
These files contain each gene as keys, and a list of all the genes in that chromosomes that are in the same gene family as the key. Testing different parsing approaches, I found that reading individual files per chromosome on each of the sampling steps is the fastest way. Mainly because a species-JSON file can get very big, and would have to be copied on all the threads, which has shown to have a negative impact in performance.  

#### Generate homology graphs  
In order to generate the homology tables used in the sampling process of CTDGFinder, graphs will be generated where two genes will have an edge if they are in the same gene family. Since a gene can be in different gene families (or have more than one domain), this procedure can significantly reduce the size of the files used in the analysis.  

#### Generate homology matrices  
Finally, homology matrices will be generated. The general format consists on csv files where genes on each chromosome (for each species) are the indices of rows and columns, and the values of each cell (0, 1) show if a pair of genes are in the same gene fmaily (or share at least one domain).  


