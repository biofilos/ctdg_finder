# CTDGFinder. A comparative genomics approach to gene clusters

 
The software package is a python3 script (`ctdg_finder.py`), which runs the homology search (via HMMeR), process coordinates information, and finds clusters of gene duplicates (ctdg).  
CTDGFinder expects a directory structure and files. A description of how to build them can be found in the [preprocessing manual](manual/preprocessing.mg)

### Installation
Before installing CTDGFinder, the following python modules should be installed:  

* numpy 
* pandas 
* HTSeq
* pybedtools (and install BedTools)
* scipy
* scikit-learn
### Running instructions
ctdgFinder is run as a python script (`python ctdg_finder.py config.json`) where `config.json` is a configuration file with the options relevant for an analysis. All the options in the configuration file have to be specified. An example configuration file is described below:  
  
```json
{
     "gff": "path to gff annotation",
  "chroms": "path to genome file (see below for details)",
  "paths": {
              "out_dir": "directory for output (will be created if necessary)"
            },
  "feat_type": "feature type to be clustered (e.g. 'gene')",
  "top_level_feat": "feature type to be considered top-level (e.g. 'chromosome')",
  "samples": "Number of regions to be sampled in the statistical assessment (WITHOUT QUOTES)",
  "percentile_threshold": "percentile cutoff: (more details in the publication)",
  "overwrite": "true or false (without quotes): if true, existing results will be overwritten"
}
```
#### Input files
##### GFF annotation
The Input GFF file should follow GFF basic specifications with some additions to the attributes column, as follows:
* The top-level sequences (e.g. chromosome) should have a *bandwidth* attribute. This should correspond to the average intergenic distance plus its standard deviation. A helper script to get this information can be found at `extras/preprocessing.py`
* The features to be clustered (e.g. gene) should have a *families* attribute. CTDGFinder can handle either features belonging to one family (e.g. orthologous groups), or many (domains, gene families, etc). If many gene families are to be used, they have to be separated with commas, with NO spaces.  

Example of valid GFf lines:  
```
chr1   Annotation         chromosome       1  1000000  .       +       .       ID=chr19;bandwidth=155348.4809  
chr1   Annotation         gene        107103   117102  .       +       .       ID=gene1;families=family1,family2  
chr1   Annotation         gene        202103   203012  .       +       .       ID=gene1;families=family3
```

#### Processing
CTDGFinder will output the progess of its run on the screen
```
Canis_familiaris   ...chr26_9: 23542658 - 23612283      D 3   P95: 0.0*
Canis_familiaris   ...chr26_10: 14264936 - 14682839     D 2   P95: 3
```
The first column is the base name of the GFF file (usually containing the species name). The second section corresponds to the coordinates of the cluster. Note that the chromosome name is followed by *\_#*, which is how each cluster is identified. The third column shows the number of duplicates, and the last column shows the percentile (in this case, the 95<sup>th</sup> percentile).  The asterisks at the left end of each line show clusters that. As shown in the example, the second cluster "chr26_10" is rejected because it has less duplicates than the percentile 95<sup>th</sup> for its samples.
### Output
CTDGFinder produces two GFF files, with the raw cluster information, and with merged clusters.  
* `clusters.gff`: This GFF file contains the features that were clustered, as well as *cluster* features, whcih are parents to the genes in that cluster.  
```
chr19   CTDGFinder      cluster 9087059 9252625 .       .       .       ID=chr19_1;length=165566;duplicates=7;percentile=0.0;families=PTHR26451
chr19   Annotation      gene    9087060 9095669 .       +       .       ID=ENSG00000170929;families=PTHR26451;biotype=protein_coding;Parent=chr19_1;order=1
chr19   Annotation      gene    9100406 9107475 .       -       .       ID=ENSG00000170923;families=PTHR26451;biotype=protein_coding;Parent=chr19_1;order=2

```
As it can be seen, the *cluster* feature has some attributes that describe the cluster. They are:  
* **Length**: Distance from the first feature in the cluster to the last one.
* **duplicates**: Number of clustered genes in the cluster.
* **percentile**: Value of the percentile taken from the random distribution for that cluster.
* **families**: Name of the families included in the clustered genes of the cluster.  

The feature GFF lines have a *Parent* attribute that connects them to the cluster they belong to, as well as an *order* attribute that describes the position of the gene in the cluster

* `merged_clusters.gff`: This file is very similar to `clusters.gff`, except that the *cluster* features do not have a *percentile* attribute. The reason for this is that a *cluster* feature in this file can be the product of merging various clusters together. Clusters of this nature exist. For example, the globin cluster on the human chromosome 11 is embeded in a larger olfaction receptor cluster ([Bulger et al. 1999](http://www.pnas.org/cgi/doi/10.1073/pnas.96.9.5129)). 