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
* The features to be clustered (e.g. gene) should have a *families* attribute. CTDGFinder can handle either features belonging to one family (e.g. orthologous groups), or many (domains, gene families,e tc). If many gene families are to be used, they have to be separated with commas, with NO spaces.  

Example of valid GFf lines:  
```
chr1   Annotation      chromosome      1       1000000        .       +       .       ID=chr19;bandwidth=155348.4809  
chr1   Annotation      gene    107103  117102  .       +       .       ID=gene1;families=family1,family2  
chr1   Annotation      gene    107103  117102  .       +       .       ID=gene1;families=family3
```

#### Processing
CTDGFinder will output the progess of its run on the screen
```
Canis_familiaris   ...chr26_9: 23542658 - 23612283      D 3   P95: 0.0*
Canis_familiaris   ...chr26_10: 14264936 - 14682839     D 2   P95: 3
```
The first column is the base name of the GFF file (usually containing the species name). The second section corresponds to the coordinates of the cluster. Note that the chromosome name is followed by *_#*, which is how each cluster is identified. The third column shows the number of duplicates, and the last column shows the percentile (in this case, the 95<sup>th</sup> percentile).  The asterisks at the left end of each line show clusters that. As shown in the example, the second cluster "chr26_10" is rejected because it has less duplicates than the percentile 95<sup>th</sup> for its samples.
### Output
CTDGFinder produces two 