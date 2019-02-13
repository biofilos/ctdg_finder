#!/bin/sh

script_dir=`dirname $0`
db_dir=$1
cpus=$2

# Create directory structure
for i in hmmer hmmers graphs matrices sp_hmmers
do
	mkdir -p $db_dir/$i
done

# Calculate bandwidth parameter
echo Calculating bandwidth
python $script_dir/bandwidth.py $db_dir $cpus
# Run HMMsearch
echo Running HMMeR
hmmsearch --cpu $cpus --noali -E 0.0001 -o $db_dir/hmmer/all.hmmer --tblout $db_dir/hmmer/all.tbl --domtblout $db_dir/hmmer/all.dom panther/panther.hmm $db_dir/seqs_hmmer.fa
# Parse hmmer output
echo Parsing HMMeR output
python $script_dir/parse_big_hmmer_strict.py $db_dir/hmmer/all.dom $db_dir/hmmer/pfam.json $db_dir/hmmer/hmmer_filtered.csv
# Generate chromosome-specific JSON files
echo Generating dictionaries
python $script_dir/hits_dictionary.py $db_dir $cpus
# Generate homology graphs
echo Generating homology graphs
python $script_dir/domains_net.py $db_dir/hmmer/pfam.json $db_dir/genes_parsed.csv $db_dir/graphs $cpus
# Generate homology matrices
echo Generating homology matrices
python $script_dir/parse_matrices.py $db_dir/graphs matrices

