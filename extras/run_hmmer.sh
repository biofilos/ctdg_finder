#!/usr/bin/env bash
set -e

script_dir=`dirname $0`
top_level=$1
proteome=$2
hmmer_db=$3
gff_in=$4
gff_out=$5
chroms_out=$6
cpus=$7



# Generate hmmer database
hmmpress "$hmmer_db"
# Run HMMsearch
echo Running HMMeR
mkdir -p hmmer
hmmsearch --cpu "$cpus" --noali -E 0.0001 -o hmmer/all.hmmer --domtblout hmmer/all.dom "$hmmer_db" "$proteome"
# Parse hmmer output
echo Parsing HMMeR output
python "$script_dir"/annotate_gff.py "$top_level" "$gff_in" hmmer/all.dom "$gff_out" "$chroms_out"