#!/bin/bash
### Script for alignment and sequence identity calculation of each intergenic region among strains

### example to run: bash Align.sh <path to fasta files of extracted intergenic regions> <number of which alignment length can be smaller than k12 promoter to be included>
### e.g.: bash Align.sh output/blasted_promoters/by_promoter/ 100

### multiple sequence alignment of all found homologous intergenic regions
for file in $1*$2*.fasta
do
  NO_SUFIX="${file%.*}"
  NAME="${NO_SUFIX#*/*/*/}"
  echo "Aligning $NAME"
  t_coffee $file -mode procoffee -run_name=output/alignments/$NAME &> /dev/null
done

### delete if the file already exists
if [ -f "output/Identities-$2.txt" ]
then
  rm output/Identities-$2.txt
fi

### calculate identities for aligned sequences
echo "Calculating identities for aligned sequences"
for align in output/alignments/*$2*.aln
do
  t_coffee -other_pg seq_reformat -in $align -output sim 2> /dev/null
done >> output/Identities-$2.txt
