#!/bin/bash
### Script for alignment and sequence identity calculation of each intergenic region among strains

### example to run: bash Align.sh <path to fasta files of extracted genes>
### e.g.: bash Align.sh output/blasted_genes/by_promoter/

### multiple sequence alignment of all found homologous intergenic regions
for file in $1*.fasta
do
  NO_SUFIX="${file%.*}"
  NAME="${NO_SUFIX#*/*/*/}"
  echo "Aligning $NAME"
  t_coffee $file -type=dna -run_name=output/alignments/genes/$NAME &> /dev/null
done

### delete if the file already exists
if [ -f "output/GeneIdentities.txt" ]
then
  rm output/GeneIdentities.txt
fi

### calculate identities for aligned sequences
echo "Calculating identities for aligned sequences"
for align in output/alignments/genes/*.aln
do
  t_coffee -other_pg seq_reformat -in $align -output sim 2> /dev/null
done >> output/GeneIdentities.txt
