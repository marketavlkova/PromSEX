#!/bin/bash
### Script incorporating all steps of promoter search and extraction
### created by Marketa Vlkova in 21-December-2018
### example to run: ./PromSEX.sh

### create input directory
if [ ! -d "./input" ]
then
  mkdir input
fi
### create output directory
if [ ! -d "./output" ]
then
  mkdir output
fi

### extract information from desired columns in RegulonDB file
printf "Extracting information from Regulon DB file\n"
./RegDBextract.awk -v COLS=2,3,4,6,8 input/MG1655_s70_promoters_RegDB.txt > output/RDBex.csv

### search and extract intergenic regions from K12 genome stated in your local Regulon DataBase
### here using only promoters with Strong confidence level
python3 SEXreference.py output/RDBex.csv input/MG1655_genome.gb Strong

### create input/Database directory if it doesn't exist
if [ ! -d "./input/Database" ]
then
  mkdir input/Database
  MAKEDB="Yes"
else
  MAKEDB="No"
fi
### create output directory for blast results
if [ ! -d "./output/blast_results" ]
then
  mkdir output/blast_results
fi

### blast for homologous intergenic regions among environmental isolates
./Blast.sh $MAKEDB ./output/*.fasta

### create output directories for intergenic regions of environmental isolates
if [ ! -d "./output/blasted_promoters" ]
then
  mkdir output/blasted_promoters
  mkdir output/blasted_promoters/by_promoter
fi

### extract blasted intergenic regions from environmental strains
### if the hits are at least 100 and 200 bp shorter then apropriate K12 intergenic region
for DIF in {1..2}
do
  printf "\nAnalysing blast hits that are at most ${DIF}00nt shorter than K12 reference"
  python3 SEXblasted.py output/*.fasta output/blast_results/ input/marketa_genomes/ "${DIF}00"
done

### create output directory for alignments
if [ ! -d "./output/alignments" ]
then
  mkdir output/alignments
fi

### align homologous intergenic regions and calculate sequence identity values
for DIF in {1..2}
do
  printf "\nAligning hits that are at most ${DIF}00nt shorter than K12 reference\n"
  ./Align.sh output/blasted_promoters/by_promoter/ "${DIF}00"
done

### extract total identity values in a separated file
for DIF in {1..2}
do
  printf "Extracting total identity values from Identities-${DIF}00.txt\n"
  ./IDextract.awk "output/Identities-${DIF}00.txt" > "output/TotalIdent-${DIF}00.csv"
done

### get files of aligned promoters having each sequence on one line only
for DIF in {1..2}
do
  for file in output/alignments/*"${DIF}"00*.aln
  do
    OUT="${file%.aln}.txt"
    ./JoinLines.awk $file > $OUT
  done
done

### get number of sequence variants for each promoter among the strains
printf "Pulling out numbers of existing promoter versions\n"
for DIF in {1..2}
do
  cat output/alignments/*"${DIF}"00*.txt | ./PromVerCount.awk > "output/VariantCounts-${DIF}00.csv"
done

### get number of segregating sites for each promoter alignment
printf "Pulling out information about segregating sites\n"
for DIF in {1..2}
do
  python3 SeqSites.py output/alignments/ "${DIF}00"
done
