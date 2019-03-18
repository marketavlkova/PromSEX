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
./Blast.sh $MAKEDB ./output/*.fasta output/blast_results/

### create output directories for intergenic regions of environmental isolates
if [ ! -d "./output/blasted_promoters" ]
then
  mkdir output/blasted_promoters
  mkdir output/blasted_promoters/by_promoter
fi

### extract blasted intergenic regions from environmental strains
### if the hits are at least 100 and 200 bp shorter then apropriate K12 intergenic region
printf "\nAnalysing blast hits that are at most 100nt shorter than K12 reference"
python3 SEXblasted.py output/*.fasta output/blast_results/ input/marketa_genomes/ output/blasted_promoters/ 100

### create output directory for alignments
if [ ! -d "./output/alignments" ]
then
  mkdir output/alignments
fi

### align homologous intergenic regions and calculate sequence identity values
printf "\nAligning hits that are at most 100nt shorter than K12 reference\n"
./Align.sh output/blasted_promoters/by_promoter/ 100

### extract total identity values in a separated file
printf "Extracting total identity values from Identities-100.txt\n"
./IDextract.awk output/Identities-100.txt > output/TotalIdent-100.csv

### get files of aligned promoters having each sequence on one line only
for file in output/alignments/*.aln
do
  OUT="${file%.aln}.txt"
  ./JoinLines.awk $file > $OUT
done

### get number of sequence variants for each promoter among the strains
printf "Pulling out numbers of existing promoter versions\n"
cat output/alignments/*.txt | ./PromVerCount.awk > output/VariantCounts-100.csv

### get number of segregating sites for each promoter alignment
printf "Pulling out information about segregating sites\n"
python3 SegSites.py output/alignments/ 100
