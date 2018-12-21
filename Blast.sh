#!/bin/bash
### script blasting intergenic sequences obtained by SEXreference.py
### againts databases of other genomes - which is can also create if needed

### example to run: ./Blast.sh [should be database created? Yes/No] <path to fasta file with K12 intergenic regions>

### set path to blast database
DB_PATH="input/Database"
### set path to query file
QUERY_FILE=$2

### check whether blast database should be made
if [ "$1" != "${1#[Yy]}" ]
then
  ### loop through all fasta (fsa) genomes provided and make database from each of them
  for FASTA in input/marketa_genomes/*/*.fsa
  do
    PRE_DB="${FASTA#*/*/*/}"
    FIN_DB="${PRE_DB%.*}"
    makeblastdb -in $FASTA -dbtype nucl -out $DB_PATH/$FIN_DB
  done
fi

### blasting
for FILE in $DB_PATH/*.nhr
do
  IN_DB="${FILE#*/*/}"
  OUT_FILE="${IN_DB%.*}_blasted.txt"
  printf "Blasting $IN_DB \n"
  blastn -db $DB_PATH/"${IN_DB%.*}" -query $QUERY_FILE -strand both -task megablast -evalue 1e-10 -out output/blast_results/$OUT_FILE -outfmt "6 std" -num_threads 4
done
