#!/bin/bash
### script blasting intergenic sequences obtained by SEXreference.py
### againts databases of other genomes - which is can also create if needed

### example to run: bash Blast.sh

### ask whether blast database should be made
printf "Do you wish to make blast databases? (y/n) \n"
read ANS
if [ "$ANS" != "${ANS#[Yy]}" ]
then
  DB_PATH="input/Database"
  printf "Enter path to fasta file/s the database/s should be made from: \n"
  read FAS
  printf "Enter type of database to be made: (nucl/prot) \n"
  read DBT
  ### loop through all fasta (fsa) genomes provided and make database from each of them
  for FASTA in $FAS/*/*.fsa
  do
    PRE_DB="${FASTA#*/*/*/}"
    FIN_DB="${PRE_DB%.*}"
    printf "Making $OUT_DB blast database: \n"
    makeblastdb -in $FASTA -dbtype $DBT -out $DB_PATH/$FIN_DB
  done
else
  printf "Enter path to your database/s: \n"
  read DB_PATH
fi

### prompt for query file (fasta output from SEXreference.py)
printf "Enter your query sequence file (including the path to it): \n"
read QUERY_FILE

### blasting
for FILE in $DB_PATH/*.nhr
do
  IN_DB="${FILE#*/*/}"
  OUT_FILE="${IN_DB%.*}_blasted.txt"
  printf "Blasting $IN_DB \n"
  blastn -db $DB_PATH/"${IN_DB%.*}" -query $QUERY_FILE -strand both -task megablast -evalue 1e-10 -out output/blast_results/$OUT_FILE -outfmt "6 std" -num_threads 4
done
