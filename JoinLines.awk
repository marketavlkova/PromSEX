#!/usr/bin/awk -f
### joint lines with the same ID
### example to run: ./JoinLines.awk <input file> > <output file>
### e.g.: ./JoinLines.awk output/alignments/alkAp_-100.included.aln > output/alkAp_-100included.txt

{
  ID=$1;
  $1="";
  SEQ[ID]=SEQ[ID]$0;
}
END{
  printf "?\n";
  for (ID in SEQ) {
    printf "%s:%s\n", ID, SEQ[ID];
  }
  printf "!\n";
}
