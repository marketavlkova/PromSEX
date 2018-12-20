#!/usr/bin/awk -f
### extract just total diversity from Diversity-?00included.fasta file
### example to run: ./IDextract.awk <path to Indentity file> > <desired output file>
### e.g.: ./IDextract.awk output/Identities-100.txt > output/TotalIdent-100.csv

BEGIN{
  ### print column description
  printf "Promoter,Indentity,NumberOfIsolates\n";
  N=0;
}
{
  ### extract and print promoter ID
  if (($0 ~ /^#/) && ($3 ~ /^U00096.3_MG1655_genome/)) {
    split($3, P, "_");
    N=1;
  } else if ($0 ~ /^#/) {
    N+=1;
  }
  ### print total percent identity of the promoter among strains
  if (($1 ~ /^TOT/) && ($4 !~ /nan/)) {
    printf "%s,%s,%s\n", P[4], $4, N-1;
  }

}
