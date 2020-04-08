#!/usr/bin/awk -f
### extract gene names from promoter names for MultiFun analysis
### example to run: ./GetGenesFromProms.awk <path to csv file with promoter name in 1st column> > <desired output file>
### e.g.: ./GetGenesFromProms.awk output/AllBoth.csv > output/AllGenes.csv

BEGIN{
  ### set field separator as ","
  FS=",";
  ### print column description
  printf "Gene\n";
}
{
  if (length($1) > 2) {
    ### split given promoter name into single letters
    N=split($1, A, "");
    ### remove the last "p" or "p[0-9]" character
    ### while saving & printing the rest
    for (I=2; I < N; I++) {
      if ((I <= 4) || ((A[I] != "p") && (A[I] !~ /^[0-9]+$/))) {
        printf "%s", A[I];
      }
    }
    printf "\n";
  }
}
