#!/usr/bin/awk -f
### count how many versions of each promoter exist among our strains
### example to run: cat <path to input files> | ./JPromVerCount.awk > <path to output file>
### cat output/alignments/*100*.txt | ./PromVerCount.awk > output/VariantCounts-100.csv

BEGIN{
  ### print column description
  printf "Promoter,NumberOfVariants,NumberOfIsolates\n";
  ### set field separator as ":"
  FS=":";
}
{
  ### find beginning of each file in line and set variables
  if ($0 ~ /^?/) {
    N=0;
    VER=1;
  }
  ### if analyzing K12 strain, extract promoter name, sequence & strain ID
  if ($0 ~ /^U00096.3/) {
    split($1, P, "_");
    ID=P[1];
    SEQ[ID]=$2;
  ### in all other strains extract just promoter sequence and strain ID
  } else if ($0 ~ /^SC1/) {
    split($1, Q, "|");
    ID=Q[1];
    SEQ[ID]=$2;
  }
  ### at the end of each file
  if ($0 ~ /^!/) {
    ### loop through all strains twice
    for (A in SEQ) {
      for (B in SEQ) {
        ### when comparing to different strains, which have distinct promoter sequences
        if ((A != B) && (SEQ[A] != SEQ[B])) {
          ### if IDENT array is not defined yet
          ### assign the promoter sequence of first strain to it
          if (length(IDENT) == 0) {
            IDENT[VER]=SEQ[A];
          }
          HIT=0;
          ### loop though all values of IDENT array and
          ### if second sequence is already asigned to one of them
          ### change HIT variable to 1
          for (I=1; I<=VER; ++I) {
            if (IDENT[I] == SEQ[B]) {
              HIT=1;
            }
          }
          ### if the second sequence is not among those in IDENT array
          ### add it there
          if (HIT == 0) {
            VER+=1;
            IDENT[VER]=SEQ[B];
          }
        }
      }
      ### count number of strains having that particular promoter
      N+=1;
    }
    ### print values into output file
    printf "%s,%s,%s\n", P[4], VER, N;
    ### delete variables before going through another input file in line
    delete SEQ;
    delete IDENT;
  }
}
