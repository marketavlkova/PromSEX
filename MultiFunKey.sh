#!/bin/bash
### extract IDs and their descriptions of ontology functions from MultiFun
### and also join BC codes with appropriate common gene names from MultiFun
### example to run: ./MultiFunKey.sh <path to file with FultiFun IDs> <path to output file with the key> <path to output file with join BC codes and common gene names>
### ./MultiFunKey.sh 'output/AllGenes_IDs?.txt' output/MutliFunKey.csv output/AllGenesMFinput.csv

### loop through all file matching 1st argument
for fileID in $1
do
  ### define file with ID descriptions
  head="${fileID%_*}_Commons"
  tail="${fileID#*IDs}"
  fileCOM="${head}${tail}"
  ### define temporary output file
  out="output/out${tail}"

  ### print the 2 input files next to each other and give it to awk
  paste $fileID $fileCOM | awk 'BEGIN{
    ### print column description
    printf "ID,function\n";
    ### set field separator as tabulator
    FS="\t";
  }
  ### exclude the header for processing
  !/Ontology/{
    ### print ID name and its description into temporary output file
    printf "%s\t%s\n", $1, $3;
  }' > $out

  ### extract the number from input file name
  number="${tail%.*}"
  ### join BC codes with common gene names if working with the first data files
  if [ "${number}" == "0" ]
  then
    paste $fileID $fileCOM | awk 'BEGIN{
      ### print column description
      printf "ID,Gene\n";
      ### set field separator as tabulator
      FS="\t";
    }
    ### exclude the header for processing
    !/Ontology/{
      ### print ID name and its description into temporary output file
      printf "%s,%s\n", $1, $4;
    }' > $3
  fi
done

### print all temporary input files and give this output to awk
cat output/out?.txt | awk 'BEGIN{
  ### print column description
  printf "ID,function\n";
  ### set field separator as tabulator
  FS="\t";
  ### set array counter
  N=1;
}
### go only through ID starting with BC
### (other ones are not really IDs)
/^BC-/{
  ### processing the first line
  if (NR == 2) {
    ### assign ID and its description to arrays
    ID[N]=$1;
    FUN[N]=$2;
    ### increase array counter
    N+=1;
  ### if processing other lines
  } else {
    ### set variable to check for presence of ID in the array
    CHECK=0;
    ### loop through the array with ID
    for (I=1; I<N; I++) {
      ### if the ID currently processing is present in the array
      if (ID[I] == $1) {
        ### set variable checking its presence to 1
        CHECK+=1;
      }
    }
    ### if the ID currently processing is not in the array
    if (CHECK == 0) {
      ### add it there together with its description
      ID[N]=$1;
      FUN[N]=$2;
      ### and increase the array counter
      N+=1;
    }
  }
}
END{
  ### print both arrays as 2 separate columns into final output file
  for (I=1; I<N; I++) {
    printf "%s,%s\n", ID[I], FUN[I];
  }
}' > $2

### remove temporary output files
rm output/out?.txt
