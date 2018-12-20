#!/usr/bin/awk -f
### extract content of desired columns in RebulonDB file
### example to run: ./RegDBextract.awk -v COLS=[columns to be extracted] <RegulonDB file> > <output csv file>
### e.g.: ./RegDBextract.awk -v COLS=2,3,4,6,8 input/MG1655_s70_promoters_RegDB.txt > output/RDBex.csv

BEGIN{
  ### set tabulator as a field separator
  FS = "\t"
  ### check which columns should be extracted
  if (COLS ~ /,/) {
    ### NC = number of columns to be extracted; C = array of column numbers to be extracted
    NC = split(COLS, C, ",");
  } else {
    ### if only one column is specified
    C[1] = COLS;
    NC = 1;
  }

}
!/^#/{
  ### if only one column is specified to be extracted,
  ### print all its fields in a column
  if (NC == 1) {
    printf "%s\n", $C[NC];
  } else {
    ### otherwise loop through all desired columns in each line
    for (I = 1; I < NC; I++) {
      printf "%s,", $C[I];
    }
    printf "%s\n", $C[NC];
  }

}
