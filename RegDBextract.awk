#!/usr/bin/awk -f
### extract defined information from RebulonDB file
### example to run: ./RegDBextract.awk <RegulonDB file> [columns to be extracted] > <output csv file>

BEGIN{

  COLS = ARGV[2];
  FS = "\t"

  if (COLS ~ /,/) {
    NC = split(COLS, C, ",");
  } else {
    C[1] = COLS;
    NC = 1;
  }

}
!/^#/{

  if (NC == 1) {
    printf "%s\n", $C[NC];
  } else {
    for (I = 1; I < NC; I++) {
      printf "%s,", $C[I];
    }
    printf "%s\n", $C[NC];
  }

}
