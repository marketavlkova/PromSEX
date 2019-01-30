### Script which counts number of segregation sites for each promoter

### example to run: python3 SegSites.py <path to dir with .aln files> <number of which alignment length can be smaller than k12 promoter to be included>
### e.g.: python3 SeqSites.py output/alignments/ 100

import sys
import glob
from Bio import AlignIO
from Bio.Align import AlignInfo

def main():

    ### save arguments
    in_dir = sys.argv[1]
    diff = sys.argv[2]

    ### read all files with blast results
    in_files = glob.glob(in_dir + '*' + diff + '*.aln')

    ### set output file
    out_file = open('output/' + 'SegSites_-{}included.csv'.format(diff), 'w')
    out_file.write('Promoter,#SegSites,#Strains,AlignLength\n')

    ### loop through all files
    for file in in_files:
        ### load promoter alignment
        align = AlignIO.read(open(file), "clustal")
        promoter = file.split("/")[2].split("_")[0]
        len = align.get_alignment_length()
        nstr = 0
        ### count number of strains having this promoter
        for seq in align:
            nstr += 1

        ### extract general summary about the alignment
        summ = AlignInfo.SummaryInfo(align)
        cons = summ.dumb_consensus()
        ### generate Position Specific Score Matrix
        pssm = summ.pos_specific_score_matrix(cons, chars_to_ignore = [ 'N' ])
        s = 0
        ### count number of segregation sites in this alignment (excluding indels)
        for nt in range(0, len):
            base = align[0][nt]
            if (pssm[nt][base] != nstr):
                if ('-' in pssm[nt]):
                    if (pssm[nt]['-'] > 0):
                        for i in pssm[nt]:
                            if all([i != '-', pssm[nt][i] > 0]):
                                if ((pssm[nt][i] + pssm[nt]['-']) != nstr):
                                    s += 1
                                    break
                    else:
                        s += 1
                else:
                    s += 1

        ### add acquired values to the output file
        out_file.write("%s,%s,%s,%s\n" % (promoter, s, nstr, len))

if __name__ == '__main__':
    main()
