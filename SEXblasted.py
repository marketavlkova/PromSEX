### Script which extracts promoter sequences according to blast search

### example to run: python3 SEXblasted.py <extracted K12 promoters - fasta file> <path to blast results> <path to env. strains' genomes in fasta> <number of which alignment length can be smaller than k12 promoter to be included>
### e.g.: python3 SEXblasted.py output/RDBex_Strong_eps.fasta output/blast_results/ input/marketa_genomes/ 100

import sys
import csv
import glob
from Bio.Seq import reverse_complement
from Bio import SeqIO, SeqFeature
from tqdm import tqdm           ### for nice progress bars

### main function
def main():

    ### save arguments
    k12_proms = sys.argv[1]
    blast_dir = sys.argv[2]
    strains_dir = sys.argv[3]
    diff = sys.argv[4]

    blast_results = glob.glob(blast_dir + '*.txt')
    # strains_fasta = glob.blog(strains_dir + '*.fsa')

    pbar = tqdm(blast_results)
    for strain in pbar:
        str = strain.split("/")[2]
        pbar.set_description("Processing %s" % str)
        with open(strain, 'r') as file:
            prom_out = open('output/blasted_promoters/' + '{}_-{}included.fasta'.format(str.split(".")[0], diff), 'w')
            for row in file:
                prom_idf = row.split("\t")[0]
                node_idf = row.split("\t")[1]
                align_length = row.split("\t")[3]
                sstart = int(row.split("\t")[8])
                send = int(row.split("\t")[9])
                evalue = row.split("\t")[10]
                with open(k12_proms, 'r') as proms:
                    for line in proms:
                        if line.startswith(">"):
                            prom_idp = line.split(">")[1].split("\n")[0]
                        else:
                            prom_length = len(line)
                            if all([prom_idf == prom_idp, int(align_length) >= prom_length-int(diff)]):
                                str_only = str.split("_")[1]
                                for fsa in SeqIO.parse(strains_dir + str_only + "/SC1_" + str_only + '.fsa', "fasta"):
                                    if node_idf in repr(fsa.id):
                                        prom_out.write(">" + fsa.id + " " + prom_idp + " alignment_length: " + align_length + " evalue: " + evalue + " strand: ")
                                        if sstart < send:
                                            prom_out.write("+\n")
                                            prom_out.write("%s\n" % fsa.seq[sstart:send])
                                        else:
                                            prom_out.write("-\n")
                                            prom_out.write("%s\n" % fsa.seq[send:sstart])


if __name__ == '__main__':
    main()
