### Script which extracts promoter sequences according to blast search

### example to run: python3 SEXblasted.py <extracted K12 promoters - fasta file> <path to blast results> <path to env. strains' genomes in fasta> <number of which alignment length can be smaller than k12 promoter to be included>
### e.g.: python3 SEXblasted.py output/RDBex_Strong_eps.fasta output/blast_results/ input/marketa_genomes/ 100

import sys
import glob
from Bio.Seq import reverse_complement
from Bio import SeqIO, SeqFeature
from tqdm import tqdm           ### for nice progress bars

### function to determine number of lines in a file
def get_num_lines(file_path):
    query = open(file_path, 'r')
    lines = 0
    for line in query:
        lines += 1
    return lines

### main function
def main():

    ### save arguments
    k12_proms = sys.argv[1]
    blast_dir = sys.argv[2]
    strains_dir = sys.argv[3]
    diff = sys.argv[4]

    ### read all files with blast results
    blast_results = glob.glob(blast_dir + '*.txt')

    print("Extracting blasted promoters.")
    pbar = tqdm(blast_results)
    ### go throught all strain files with blast results
    for strain in pbar:
        str = strain.split("/")[2]
        pbar.set_description("Processing %s" % str)
        with open(strain, 'r') as file:
            ### make a separate file for each strain
            prom_out = open('output/blasted_promoters/' + '{}_-{}included.fasta'.format(str.split(".")[0], diff), 'w')
            ### loop through all line in strain's blast results
            for row in file:
                ### save important variables from blast results
                prom_idf = row.split("\t")[0]
                node_idf = row.split("\t")[1]
                align_length = row.split("\t")[3]
                sstart = int(row.split("\t")[8])
                send = int(row.split("\t")[9])
                evalue = row.split("\t")[10]
                ### now loop through all extracted k12 intergenic regions
                with open(k12_proms, 'r') as proms:
                    for line in proms:
                        ### save promoter name if line starts with ">"
                        if line.startswith(">"):
                            prom_idp = line.split(">")[1].split("\n")[0]
                        ### if not
                        else:
                            ### save length of k12 extracted intergenic region
                            prom_length = len(line)
                            ### continue only if promoter names are identical for k12 and tested strain
                            ### and at the same time the alighment length is not shorter then k12 intergenic region reduced by diff (value stated when running the script)
                            if all([prom_idf == prom_idp, int(align_length) >= prom_length-int(diff)]):
                                ### extract just the strain name
                                str_only = str.split("_")[1]
                                ### loop through all sequences of that particular strain in its fasta file
                                for fsa in SeqIO.parse(strains_dir + str_only + "/SC1_" + str_only + '.fsa', "fasta"):
                                    ### find the sequence where blast found a hit with the k12 intergenic region
                                    if node_idf in repr(fsa.id):
                                        ### and extract intergenic region of this particular strain
                                        prom_out.write(">" + fsa.id + "|" + prom_idp + "|alignment_length:_" + align_length + "|evalue:_" + evalue + "|strand:_")
                                        if sstart < send:
                                            prom_out.write("+\n%s\n" % fsa.seq[sstart-1:send])
                                        else:
                                            prom_out.write("-\n%s\n" % reverse_complement(fsa.seq[send-1:sstart]))

    ### read all files from previous step
    blast_extracted = glob.glob('output/blasted_promoters/' + '*{}included.fasta'.format(diff))

    print("Creating separated file for each promoter.")
    pbar = tqdm(SeqIO.parse(k12_proms, "fasta"), total=get_num_lines(k12_proms)//2)
    ### loop through all lines in file with extracted k12 intergenic regions
    for entry in pbar:
        ### save promoter name only
        promoter = entry.id.split("_")[3]
        pbar.set_description("Processing %s" % promoter)
        ### create a separate file for each promoter
        out_by_prom = open('output/blasted_promoters/by_promoter/' + '{}_-{}included.fasta'.format(promoter, diff), 'w')
        ### write k12 annotation and promoter sequence into it
        out_by_prom.write(">%s\n%s\n" % (entry.id, entry.seq))
        ### loop through all promoters filtered in previous step for all strains
        for file in blast_extracted:
            for line in SeqIO.parse(file, "fasta"):
                ### if there is such promoter present in that file, add it to the promoter output file
                if "_{}_".format(promoter) in repr(line.id):
                    out_by_prom.write(">%s|SC1_%s\n%s\n" % (line.id.split("\n")[0], file.split("_")[2], line.seq))


if __name__ == '__main__':
    main()
