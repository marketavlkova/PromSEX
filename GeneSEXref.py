### Script which searches and extracts promoter sequences
### from K12 genome according to a given RegulonDB promoter
### database

### example to run: python3 GeneSEXref.py <path to filtered RegulonDB csv file> <path to genome GenBank file> [confidence level - optional]
### e.g.: python3 GeneSEXref.py output/RDBex.csv input/MG1655_genome.gb Strong

import sys
import csv
import os
from Bio.Seq import reverse_complement
from Bio import SeqIO, SeqFeature
from tqdm import tqdm           ### for nice progress bars

### function to flatten list of lists
def flatten(listOfLists):
    flattened = []
    for item in listOfLists:
        if isinstance(item, str):
            flattened.append(item)
        else:
            flattened.extend(flatten(item))
    return flattened

### main function
def main():

    ### save & read arguments
    reg_db = sys.argv[1]
    reference = sys.argv[2]
    if len(sys.argv) > 3:
        confidence = sys.argv[3]
    else:
        confidence = None
        print("Confidence level not stated - taking all proposed promoters")

    ### define output files
    db_name = reg_db.split("/")[1].split(".")[0]
    genes_out = open('output/' + '{}_{}_gene_regions.csv'.format(db_name, confidence), 'w')
    gene_seqs_out = open('output/' + '{}_{}_extracted_genes.fasta'.format(db_name, confidence), 'w')

    ### load reference genome
    for genome_k12 in SeqIO.parse(reference, "genbank"):
        print('')

    ### filter promoters with desired confidence level
    print('\n' + "Filtering promoters with", confidence, "confidence level.")
    filtered_dict = {}
    with open(reg_db, 'r') as proms:
        for line in proms:
            if line.split(",")[4] == confidence + '\n':
                pos = int(line.split(",")[2])
                ### remove promoters without defined TSS position
                if pos > 0:
                    key = str(line.split(",")[0])                                           ### get promoter name
                    filtered_dict[key] = [line.split(",")[1], pos, line.split(",")[3]]      ### store strand orientation, TSS position and sequence including this TSS

    ### if file with sequence annotations of tss exists, read it
    if os.path.isfile('output/' + '{}_{}_promoter_annotations.csv'.format(db_name, confidence)):
        feats_fin = {}
        with open('output/' + '{}_{}_promoter_annotations.csv'.format(db_name, confidence)) as in_file:
            in_line = in_file.readline()
            while in_line:
                # print(in_line)
                in_line = in_file.readline()
                id = in_line.split(",")[0]
                input = []
                for i in range(1, (len(in_line.split(",")) - 1)):
                    input += [''.join(in_line.split(",")[i])]
                feats_fin[id] = input
        del feats_fin['']
    ### if not make it
    else:
        feats_out = open('output/' + '{}_{}_promoter_annotations.csv'.format(db_name, confidence), 'w')

        ### get annotations for all filtered promoter TSS
        print('\n' + "Searching sequence annotations of TSS in filtered promoters:")
        feats = {}
        pbar = tqdm(filtered_dict)
        for id in pbar:
            pbar.set_description("Processing %s" % id)
            for feature in genome_k12.features:
                if filtered_dict[id][1] in feature:
                    if id not in feats:
                        feats[id] = feature.type
                    else:
                        feats[id] = [feats[id], feature.type]

        ### iterate through the dictionary to flatten the output file
        feats_clear = {}
        for id in feats:
            feats_clear[id] = flatten(feats.get(id))

        ### clean after flattening & save into a file
        feats_fin = {}
        repaired = []
        for id in feats_clear:
            if len(feats_clear[id][0]) == 1:
                item = feats_clear.get(id)
                repaired = [''.join(item)]
            else:
                repaired = feats_clear[id]
            feats_fin[id] = repaired
            ### writting output file
            key = [id]
            key += repaired
            key.append('\n')
            row = ','.join(key)
            feats_out.write(row)

    ### search borders of genes using GenBank annotations
    print('\n' + "Identifying gene regions. Genes having TSS located within CDS or similar excluded:")
    gene_regs = {}
    genes_out.write("Gene,TSS,Gene_Start,Gene_End,Gene_length,strand" + '\n')
    pbar = tqdm(feats_fin)
    for id in pbar:
        if len(feats_fin[id]) == 1:
            pbar.set_description("Processing %s" % id)
            # feat = feats_fin[id][0]
            strand = filtered_dict[id][0]
            tss = filtered_dict[id][1]
            codon = None
            gene_end = None
            gene_len = None
            end1 = None
            end2 = None
            end3 = None
            end4 = None
            if 'forward' in strand:
                ### shorten searching range towards 3' end of + strand
                for nt in range(tss, tss+1501, 5):
                    for feature in genome_k12.features:
                        ### correction for genes spanning through ori in circular chromosome
                        if nt > len(genome_k12):
                            nt0 = nt - len(genome_k12)
                            if all([nt0 in feature, 'gene' in feature.type]):
                                if nt0-6 < 0:
                                    end1 = abs(len(genome_k12) + nt0-6)
                                else:
                                    end1 = nt0-6
                                end2 =nt0+1
                                break
                        if all([nt in feature, 'gene' in feature.type]):
                            end2 = nt+1
                            end1 = nt-6
                            break
                    if isinstance(end2, int):
                        break
                ### find exact 3' position of intergenic region - considering '+ strand'
                if isinstance(end2, int):
                    for nt2 in range(end1, end2, 1):
                        for feature in genome_k12.features:
                            if all([nt2 in feature, 'gene' in feature.type]):
                                codon = nt2
                                break
                        if isinstance(codon, int):
                            break
                ### shorten searching range towards 5' end of + strand
                point = None
                for feature in genome_k12.features:
                    if all(['gene' in feature.type, codon in feature.location]):
                        for bp in range(codon, codon+15001, 5):
                            ### correction for genes spanning through ori in circular chromosome
                            if bp > len(genome_k12):
                                bp0 = bp - len(genome_k12)
                                if bp0 in feature.location:
                                    point = bp0
                                if point-1 < 0:
                                    end3 = abs(len(genome_k12) + point-1)
                                else:
                                    end3 = point-1
                                end4 = point+6
                            if bp in feature.location:
                                point = bp
                            if isinstance(point, int):
                                end4 = point+6
                                end3 = point-1
                        ### find exact 5' position of intergenic region - considetin '+ strand'
                        point = None
                        if isinstance(end4, int):
                            for bp2 in range(end3, end4, 1):
                                if bp2 in feature.location:
                                    point = bp2
                                gene_end = point
            else:
                ### shorten searching range towards 3' end of - strand
                for nt in range(tss, tss-1501, -5):
                    for feature in genome_k12.features:
                        ### correction for genes spanning through ori in circular chromosome
                        if nt < 0:
                            nt0 = len(genome_k12) + nt
                            if all([nt0 in feature, 'gene' in feature.type]):
                                if nt0+6 > len(genome_k12):
                                    end1 = abs(nt0+6 - len(genome_k12))
                                else:
                                    end1 = nt0+6
                                end2 =nt0-1
                                break
                        if all([nt in feature, 'gene' in feature.type]):
                            end2 = nt-1
                            end1 = nt+6
                            break
                    if isinstance(end2, int):
                        break
                ### find exact 3' position of intergenic region - considering '- strand'
                if isinstance(end2, int):
                    for nt2 in range(end1, end2, -1):
                        for feature in genome_k12.features:
                            if all([nt2 in feature, 'gene' in feature.type]):
                                codon = nt2
                                break
                        if isinstance(codon, int):
                            break
                ### shorten searching range towards 5' end of - strand
                point = None
                for feature in genome_k12.features:
                    if all(['gene' in feature.type, codon-1 in feature.location]):
                        for bp in range(codon-1, codon-15001, -5):
                            ### correction for genes spanning through ori in circular chromosome
                            if bp < 0:
                                bp0 = len(genome_k12) + bp
                                if bp0 in feature.location:
                                    point = bp0
                                if point+1 > len(genome_k12):
                                    end3 = abs(point+1 - len(genome_k12))
                                else:
                                    end3 = point+1
                                end4 = point-6
                            if bp in feature.location:
                                point = bp
                            if isinstance(point, int):
                                end4 = point-6
                                end3 = point+1
                        ### find exact 5' position of intergenic region - considetin '- strand'
                        point = None
                        if isinstance(end4, int):
                            for bp2 in range(end3, end4, -1):
                                if bp2 in feature.location:
                                    point = bp2
                                gene_end = point

            ### get length of the gene
            if all([isinstance(codon, int), isinstance(gene_end, int)]):
                gene_len = abs(codon-gene_end)
                ### correction for genes spanning through ori in circular chromosome
                if gene_len > 10000:
                    if codon < gene_end:
                        gene_len = abs(gene_end-len(genome_k12)+codon)
                    else:
                        gene_len = abs(codon-len(genome_k12)+gene_end)
            ### store positions in dictionary
            gene_regs[id] = [codon, gene_end, strand]
            ### writting output file
            entry = [id, str(tss), str(codon), str(gene_end), str(gene_len), strand, '\n']
            row = ','.join(entry)
            genes_out.write(row)

    ## extract genes
    print('\n' + "Extracting genes from", reference, ":")
    pbar = tqdm(gene_regs)
    for key in pbar:
        pbar.set_description("Processing %s" % key)
        down = gene_regs[key][0]
        up = gene_regs[key][1]
        strand = gene_regs[key][2]
        gene_seqs_out.write(">" + genome_k12.id + "_" + reference.split("/")[1].split(".")[0] + "_" + key + "_" + strand + "\n")
        ### correction for genes spanning through ori in circular chromosome
        if abs(down-up) > 3000:
            if down < up:
                gene_seqs_out.write("%s" % genome_k12.seq[up:len(genome_k12.seq)+1])
                gene_seqs_out.write("%s\n" % genome_k12.seq[0:down+1])
            else:
                gene_seqs_out.write("%s" % genome_k12.seq[down:len(genome_k12.seq)+1])
                gene_seqs_out.write("%s\n" % genome_k12.seq[0:up+1])
        ### all other genes
        else:
            if down < up:
                gene_seqs_out.write("%s\n" % genome_k12.seq[down:up+1])
            else:
                gene_seqs_out.write("%s\n" % genome_k12.seq[up:down+1])

if __name__ == '__main__':
    main()
