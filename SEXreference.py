### Script which searches and extracts promoter sequences
### from K12 genome according to a given RegulonDB promoter
### database

### example to run: python3 SEXreference.py <path to filtered RegulonDB csv file> <path to genome GenBank file> [confidence level - optional]
### e.g.: python3 SEXreference.py output/RDBex.csv input/MG1655_genome.gb Strong

import sys
import csv
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
    feats_out = open('output/' + '{}_{}_promoter_annotations.csv'.format(db_name, confidence), 'w')
    int_regs_out = open('output/' + '{}_{}_intergenic_regions.csv'.format(db_name, confidence), 'w')
    int_reg_seqs_out = open('output/' + '{}_{}_extracted_promoter_sequences.fasta'.format(db_name, confidence), 'w')

    ### load reference genome
    for genome_k12 in SeqIO.parse(reference, "genbank"):
        print(repr(genome_k12.seq))

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

    ### search borders of intergenic regions using GenBank annotations
    print('\n' + "Identifying intergenic regions. Promoters with TSS located in CDS or similar excluded:")
    int_regs = {}
    int_regs_out.write("Promoter,TSS,Intergenic_Start,Intergenic_End,Intergenic_length,strand" + '\n')
    pbar = tqdm(feats_fin)
    for id in pbar:
        if len(feats_fin[id]) == 1:
            pbar.set_description("Processing %s" % id)
            feat = feats_fin[id][0]
            strand = filtered_dict[id][0]
            tss = filtered_dict[id][1]
            codon = None
            int_reg_end = None
            int_reg_len = None
            end1 = None
            end2 = None
            end3 = None
            end4 = None
            ### shorten searching range towards 3' end of + strand
            for nt in range(tss, tss+1501, 10):
                for feature in genome_k12.features:
                    ### correction for intergenic regions spanning through ori in circular chromosome
                    if nt > len(genome_k12):
                        nt0 = nt - len(genome_k12)
                        if all([nt0 in feature, feature.type not in feat]):
                            if nt0-11 < 0:
                                end1 = abs(len(genome_k12) + nt0-11)
                            else:
                                end1 = nt0-11
                            end2 =nt0+1
                            break
                    if all([nt in feature, feature.type not in feat]):
                        end2 = nt+1
                        end1 = nt-11
                        break
                if isinstance(end2, int):
                    break
            ### shorten searching range towards 5' end of + strand
            for bp in range(tss, tss-1501, -10):
                for feature in genome_k12.features:
                    ### correction for intergenic regions spanning through ori in circular chromosome
                    if bp < 0:
                        bp0 = len(genome_k12) + bp
                        if all([bp0 in feature, feature.type not in feat]):
                            if bp0+11 > len(genome_k12):
                                end4 = abs(bp0+11 - len(genome_k12))
                            else:
                                end4 = bp0+11
                            end3 = bp0-1
                            break
                    if all([bp in feature, feature.type not in feat]):
                        end4 = bp+11
                        end3 = bp-1
                        break
                if isinstance(end4, int):
                    break
            ### find exact 3' position of intergenic region - considering '+ strand'
            if isinstance(end2, int):
                for nt2 in range(end1, end2, 1):
                    for feature in genome_k12.features:
                        if all([nt2 in feature, feature.type not in feat]):
                            if 'forward' in strand:
                                codon = nt2+1
                            else:
                                int_reg_end = nt2+1
                            break
                    if any([isinstance(codon, int), isinstance(int_reg_end, int)]):
                        break
            ### find exact 5' position of intergenic region - considetin '+ strand'
            if isinstance(end4, int):
                for bp2 in range(end4, end3, -1):
                    for feature in genome_k12.features:
                        if all([bp2 in feature, feature.type not in feat]):
                            if 'forward' in strand:
                                int_reg_end = bp2+1
                            else:
                                codon = bp2+1
                            break
                    if all([isinstance(codon, int), isinstance(int_reg_end, int)]):
                        break
            ### get length of intergenic region
            if all([isinstance(codon, int), isinstance(int_reg_end, int)]):
                int_reg_len = abs(codon-int_reg_end)
                ### correction for intergenic regions spanning through ori in circular chromosome
                if int_reg_len > 3000:
                    if codon < int_reg_end:
                        int_reg_len = abs(int_reg_end-len(genome_k12)+codon)
                    else:
                        int_reg_len = abs(codon-len(genome_k12)+int_reg_end)
            ### store positions in dictionary
            int_regs[id] = [codon, int_reg_end, strand]
            ### writting output file
            entry = [id, str(tss), str(codon), str(int_reg_end), str(int_reg_len), strand, '\n']
            row = ','.join(entry)
            int_regs_out.write(row)

    ### extract intergenic regions with extra 100bp upstream and downstream
    print('\n' + "Extracting intergenic regions from", reference, ":")
    pbar = tqdm(int_regs)
    for key in pbar:
        pbar.set_description("Processing %s" % id)
        down = int_regs[key][0]
        up = int_regs[key][1]
        strand = int_regs[key][2]
        int_reg_seqs_out.write(">" + genome_k12.id + "_" + reference.split("/")[1].split(".")[0] + "_" + key + "_" + strand + "\n")
        ### correction for intergenic regions spanning through ori in circular chromosome
        if abs(down-up) > 3000:
            if down < up:
                int_reg_seqs_out.write("%s" % genome_k12.seq[up-100:len(genome_k12.seq)+1])
                int_reg_seqs_out.write("%s\n" % genome_k12.seq[0:down+101])
            else:
                int_reg_seqs_out.write("%s" % genome_k12.seq[down-100:len(genome_k12.seq)+1])
                int_reg_seqs_out.write("%s\n" % genome_k12.seq[0:up+101])
        ### all other intergenic regions
        else:
            if down < up:
                int_reg_seqs_out.write("%s\n" % genome_k12.seq[down-100:up+101])
            else:
                int_reg_seqs_out.write("%s\n" % genome_k12.seq[up-100:down+101])

if __name__ == '__main__':
    main()
