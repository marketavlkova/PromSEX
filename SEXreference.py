### Script which searches and extracts promoter sequences
### from K12 genome according to a given RegulonDB promoter
### database

### example to run: python3 SEXreference.py <RegulonDB csv file> <genome GenBank file> [confidence level - optional]

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
        print("Confidence level not stated - taking all proposed promoters")

    ### define output files
    db_name = reg_db.split(".")[0]
    feats_out = open('output/' + '{}_{}_promoter_annotations.csv'.format(db_name, confidence), 'w')
    cds_start_out = open('output/' + '{}_{}_start_codons.csv'.format(db_name, confidence), 'w')

    ### load reference genome
    for genome_k12 in SeqIO.parse(reference, "genbank"):
        print(repr(genome_k12.seq))

    ### filter promoters with desired confidence level
    print('\n' + "Filtering promoters with", confidence, "confidence level.:")
    filtered_dict = {}
    for tb in tqdm(range(sum(1 for line in open(reg_db))//100)):           ### this is progress bar
        with open(reg_db, 'r') as proms:
            for line in proms:
                if line.split(",")[4] == confidence + '\n':
                    pos = int(line.split(",")[2])
                    ### remove promoters without defined TSS position
                    if pos > 0:
                        key = str(line.split(",")[0])                                           ### get promoter name
                        filtered_dict[key] = [line.split(",")[1], pos, line.split(",")[3]]      ### store strand orientation, TSS position and sequence including this TSS
                        ### get position 500nt upstream or downstream of TSS for start codon search
                        if filtered_dict[key][0] == 'reverse':
    #                        tss_codon = reverse_complement(genome_k12.seq[pos-101:pos])
                            cds_search = pos-501
                        else:
    #                        tss_codon = genome_k12.seq[pos-1:pos+100]
                            cds_search = pos+500
                        filtered_dict[key].append(cds_search)

    ### get annotations for all filtered promoter TSS
    print('\n' + "Searching sequence annotations of TSS in filtered promoters.:")
    for tb in tqdm(range(len(filtered_dict)//100)):           ### this is progress bar
        feats = {}
        for feature in genome_k12.features:
            for id in filtered_dict:
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

    ### find start codons using GenBank annotations
    print('\n' + "Looking for start codons. TSS located in CDS or similar excluded.:")
    cds_hits = {}
    cds_start_out.write("Promoter,TSS,CDS_start,strand" + '\n')
    for id in feats_fin:
        if len(feats_fin[id]) == 1:
            feat = feats_fin[id][0]
            tss = filtered_dict[id][1]
            stop = filtered_dict[id][3]
            codon = 0
            print("Processing:", id)
            for tb in tqdm(range(5)):           ### this is progress bar
                for feature in genome_k12.features:
                    if tss < stop:
                        strand = '+'
                        for nt in range(tss, stop+1, 1):
                            if all([nt in feature, feature.type not in feat]):
                                codon = nt
                                codon += 1
                                break
                            if all([nt == stop, codon == 0]):
                                codon = None
                    else:
                        strand = '-'
                        for nt in range(tss, stop-1, -1):
                            if all([nt in feature, feature.type not in feat]):
                                codon = nt
                                codon += 1
                                break
                            if all([nt == stop, codon == 0]):
                                codon = None
            cds_hits[id] = codon
            ### writting output file
            entry = [id, str(tss), str(codon), strand]
            entry.append('\n')
            row = ','.join(entry)
            cds_start_out.write(row)

if __name__ == '__main__':
    main()
