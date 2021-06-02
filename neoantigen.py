import os
import sys
import argparse
import pandas as pd
import re
from Bio import SeqIO


def parse_Args():
    parser = argparse.ArgumentParser(
        description='Process neoantigens', add_help=True)
    parser.add_argument('--MAF',  type=str,
                        help='MAF file which has mutation annotations')
    parser.add_argument('--config',  type=str,
                        help='config file which has all system requirements')
    parser.add_argument('--cdna',  type=str,
                        help='cDNA file which has nucleotide sequences of genes')
    parser.add_argument('--cds',  type=str,
                        help='coding regions for the cDNA file')
    parser.add_argument('--out',  type=str,
                        help='Directory where all the input will be written')
    parser.add_argument('--sample',  type=str,
                        help='Sample name of the run')

    args = parser.parse_args()

    # Intitialize arguments
    maf_file = args.MAF
    config_file = args.config
    cdna_file = args.cdna
    cds_file = args.cds
    out_dir = args.out
    sample = args.sample
    return (maf_file, cdna_file, cds_file, out_dir, sample)


def runNetMHC_pan():
	cmd = "/gpfs/data/krogsgaardlab/Pamela/krogsgaardlab/Neoantigen_pipeline/netMHCpan-4.1/netMHCpan -s -BA  -a HLA-A03:01,HLA-A29:02:01,HLA-B35:02:01,HLA-B58:01:01,HLA-C07:18,HLA-C04:01:01  -f TNBC_72/TNBC_72.mutated_sequences.fa -l 9,10,11 -inptype 0  -xls  -xlsfile TNBC_72/TNBC_72.netmhcpan.output.xls  > TNBC_72/TNBC_72.netmhcpan.output.txt"


def make_dir(out_dir):
    try:
        os.makedir(out_dir)
    except:
        print('Directory already exists')


def match_hgsc(hgvsc):

    position, ref_allele, alt_allele, sequence, hgvsc_type = [
        0, '', '', '', '']
    print(hgvsc)

    # matches substitution
    if re.match('^c\.(\d+).*([ATCG]+)>([ATCG]+)$', hgvsc):
        position, ref_allele, alt_allele = re.match(
            r'^c\.(\d+).*(\w+)>(\w+)', hgvsc).groups()

    # Example Frame shift insertion c.152_153insAGCTG
        # Also something like this c.167_167+1insTGCTGACAATACTT
    elif re.match('^c\.(\d+).*.(ins).([ATCG]+)$', hgvsc):
        position,  hgvsc_type, sequence = re.match(
            r'^c\.(\d+).*.(ins)([ATCG]*)$', hgvsc).groups()

    # In_frame_del and Frame shift deletion have the same kind of format
    elif re.match('^c\.(\d+).*(del)$', hgvsc):
        position, hgvsc_type, sequence = re.match(
            '^c\.(\d+)._.(\d+).*(del)$', hgvsc).groups()
        entries = re.split("_", hgvsc)
        hgvsc_type = "del"

        #  Frame shift Insertion specifically/Should not pick up/Have to think about this. its not right
    elif re.match('^c\.(\d+).*(dup)$', hgvsc):
        position, hgvsc_type = re.match(r'^c\.(\d+).*(dup)$', hgvsc).groups()
        #hgvsc_type = "novel"

    else:
        sys.exit('Error: Does not match to the format: ' + hgvsc)

    position = int(position) - 1

    if hgvsc_type == 'del':
        ref_allele = sequence
    elif hgvsc_type == 'ins':
        alt_allele = sequence
    elif hgvsc_type == 'dup':
        alt_allele = sequence

    return (position, ref_allele, alt_allele)


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])


def cds_to_aa(cds):
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_', 'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
        'NAC': '_', 'NAT': '_', 'NAA': '_', 'NAG': '_', 'NGC': '_', 'NGT': '_', 'NGA': '_', 'NGG': '_',
        'NTC': '_', 'NTA': '_', 'NTT': '_', 'NTG': '_', 'NCC': '_', 'NCA': '_', 'NCT': '_', 'NCG': '_'

    }
    protein = ''
    for i in range(0, len(cds), 3):

        codon = cds[i:i + 3]

        if len(codon) != 3:
            break

        if codon_table[codon] == '_':  # stop codon reached
            break
        protein += codon_table[codon]

    return (protein)


def variant_approve(variant):
    variant_type = ['SNP', 'DEL', 'INS', 'DUP']
    if (variant in variant_type):
        return 1
    else:
        return 0


def type_non_syn(type):
    non_syn_types = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del',
                     'In_Frame_Ins', 'Missense_Mutation', 'Nonstop_Mutation']
    if (type in non_syn_types):
        return 1
    else:
        return 0


def read_dataframe(filename):
    data = pd.read_table(filename, low_memory=False, header=1)
    return (data)


def load_cds_fasta(fasta_file):
    seqs = dict()
    if fasta_file[-3:len(fasta_file)] == '.gz':
        lines = gzip.open(fasta_file, 'rb').readlines()
    else:
        lines = open(fasta_file).readlines()
    idx = 0
    while idx < len(lines):
        line = lines[idx]

        # Extract the line which has the header and starts with the transcript ID
        begin = re.search('^>ENST\d+', line)
        transcript_id = ''
        if not begin:
            sys.exit('Error parsing fasta file: ' +
                     fasta_file + ' at line: ' + line)
        else:
            transcript_id = begin.group()

        idx = idx + 1
        seq_str = ''

        # Extract the rest of the lines into a dictionary
        while idx < len(lines) and not re.match('^>ENST', lines[idx]):
            seq_str = seq_str + lines[idx].strip()
            idx = idx + 1
            new_trans_id = transcript_id.replace(">", "")
            seqs[new_trans_id] = seq_str

    return (seqs)


def main():
    (maf_file, cdna, cds, out_dir, sample) = parse_Args()
    maf_df = read_dataframe(maf_file)

    make_dir(out_dir)

    # Write the mutated sequence to this file
    mutated_sequences_fa = out_dir + '/' + sample + '.mutated_sequences.fa'
    out_fa = open(mutated_sequences_fa, 'w')

    cdna_seq_file = load_cds_fasta(cdna)
    cds_seq_file = load_cds_fasta(cds)

    # for key,value in cds_seq_file.items():
    # print ("Trans:" +key)
    # print ("Value :" + value)
    # Extract the sequences (initialize)
    cds_seq = ''
    cdna_seq = ''

    # Extract non-synonymous and their Transcript ID
    for i, j in maf_df.iterrows():
        transcript_id = j['Transcript_ID']
        variant_class = j['Variant_Classification']
        # has to be a SNP
        variant_type = j['Variant_Type']
        chr = j['Chromosome']
        start = j['Start_Position']
        end = j['End_Position']
        hugo_symbol = j['Hugo_Symbol']

        # Find the sequences for non-synonymous variants
        if transcript_id in cds_seq_file:
            cds_seq = cds_seq_file[transcript_id]

        if transcript_id in cdna_seq_file:
            cdna_seq = cdna_seq_file[transcript_id]
            # new_protein = cds_to_aa(cdna_seq)
            # print (new_protein)

        if cds_seq == '':
            print("This mutation of variant class " + variant_class +
                  " and " + transcript_id + " is missing")

        if (type_non_syn(variant_class) and variant_approve(variant_type)):
            # Convert this to protein
            if cds_seq != '' and type_non_syn(variant_class):
                # print ("This is complete with sequence")
                hgvsc = j['HGVSc']

            # print ("Variant " + str(hgvsc) + " is being worked on:")
                (position, ref_allele, alt_allele) = match_hgsc(hgvsc)
                # cds = re.search(cds_seq + '.*', cdna_seq).group()
                seq_5p = cds_seq[0:position]
                # print ("I am printing" + str(position))
                # print ("This is after" + str(position))
                seq_3p = cds_seq[position:len(cds_seq)]
                ref = cds_seq[position]
                # print ("Reference is " + ref )
                # print (transcript_id)
                # print ("Positin " + str(position) + "Reference " + ref_allele)
                # print ("Length of 3p" + str(len(seq_3p)))
                # print(seq_5p)
                wt_cds = seq_5p + ref_allele + \
                    seq_3p[len(ref_allele):len(seq_3p)]
                to_pos = seq_5p + ref_allele
                mut_cds = seq_5p + alt_allele + \
                    seq_3p[len(ref_allele):len(seq_3p)]
                # print ("Length of the cds is " + str(len(cds)))
                # print ("Wild type")
                # print (wt_cds)
                # print (mut_cds)
                # print (wt_cds)
                wt_aa_peptide = cds_to_aa(wt_cds)
                # wt  = cds_to_aa(wt_cds)
                mut_aa_peptide = cds_to_aa(mut_cds)

                # This one is copied from Chai's script:
                len_from_start = len_from_end = 0
                for i in range(0, min(len(wt_aa_peptide), len(mut_aa_peptide))):
                    len_from_start = i
                    if wt_aa_peptide[i:i + 1] != mut_aa_peptide[i:i + 1]:
                        break

                # from end
                wt_rev = wt_aa_peptide[::-1]
                mt_rev = mut_aa_peptide[::-1]
                for i in range(0, min(len(wt_aa_peptide), len(mut_aa_peptide))):
                    len_from_end = i
                    if len_from_end + len_from_start >= min(len(wt_aa_peptide), len(mut_aa_peptide)) or \
                            wt_rev[i:i + 1] != mt_rev[i:i + 1]:
                        break

                wt_start = len_from_start
                wt_end = len(wt_aa_peptide) - len_from_end
                mt_start = len_from_start
                mt_end = len(mut_aa_peptide) - len_from_end
                paddinglen = 10
                # change paddinglen = 6 for 9 mers and even number n=5
                wt_altered_aa = wt_aa_peptide[max(
                    0, wt_start - paddinglen+1):min(len(wt_aa_peptide), wt_end + paddinglen-1)]
                # mut_altered_aa = mut_aa_peptide[max(0, mt_start - paddinglen + 1):min(len(mut_aa_peptide), mt_end + paddinglen-1)]
                mut_altered_aa = mut_aa_peptide[max(
                    0, mt_start - paddinglen+1):min(len(mut_aa_peptide), mt_end + paddinglen-1)]
                header = "" + str(transcript_id) + "_" + str(start) + "_" + \
                    str(end) + "_" + variant_type + "_" + str(hugo_symbol)
                # Down
                down = mt_end + paddinglen-1 + 30
                up = mt_start - paddinglen + 1 - 30
                up_mut_aa = mut_aa_peptide[up:max(
                    0, mt_start - paddinglen + 1)]
                down_mut_aa = mut_aa_peptide[mt_end + paddinglen-1:down]
                start_mut = max(0, mt_start - paddinglen + 1)
                # min_mut = min(len(mut_aa_peptide))
                if len(mut_aa_peptide) > 10:
                    out_fa.write('>' + header + '\n')
                    out_fa.write(mut_altered_aa + '\n')
                    # out_fa.write(up_mut_aa + '\t')
                    # out_fa.write(down_mut_aa + '\n')
                    # out_fa.write(str(max(0, mt_start - paddinglen + 1)) + '\n')
                    # out_fa.write(str(mt_end + paddinglen-1) + '\n')
                    # out_fa.write(up_mut_aa + '\n')
                    # out_fa.write(str(down) + '\n')

                # out_fa.close()


if __name__ == '__main__':
    main()
