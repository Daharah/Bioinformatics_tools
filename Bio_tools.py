# Conor's Bioinformatics DNA anaylsis toolkit

#Definitions: C_ = conversion tool

import random 
import math as m
import Bio_structures as bio_s

Nucleotides = ['A', 'C', 'G', 'T']

test_seq = ('TTGAAGATCCTTTCACCTTCCTTCCGTCATATTCCAGCCCTACACGCAGGTACTGTTTCAACTCGAATCTAGTTTGAAACTTCATCGGTCTAGCCTCTGATGCCTAGCCCGCCAATAAAATACTAGATATGACCTATCTTTTAAGATTCGCGTCACGCTGTGTGAGTGAAAACACCCCCAAATATTTAATCCATGAAAGGTGGACTCTGCCCCGGCTAAGTTGTATTCTCCAGCGCGACGGTGCTCTTACGCGGCGGAGTCTCTCGTCGGTAAAGCTCCCAGGGGACCTTCGGCCGCACAAGTATACAGGCTCCAACCAGTGTGGATATGAAGTTTGTCGCCCTTGGCAACCTCGAACCGGGTTTGCCGCGAATGTACGGCTGACACCGTTCTCGCCGTGCCCGCCACGCGACATCGTCAAGTGTTGTATATTAGCCGACCCCCTCATTTATAATGAGGGAACACTTTTACAGTGGCTAAGTATGAGTGTCAAGGTGGTGCCTGCCGAGTTTCTGCGACCCCAACCGCTCCCTCTAAATTCATACCACCAAAGGTGACATGTACAGCAGCTCATTGCTACGCAATAGTGGCGTCTATTGCATTTAGTGATACGAGGAACCTCGAACCATCGGGGCGTCGGCCCAATCAGCCGCACGTTCCCCGGACCGTCACCGGACAGAGTCTCGAGACCGCTGGAAGAATCGGGTTGGGGGCAATGTCGTCGCTTCCGCGTCGGTCATACGATACTCTCGCCGTGCCGAGATACTGGGTGCATACCATCTTGAGGCTCTCGGCCTCAGCAGTTGAGTAGGATAGTTTCACCGGATTTGGCCTCCGAAATCCTGGTTACCTAGAGCTGGCCAACTTTGGAGGGAGTCCCTACTATATCCTTAGGGGAGGTCGGTCGATGGCGGCGGACGAGTGGTTACGTATCCGGGTTTGCCGATGCTGTAGCAGGCATTAGGGTTAACGTTATACGCGTGCATCAATTCGCGTCAT')

############ Conversion tools ###########

def tscribe_cdna(dna_seq0):
    """Coding DNA to RNA. Thymine to Uracil"""
    return dna_seq0.replace("T", "U")

def rev_tscribe_2cdna(rna_seq0):
    """Reverse transcription of RNA to DNA. Uracil to Thymine."""
    return rna_seq0.replace("U", "T")

def comp_dna_seq(dna_seq):
    """Generate a complimentary DNA sequence."""
    comp_nuc = bio_s.DNA_comp_nuc
    return ''.join(comp_nuc[nuc] for nuc in dna_seq)

def rev_comp_dna_seq(dna_seq1):
    """Generate a reverse complimentary DNA sequence."""
    comp_nuc = bio_s.DNA_comp_nuc
    return ''.join(comp_nuc[nuc] for nuc in dna_seq1[::-1])

def tscribe_tdna(dna_seq2):
    """Transcribe an RNA sequence from template DNA."""
    comp_nuc = {'A':'U', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(comp_nuc[nuc] for nuc in dna_seq2[::-1])

def tlate_cdna(dna_seq3):
    """Translate a coding DNA sequence into an amino sequence."""
    nuc2amino = bio_s.DNA_Codons
    amino_seq = ""
    for i in range (0, len(dna_seq3), 3): 
        codon = dna_seq3[i:i+3]
        if len(codon) == 3:
            amino_seq += nuc2amino.get(codon, 'Invalid sequence')
    return(amino_seq)

def tlate_tdna(dna_seq4):
    """Translate a template DNA sequence into an amino sequence."""
    nuc2amino = bio_s.DNA_Codons
    amino_seq = ""
    for i in range (0, len(dna_seq4), 3): 
        codon = dna_seq4[i:i+3]
        if len(codon) == 3:
            amino_seq += nuc2amino.get(codon, 'Invalid sequence')
    return(amino_seq[::-1])

def tlate_mrna(rna_seq):
    """Translate an mRNA sequence into an amino sequence."""
    nuc2amino = bio_s.RNA_Codons
    amino_seq = ""
    for i in range (0, len(rna_seq), 3): 
        codon = rna_seq[i:i+3]
        if len(codon) == 3:
            amino_seq += nuc2amino.get(codon, 'Invalid sequence')
    return(amino_seq)

# def conv_short2long_amino(dna_seq5):

# def conv_long2short_amino(dna_seq6):

########## Calculations (calc_) ##########

def calc_count_nucleotides(dna_rna_seq):
    """Counts the individual number of each nucleotide in a DNA/RNA sequence, and sequence length."""
    dna_rna_seq = dna_rna_seq.upper()
    nuc_freq_dict = {'A': 0, 'C': 0, "G": 0, 'T': 0, 'U': 0, "Total:": 0}
    rna_seq = 'U' in dna_rna_seq
    for nuc in dna_rna_seq:
        if nuc in nuc_freq_dict:
            nuc_freq_dict[nuc] += 1
            nuc_freq_dict['Total:'] += 1
    if rna_seq:
        del nuc_freq_dict['T']
    else: 
        del nuc_freq_dict['U']
    return nuc_freq_dict

def calc_melting_temp(dna_seq1):
    AT_count = dna_seq1.count('A') + dna_seq1.count('T')
    CG_count = dna_seq1.count('C') + dna_seq1.count('G')
#    total_count = int(AT_count + CG_count)
#    CG_percent = (round(CG_count / total_count * 100))
    return (f' 2°C * {AT_count} + 4°C * {CG_count} = {AT_count * 2 + CG_count * 4}°C ')
#    if total_count <= 13:     
#        return(f'GC% = {CG_percent}% Tm = {AT_count * 2 + CG_count * 4}°C ')
#    elif total_count <=50:   
#         return(f'GC% = {CG_percent}% Tm = {64.9 = 41 * CG_count - 16.4 / total_count}')  #tm = 64.9 + 41 * ((G + C - 16.4) / N)
#    else:                   
#        return(f'GC% = {CG_percent}% Tm = {81.5 + 16.6 * m.log10(Na) + 0.41 * CG_percent - 500 / total_count}') #tm = 81.5 + 16.6 * math.log10(Na(assumed to be 0.05)) + 0.41 * GC_percent - (500 / N)


########## Utility ##########

# Random DNA string generator
u_randomdna = "".join([random.choice(Nucleotides) 
                     for nuc in range(100)])

# DNA sequence validator (Reject if not A, C, G, or T)
def u_validate_dna_seq(valseq):
    cap_seq = valseq.upper()
    for nuc in cap_seq:
        if nuc not in Nucleotides:
            return False
    return cap_seq

# def mass(pro_seq):
#     mass_table = bio_s.protein_mass
#     for letter, weight in mass_table:
#             letter