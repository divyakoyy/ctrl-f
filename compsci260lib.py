'''
COMPSCI 260 Python library
version 1.0
@author: YM
@date: 2-5-13

A basic library module that contains some useful features.
'''


import sys, re

codon_to_aa_dict = {}
codon_to_aa_dict[r"GC."] = "A"                          #Alanine
codon_to_aa_dict[r"[TU]G[TUC]"] = "C"                   #Cysteine
codon_to_aa_dict[r"[TUC]"] = "D"                        #Aspartic Acid
codon_to_aa_dict[r"GA[AG]"] = "E"                       #Glutamic Acid
codon_to_aa_dict[r"[TU][TU][TUC]"] = "F"                #Phenylalanine
codon_to_aa_dict[r"GG."] = "G"                          #Glycine
codon_to_aa_dict[r"CA[TUC]"] = "H"                      #Histidine
codon_to_aa_dict[r"A[TU][TUCA]"] = "I"                  #Isoleucine
codon_to_aa_dict[r"AA[AG]"] = "K"                       #Lysine
codon_to_aa_dict[r"[TU][TU][AG]|C[TU]"] = "L"           #Leucine
codon_to_aa_dict[r"A[TU]G"] = "M"                       #Methionine
codon_to_aa_dict[r"AA[TUC]"] = "N"                      #Asparagine
codon_to_aa_dict[r"CC."] = "P"                          #Proline
codon_to_aa_dict[r"CA[AG]"] = "Q"                       #Glutamine
codon_to_aa_dict[r"CG.|AG[AG]"] = "R"                   #Arginine
codon_to_aa_dict[r"[TU]C.|AG[TUC]"] = "S"               #Serine
codon_to_aa_dict[r"AC."] = "T"                          #Threonine
codon_to_aa_dict[r"G[TU]."] = "V"                       #Valine
codon_to_aa_dict[r"[TU]GG"] = "W"                       #Tryptophan
codon_to_aa_dict[r"[TU]A[TUC]"] = "Y"                   #Tyrosine
codon_to_aa_dict[r"[TU]A[AG]|[TU]GA"] = "*"             #Stop  

def codon_2_amino_acid(codon):
    """Uses regular expression encodes a codon to the respective amino acid.
       (adapted from Pallavi's code)"""
    
    for codon_regex in codon_to_aa_dict:
        if re.search(codon_regex, codon.upper()) is not None:
            return codon_to_aa_dict[codon_regex]
    
    print "unresolved codon " + codon
    sys.exit('Error!')  #if failed, then output "Error!" to console


def translate(dna):
    """A subroutine to translate an input DNA sequence into a peptide."""
    
    protein = ""
    
    for i in range(0, len(dna)-2, 3):
        protein += codon_2_amino_acid(dna[i:i+3])
     
    return protein
    
def get_fasta_dict(filename):
    """Given a fasta input file, return a dictionary with the name of each sequence as a key, and the
       actual sequence as the value.  Thus, if you wrote dict = get_fasta_dict(filename), you could 
       see the sequence named read1 (given that it exists in the fasta file) using dict["read1"].
       To see all of the sequence names, you could use dict.keys().  Etc."""
    
    filename.rstrip();
    fasta_dict = {}
    
    try:
        f = open(filename, "r")
    except IOError: 
        print "Error: Unable to open file: " + filename 
    
    try:
        curr_seq = ''
        curr_seq_key = ''
        
        for line in f.readlines():
            #get rid of leading and trailing line
            line = line.strip()
            
            if line is '' or line.isspace() or line[0] == '#':
                continue
            elif line[0] == '>':
                if curr_seq_key is not '' and curr_seq is not '':
                    fasta_dict[curr_seq_key] = curr_seq
                curr_seq = ''
                curr_seq_key = line[1:]
            else:
                curr_seq += line
            
       
            #update the latest line
        fasta_dict[curr_seq_key] = curr_seq
        f.close()
    
    except:
        print "Error reading fasta file: make sure the file format is correct"
        raise
    
    return fasta_dict

def max_over_indices(array, indices):
    """Return the maximum value over the indices (list) given. If indices
       are out of range, a None object will be returned."""
    
    if any(t < 0 for t in indices):
        return None
    if max(indices) >= len(array):
        return None
    
    temp = []
    for i in indices:
        temp.append(array[i])
    
    return max(temp)

def print_abbrv(string):
    """Abbreviate long strings for printing to the console: certain versions/configurations of Eclipse
       have a bug where long strings will be printed "invisibly" for some mysterious reason.  That is, 
       a variable may actually contain a long string, but when you ask Eclipse to print it to the screen,
       it doesn't seem to appear.  Hence, this function can help you ensure your variable really contains
       what you think it contains."""
    
    l = len(string)
    if l < 100:
        print string
    else:
        ellipsis = "...%d more characters..." % (l-70)
        print string[:50] + ellipsis + string[-20:]

