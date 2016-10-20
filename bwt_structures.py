'''
@author: Divya Koyyalagunta 
@netid: dk160
@date: October 14, 2016

@note: CS260 PS3

'''
from collections import Counter
from compsci260lib import *

# Note that the '$' character will be used to designate the end of a given string.

def forward_bwt(seq):
    """
    forward_bwt(seq) takes as input a string containing the EOF
    character to which the BWT must be applied. The method should
    then return the result of the BWT on the input string.
    
    For example:
    print forward_bwt('GATTACA$') --> 'ACTGA$TA'
    """
    n = len(seq)
    suffixes = []
    x = n - 1

    while(x >= 0):
        suffixes.append(seq[x:] + seq[0:x])
        x -= 1
    suffixes_sorted = sorted(suffixes)
    
    bwt = ""
    for str in suffixes_sorted:
        bwt += str[n-1]
        
 
    return bwt
    
def reverse_bwt(seq):
    """
    reverse_bwt(seq) takes as input a string containing the EOF
    character to which the reverse of the BWT must be applied.
    The method should then return the result of the reversal on
    the input string.
    
    For example:
    print reverse_bwt('ACTGA$TA') --> 'GATTACA$'
    """    
    lst = []
    for char in seq:
        lst.append(char)
    lst = sorted(lst)   
    i = len(seq) - 1
    while (i > 0):
        for x in range(0,len(seq)):
            lst[x] = seq[x] + lst[x]
        lst = sorted(lst)
        print lst
        i -= 1
    for elem in lst:
        if elem[-1] == '$':
            return elem
    
# It may be helpful to read the documentation for the methods
# given below, but you will NOT have to make any changes to
# them in order to complete the problem set.

def make_suffix_array(seq):
    """
    Returns the suffix array of the input string.
    
    For example:
    print make_suffix_array('GATTACA$') --> [7, 6, 4, 1, 5, 0, 3, 2]
    """
    suffixes = {}
    for x in range(len(seq)):
        suffixes[seq[x:]] = x
    suffix_array = [suffixes[suffix] for suffix in sorted(suffixes.keys())]
    return suffix_array

def rank(bwt_seq):
    """
    Takes as input a string transformed by the BWT. Returns a
    dictionary with characters as keys and lists as values.
    Each list contains the total number of occurrences for the 
    corresponding character up until each position in the
    BWT-transformed string (i.e., its rank).
    
    For example:
    print rank('ACTGA$TA')['$'] --> [0, 0, 0, 0, 0, 1, 1, 1]
    print rank('ACTGA$TA')['A'] --> [1, 1, 1, 1, 2, 2, 2, 3]
    print rank('ACTGA$TA')['C'] --> [0, 1, 1, 1, 1, 1, 1, 1]
    print rank('ACTGA$TA')['G'] --> [0, 0, 0, 1, 1, 1, 1, 1]
    print rank('ACTGA$TA')['T'] --> [0, 0, 1, 1, 1, 1, 2, 2]
    """
    rank = {}
    characters = set(bwt_seq)
    for character in characters:
        rank[character] = [0]
    rank[bwt_seq[0]] = [1]
    for letter in bwt_seq[1:]:
        for k, v in rank.items():
            v.append(v[-1] + (k == letter))
    return rank

def count_smaller_chars(seq):
    """
    Takes as input a string. Returns a dictionary with characters
    as keys and integers as values. The integers track the
    number of characters in the input string which are
    lexicographically smaller than the corresponding character
    key.
    
    For example, using an input DNA sequence like 'GATTACA':
    print count_smaller_chars('GATTACA')['A'] --> 0  (A, being lexicographically first, should always return 0)
    print count_smaller_chars('GATTACA')['C'] --> 3  (C, being second, should return the number of A's, which is 3)
    print count_smaller_chars('GATTACA')['G'] --> 4  (G, being third, should return the number of A's or C's, which is 4)
    print count_smaller_chars('GATTACA')['T'] --> 5  (T, being fourth, should return the number of A's or C's or G's, which is 5)
    """
    characters = set(seq)
    cntr = Counter(seq)
    total = 0
    counts = {}
    for character in sorted(characters):
        counts[character] = total
        total += cntr[character]
    return counts

def make_all(reference):
    """
    Takes as input a reference string. Returns the data
    structures necessary to perform efficient exact string
    matching searches.
    """
    counts = count_smaller_chars(reference)
    reference = reference + '$'
    suffix_array = make_suffix_array(reference)
    bwt = forward_bwt(reference)
    ranks = rank(bwt)
    return bwt, suffix_array, ranks, counts

if __name__ == '__main__':
    print forward_bwt("GAGCTGAT$")
    print reverse_bwt('C$TTTACAG')