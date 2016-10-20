'''
@author: Divya Koyyalagunta 
@netid: dk160
@date: October 14, 2016

@note: CS260 PS3

'''
from bwt_structures import *
from compsci260lib import *

def find(query, bwt_data):
    """
    Given a query sequence and a series of data structures
    containing various information about the reference genome,
    return a list containing all the locations of the query
    sequence in the reference genome. 
    """
    
    bwt, suffix_array, ranks, counts = bwt_data
      
    length = len(bwt)
    results = []
    new_beginning = 0
    new_end = 0
    
    n = len(query)
    x = n - 1
    character = query[x]
    beginning, end = find_range(character, counts, length)
    x -= 1

    while(x >= 0):
        character = query[x]
        if character in ranks:
            if beginning != 0 and ranks[character][beginning] != ranks[character][beginning - 1]:
                beginning = counts[character] + ranks[character][beginning] 
            else:
                beginning = counts[character] + ranks[character][beginning] + 1
                if beginning > len(bwt) - 1:
                    return []
            end = counts[character] + ranks[character][end]
        else:
            return []
        x -= 1

    results = suffix_array[beginning:end + 1]
    return sorted(results)

def find_range(character,counts, length):
    
    if character == 'A':
        return 1, counts['C']
    if character == 'C':
        return counts['C'] + 1, counts['G']
    if character == 'G':
        return counts['G'] + 1, counts['T'] 
    if character == 'T':
        return counts['T'] + 1, length - 1
    else:
        return -1
        
if __name__ == '__main__':
    query = "CAG"
    genome = "TACACAGTTACGACAG"
    bwt_data = make_all(genome)
    print find(query, bwt_data)