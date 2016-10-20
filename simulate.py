'''
@author: Divya Koyyalagunta 
@netid: dk160
@date: October 14, 2016

@note: CS260 PS3

'''
import sys, random
from compsci260lib import *
from math import exp 
import math

def simulate():
      
    G = 3000000
    R = 45000
    L = 500
    iterations = 20

    C_list = []
    contig_list = []
    contig_length_averages_list = []
    nuc_list = []
    for x in range(0, iterations):
        run = x + 1
        print "----------------------------------------------------------------------"
        print "Run number", run
        list = [0]*G
    
        randStarts = []
    
        for x in range (0, R):
            randStarts.append(random.randint(0,G - L))
        for randStart in randStarts:
            for x in range(randStart, randStart + L):
                list[x] += 1
        sum = 0.0  
        nuc_not_covered = 0
        C = 0
        contig_count = 0
        started_contig = False
        zero_positions = []
        for x in range(0, len(list)):
            sum += list[x]
            if list[x] == 0:
                started_contig = False
                nuc_not_covered += 1
                zero_positions.append(x)
            if list[x] != 0 and started_contig == False:
                contig_count += 1
                started_contig = True

        contig_lengths = []
        
        # handle cases when first index is not 0
         
        if list[0] != 0:
            contig_lengths.append((zero_positions[0]))
            
            
        x = 0
        
        while(x < len(zero_positions) - 1):
            if zero_positions[x] + 1 != zero_positions[x + 1]:
                curr_contig_length =  zero_positions[x + 1] - zero_positions[x] - 1 
                contig_lengths.append(curr_contig_length)
            x += 1
                
        # handle cases when last index is not 0
    
        if list[len(list) - 1] != 0:
            last_zero_pos = len(zero_positions) - 1
            contig_lengths.append(len(list) - zero_positions[last_zero_pos] - 1)
                                  
        # print info     
        C = sum/len(list)      
        print "Number of contigs:",contig_count
        print "Average length of a contig:", math.fsum(contig_lengths)/contig_count
        print "C:", C
        print "Number of nucleotides not covered by any read:", nuc_not_covered
        contig_length_averages_list.append(math.fsum(contig_lengths)/contig_count)
        C_list.append(C)
        contig_list.append(contig_count)
        nuc_list.append(nuc_not_covered)
    
    # print average info
    print "----------------------------------------------------------------------"
    print "AVERAGES:"
    print "Average number of contigs:", math.fsum(contig_list)/len(contig_list)
    print "Average length of contigs:", math.fsum(contig_length_averages_list)/len(contig_length_averages_list)
    print "Average value of C:", math.fsum(C_list)/len(C_list)
    print "Average number of nucleotides not covered by any read:", math.fsum(nuc_list)/len(nuc_list)
    
if __name__ == '__main__':
    simulate()
  