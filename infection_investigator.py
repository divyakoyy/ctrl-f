'''
@author: Divya Koyyalagunta 
@netid: dk160
@date: October 14, 2016

@note: CS260 PS3

'''

from bwt_structures import *
from read_aligner import *
from compsci260lib import *
from operator import pos

def reverse_complement(seq):
    """
    Returns the reverse complement of the input string.
    """
    comp_bases = {'A': 'T',
                  'C': 'G',
                  'G': 'C',
                  'T': 'A'} 
    rev_seq = list(seq)
    rev_seq = rev_seq[::-1]
    rev_seq = [comp_bases[base] for base in rev_seq] 
    return ''.join(rev_seq)
    
def align_patient_reads():
    """
    Given the genomes of different bacteria and 100,000 reads from each patient,
    uses the find method from read_aligner to map each read to a genome. Reports
    the prevalence of each microbe in each patient's genome 
    """
    
    patients = ['patient1','patient2','patient3']
    bacteria = [ 'Bacteroides_ovatus', 'Bacteroides_thetaiotaomicron','Bifidobacterium_longum','Eubacterium_rectale','Lactobacillus_acidophilus','Peptoniphilus_timonensis','Prevotella_copri','Roseburia_intestinalis','Ruminococcus_bromii','Vibrio_cholerae']

    
    bacteria_info = {}
    for bacterium in bacteria:
        reference = get_fasta_dict('reference_genomes/' + bacterium + '.fasta')
        genome =  reverse_complement(reference[' '.join(bacterium.split('_'))])
        bwt_data = make_all(genome)
        bacteria_info[bacterium] = bwt_data
   
    count_vectors_list = []

    for patient in patients:
        reads_dict = get_fasta_dict('patients/' + patient + '.fasta')
        reads_mapped = {}
        count_vector = 15000*[0]
        num_unique_reads = 0
        ct = 0
        genomes_mapped = {}
        for read in reads_dict:
            already_found = False
            for bacterium in bacteria:
                results = find(reads_dict[read], bacteria_info[bacterium])
                # add result to reads_mapped dictionary if there is a unique match found
                if len(results) != 0:                     
                    num_unique_reads += 1
                    if bacterium in genomes_mapped and already_found == False:
                        genomes_mapped[bacterium] += 1
                    elif already_found == False: 
                        genomes_mapped[bacterium] = 1
                        
                    if(bacterium == 'Vibrio_cholerae'):
                        for pos in results:
                            x = pos
                            while(x < pos + len(reads_dict[read])):
                                count_vector[x] += 1
                                x += 1
                    already_found = True
                    
        # print info
        print "------------------------------------------------------------------------------------"
        
        print patient
        print "Number of unique reads: " + str(num_unique_reads)
        for genome_name in genomes_mapped:
            print genome_name + ":" + str(float(genomes_mapped[genome_name])/num_unique_reads)
    
        range = longest_internal_zero(count_vector)
        print "Start and stop position of longest stretch of internal zeroes: " + str(range)


def longest_internal_zero(count_vector):

    nonzero_positions = []
    for x in range(0,len(count_vector)):
        if count_vector[x] != 0:
            nonzero_positions.append(x)

    x = 0
    internal_strings_starts_stops = []
    while(x < len(nonzero_positions) - 1):
        if nonzero_positions[x] + 1 != nonzero_positions[x + 1]:
            curr_string_length =  nonzero_positions[x + 1] - nonzero_positions[x] - 1 
            internal_strings_starts_stops.append((curr_string_length,nonzero_positions[x] + 1,nonzero_positions[x + 1] - 1))
        x += 1
    sorted_list = sorted(internal_strings_starts_stops)
      
    if len(sorted_list) == 0:
        return "no internal zeroes"
    return sorted_list[-1][1],sorted_list[-1][2]

if __name__ == '__main__':
    align_patient_reads()
   