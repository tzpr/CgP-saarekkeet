#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Bio.Entrez
import Bio.SeqIO
from optparse import OptionParser  # obs! deprecated module (but good for now)



# global variable that holds the messages shown to user
TEXTS_DICT = {
    'organism_input_1': 'Anna eliö 1: ',
    'organism_input_2': 'Anna eliö 2: ',
    'email_input': 'Anna email: ',
    'empty_organism_input': 'Anna jokin eliö!',
    'empty_email_input': 'Sähköpostiosoite tarvitaan!'
}
# global variable that holds command options
OPTIONS_DICT = {}


class ResearchSubject:
    name = ''
    genome_sequence = ''
    cpg_islands = [] # list holding indexes of found islands
    cpg_ratio_threshold = 0
    gc_percentage_threshold = 0
    researcher_email = ''
    
    
    def populateSubjectFromInput(self):
        self.name = input('Anna tutkittavan eliön nimi (tai Q niinkuin quit): ')
        
        if self.name == 'Q':
            return
        if self.name == '' or self.name == None:
            raise Exception('Eliön nimi tarvitaan!')
        
        self.researcher_email = input('Anna sähköposti: ')
        
        if self.researcher_email == None or self.researcher_email == '':
            raise Exception('Email tarvitaan!')
        
        # get the genomes
        if self.name == 'test':
            self.genome_sequence = self.__read_test_fasta_file_from_filesystem()
            self.name = 'SH1' # this is the name of the test subject
        else:
            self.genome_sequence = self.__read_fasta_file_from_the_internet(
                    self.name, self.researcher_email)
            
                
    def genome(self):
        return self.genome_sequence 
    
    
    def set_cpg_ratio_threshold(self, cpg_ratio):
        self.cpg_ratio_threshold = ratio
        
        
    def set_gc_percentage_threshold(self, gc_percentage):  
        self.gc_percentage_threshold = gc_percentage
        
        
    def cpg_ratio_threshold(self):
        return self.cpg_ratio_threshold
        
        
    def gc_percentage_threshold(self):  
        return self.gc_percentage_threshold   
      
            
    def __read_test_fasta_file_from_filesystem(self):
        '''
        Mainly for testing purposes.
        '''
        dummy_test_seq = 'CTGGACACCAGCGTAGACCTGCGGTTCAAGTGACCATGCCGGGAATCGTCTCACAGTACGTGCTCCCCGT'
        
        test_fasta_file = 'SH1_genome.fasta'
        genome_seq = []
        
        with open(test_fasta_file) as f: 
            header_line = f.readline()
            for line in f: 
                genome_seq.append(line.strip())
        
        return ''.join(genome_seq) 
    

    def __read_fasta_file_from_the_internet(self, organism_name, email):
        '''
        Try to find genome for given organism from the internet in a fasta format.
        (http://biopython.org/DIST/docs/api/Bio.Entrez-module.html)
        
        '''   
        Bio.Entrez.email = email
        handle = Bio.Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", 
                                   id=organism_name)
        seq_record = Bio.SeqIO.read(handle, "gb")
        handle.close()
    
        return seq_record.seq       
            

def read_command_line_arguments():
    ''' 
    Read and store predefined optional commandline arguments. Uses
    optparse module. 
    Available options: 
        -s (sample window size, default value 200 nukleotides)
        -p (GC persentage, default value 0,5)
        -r (CpG ratio, default value 0,6)
        -h (Help)
        -g (Show graphics of the results, default value is 0, no graphics)
        
    '''
    options_dict = {}
    # initialize options_dict with default values
    options_dict['sample_window_size'] = 200
    options_dict['gc_percentage_threshold'] = 50
    options_dict['cpg_ratio_threshold'] = 60
    options_dict['result_graphics'] = 0

    # https://docs.python.org/2/library/optparse.html
    parser = OptionParser()
    parser.add_option('-s', '--size', dest='sample_window_size',
                      help='set the sample window size.', type=int)
    parser.add_option('-p', '--gcp', dest='gc_percentage',
                      help='set the GC persentage treshold.', type=int)
    parser.add_option('-r', '--ratio', dest='cpg_ratio',
                      help='set the CpG ratio. ' +
                      'Note: computer players continue the game till the end! ',
                      type=int)
    parser.add_option('-g', '--graph', dest='show_graphics',
                      help='activate graphics.', type=int)
    # get the options, discard the leftover arguments with _
    (options, _) = parser.parse_args()

    if options.sample_window_size is not None:
        options_dict['sample_window_size'] = options.sample_window_size

    if options.gc_percentage is not None:
        options_dict['gc_percentage_threshold'] = options.gc_percentage

    if options.cpg_ratio is not None:
        options_dict['cpg_ratio_threshold'] = options.cpg_ratio
        
    if options.show_graphics is not None:
        options_dict['result_graphics'] = options.show_graphics    

    return options_dict


def island_rule_1_ok(nucleotide_seq_str):
    '''
    Check if the CpG sites rule #1 is satisfied
    A GC percentage greater than 50% (https://en.wikipedia.org/wiki/CpG_site)
    '''      
    if(len(nucleotide_seq_str) == 0):
        return False
    
    default_value = 50
    gc_percentage_threshold = OPTIONS_DICT['gc_percentage_threshold'] or default_value
    
    sample_seq_length = len(nucleotide_seq_str)
    g_nucleotides_count = nucleotide_seq_str.count('G')
    c_nucleotides_count = nucleotide_seq_str.count('C')
    
    g_percentage = (g_nucleotides_count/sample_seq_length) * 100
    c_percentage = (c_nucleotides_count/sample_seq_length) * 100
    gc_percentage = g_percentage + c_percentage
    
    if(gc_percentage >= gc_percentage_threshold):
        return True
    
    return False


def island_rule_2_ok(nucleotide_seq_str):
    '''
    Check if the CpG sites rule #2 is satisfied
    An observed-to-expected CpG ratio greater than 60 % 
    (https://en.wikipedia.org/wiki/CpG_site)
    '''    
    default_value = 60
    
    cpg_ratio_threshold = OPTIONS_DICT['cpg_ratio_threshold'] or default_value
    
    gc_pairs_obs = nucleotide_seq_str.count('CG')
    gc_pairs_exp = (nucleotide_seq_str.count('C') * 
                    nucleotide_seq_str.count('G')) / len(nucleotide_seq_str)
    if(gc_pairs_exp == 0):
        return False
    
    cpg_ratio = (gc_pairs_obs/gc_pairs_exp) * 100
    
    if(cpg_ratio > cpg_ratio_threshold):
        return True
    
    return False


def island_conditions_ok(seq):
    '''
    Check if island conditions are met
    '''
    return (island_rule_1_ok(seq) and island_rule_2_ok(seq))


def print_comparison_of_cpg_islands(islands_list_1, islands_list_2):
    '''
    Some comparison between two genomes
     
    Lists is constructed in a following way:
    
    list[0] first elements contains the full sequence of the genome
    list[1]...[n] contain the possible CpG islands (tuples containing islands
                  start and end indexes )
    list[n+1] the last element contains the name of the organism
    
    '''
    
    print('*****************************************************************')
    print('')
    print('Comparison of CpG islands in ', islands_list_1[-1], 'and', 
          islands_list_2[-1])
    print('')
    print('Organism: ', islands_list_1[-1])  
    print('Number of possible CpG islands: ', len(islands_list_1[1:-1]))
    print('')
    print('Organism: ', islands_list_2[-1])  
    print('Number of possible CpG islands: ', len(islands_list_2[1:-1]))
    print('')
    print('')
    print('*****************************************************************')


def find_islands(subject):
    '''
    Find CpG islands from the given nucleotide sequence
    
    '''  
    nucleotide_seq = subject.genome()
    default_window = 200
    SAMPLE_LENGTH = OPTIONS_DICT['sample_window_size'] or default_window
    start_index = 0
    end_index = SAMPLE_LENGTH
    sample_seq = nucleotide_seq[start_index:end_index]
    cpg_islands_list = []
    genome_length = len(nucleotide_seq)
    # let's put the original genome sequence first
    cpg_islands_list.append(nucleotide_seq)
    
    print('genome_length: ', genome_length)
     
    while(True):
        if(island_conditions_ok(sample_seq)):
            end_index = end_index + 1
            if(end_index > genome_length):
                cpg_islands_list.append((start_index, end_index-1))
                print('break 1:', end_index)
                break
            
            sample_seq = nucleotide_seq[start_index:end_index]
        else:
            if(len(sample_seq) == SAMPLE_LENGTH):
                # no island in default sample, lets move on
                start_index = start_index + 1
                end_index = start_index + SAMPLE_LENGTH
                if(end_index > genome_length):
                    print('break 2:', end_index)
                    break
                
                sample_seq = nucleotide_seq[start_index:end_index]
            else:
                cpg_islands_list.append((start_index, end_index))
                start_index = start_index + 1 # move start index by one
                end_index = start_index + SAMPLE_LENGTH
                if(end_index > genome_length):
                    print('break 3:', end_index)
                    break
                sample_seq = nucleotide_seq[start_index:end_index]
                
    return cpg_islands_list

    
def create_study_subjects():
    subjects = []
    
    while 1:
        subject = ResearchSubject()
        subject.populateSubjectFromInput()
	   
        if subject.name == 'Q':
            break

        subjects.append(subject)
       
    return subjects
     
        
def search_islands(subjects):
    '''
    Gets list of organisms whos genomes are studied for islands
    
    '''
    for subject in subjects:
        find_islands(subject)
    
    
def visualize_results(subject_list):
    print('visualize_results NOT IMPLEMENTED YET')
    
    
def compare_results(subject_list):
    print('compare_results NOT IMPLEMENTED YET')
    
    
def start():
    '''
    The main thing. Asks user input.
    
    If user gives "test" as the name of the organism then test fasta file is 
    used in island search.
    
    If user gives some other name then organism's genome sequence is fetched
    from the internet and used in the island search.
    
    The organism 1 is mandatory and also email if not "test" given as name.
    The organism 2 is optional. 
    
    If both organism are given then the results are compared in the end. 
    
    Some visualizations are shown if user has activated the graphics using
    options.
    
    '''
    global OPTIONS_DICT 
    OPTIONS_DICT= read_command_line_arguments()
    
    subjects = []
    
    subjects = create_study_subjects()
    
    search_islands(subjects)
    
    compare_results(subjects)

    visualize_results(subjects)


# guard to only execute code when a file is invoked as a script
if __name__ == '__main__':
    start()
