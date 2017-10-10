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
    'empty_email_input': 'Sähköpostiosoite tarvitaan internet-hakuun!'
}
# global variable that holds command options
OPTIONS_DICT = {}


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


def read_fasta_file_from_the_internet(organism_name, email):
    '''
    Try to find genome for given organism from the internet in a fasta format.
    (http://biopython.org/DIST/docs/api/Bio.Entrez-module.html)
    
    '''
    
    Bio.Entrez.email = email
    handle = Bio.Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", 
                               id=organism_name)
    seq_record = Bio.SeqIO.read(handle, "gb")
    handle.close()
    #print("Genbank ID:", seq_record.id)
    #print("Annotations:", seq_record.annotations)
    #print("Features:", seq_record.features)
    #print("Sekvenssi:", seq_record.format("fasta"))

    return seq_record.seq
    

def read_test_fasta_file_from_filesystem():
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


def find_islands(nucleotide_seq):
    '''
    Find CpG islands from the given nucleotide sequence
    
    '''  
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


def two_genomes(islands_list_1, islands_list_2):
    '''
    Return true if both genomes contain islands
    
    '''
    return (len(islands_list_1) >= 2 and len(islands_list_2 >= 2))


def do_visuals(island_list):
    '''
    Show some visualization of the CpG data
    
    '''
    if(len(island_list) < 3):
        return
    
    print('visualization...')
    print('Organism: ', island_list[-1])
    

def print_island_details(island_list):
    '''
    Print some details of the data
    
    '''
    if(len(island_list) < 2):
        return
    
    print('Organism: ', island_list[-1])
    print('Number of possible CpG islands: ', len(island_list[1:-1]))
    
    
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
    
    organism1_id = input(TEXTS_DICT['organism_input_1'])
    organism2_id = input(TEXTS_DICT['organism_input_2'])
    email = input(TEXTS_DICT['email_input'])
    islands_list_1 = []
    islands_list_2 = []
    islands_test_list = []
    
    if(organism1_id == 'test'):
        islands_test_list = find_islands(read_test_fasta_file_from_filesystem())
        islands_test_list.append('SH1')
    elif(organism1_id is not '' and email is not ''):
        islands_list_1 = find_islands(read_fasta_file_from_the_internet(
                organism1_id, email))
        islands_list_1.append(organism1_id)
        if(organism2_id is not ''):
            islands_list_2 = find_islands(read_fasta_file_from_the_internet(
                    organism2_id, email))
            islands_list_2.append(organism2_id)
    else:
        if(organism1_id is ''):
            print(TEXTS_DICT['empty_organism_input'])
        if(email is ''):
            print(TEXTS_DICT['empty_email_input'])
             
    #print('islands for test found:', islands_test_list)
    #print('islands for 1 found:', islands_list_1)
    #print('islands for 2 found:', islands_list_2)
    
    print_island_details(islands_test_list)

    # TODO: compare CpG-islands if two genomes were given
    if(two_genomes(islands_list_1, islands_list_2)):
        print_comparison_of_cpg_islands(islands_list_1, islands_list_2)
    else:
        print_island_details(islands_list_1)
    
        
    # TODO: visualize CpG-islands
    if(OPTIONS_DICT['result_graphics'] == 1):
        do_visuals(islands_test_list)
        do_visuals(islands_list_1)
        do_visuals(islands_list_2)


# guard to only execute code when a file is invoked as a script
if __name__ == '__main__':
    start()
