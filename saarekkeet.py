#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Bio.Entrez
import Bio.SeqIO



# global variable that holds the messages shown to user
TEXTS_DICT = {
    'invalid_condition': 'ERROR ERROR!',
    'organism_input_1': 'Anna eliö 1: ',
    'organism_input_2': 'Anna eliö 2: ',
    'email_input': 'Anna email: ',
    'empty_organism_input': 'Anna jokin eliö!',
    'empty_email_input': 'Sähköpostiosoite tarvitaan internet-hakuun!'
}

def read_fasta_file_from_the_internet(organism_name, email):
    '''
    Try to find genome for given organism from the internet in a fasta format.
    (http://biopython.org/DIST/docs/api/Bio.Entrez-module.html)
    
    Throws exception in case of an unknown organism.
    '''
    # TODO verify email and organism_name
    
    # working implementation
    # Bio.Entrez.email = email
    #handle = Bio.Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=organism_name)
    #seq_record = Bio.SeqIO.read(handle, "gb")
    #handle.close()
    #print("Genbank ID:", seq_record.id)
    #print("Annotations:", seq_record.annotations)
    #print("Features:", seq_record.features)
    #print("Sekvenssi:", seq_record.seq)    
    
    
    print('read_fast_from_the_internet')
    

def read_fasta_file_from_filesystem(file_name):
    '''
    Mainly for testing purposes.
    '''
    # TODO: read fasta file
    
    dummy_test_seq = 'CTGGACACCAGCGTAGACCTGCGGTTCAAGTGACCATGCCGGGAATCGTCTCACAGTACGTGCTCCCCGT'
    
    return dummy_test_seq
    

def parse_genome_from_fasta_file(file):
    '''
    Read genome from fasta file to an immutable list (tuple) of nucleotides.
    '''
    nucleotides_str = '' 
    
    # TODO: parse file and populate the nucleotides_str
    
    dummy_test_seq = 'CTGGACACCAGCGTAGACCTGCGGTTCAAGTGACCATGCCGGGAATCGTCTCACAGTACGTGCTCCCCGT'
    
    # return nucleotides_str
    return dummy_test_seq


def island_rule_1_ok(nucleotide_seq_str):
    '''
    Check if the CpG sites rule #1 is satisfied
    A GC percentage greater than 50% (https://en.wikipedia.org/wiki/CpG_site)
    '''    
    # calculations could be in a separate functions
    
    sample_seq_length = len(nucleotide_seq_str)
    g_nucleotides_count = nucleotide_seq_str.count('G')
    c_nucleotides_count = nucleotide_seq_str.count('C')
    
    g_percentage = (g_nucleotides_count/sample_seq_length) * 100
    c_percentage = (c_nucleotides_count/sample_seq_length) * 100
    gc_percentage = g_percentage + c_percentage
    
    if(gc_percentage >= 50):
        return True
    
    return False


def island_rule_2_ok(nucleotide_seq_str):
    '''
    Check if the CpG sites rule #2 is satisfied
    An observed-to-expected CpG ratio greater than 60 % 
    (https://en.wikipedia.org/wiki/CpG_site)
    '''    
    # calculations could be in a separate functions
    
    gc_pairs_obs = nucleotide_seq_str.count('GC')
    gc_pairs_exp = (nucleotide_seq_str.count('C') * 
                    nucleotide_seq_str.count('G')) / len(nucleotide_seq_str)
    if(gc_pairs_exp == 0):
        return False
    
    cpg_ratio = (gc_pairs_obs/gc_pairs_exp) * 100
    
    if(cpg_ratio > 60):
        return True
    
    return False


def island_conditions_ok(seq):
    '''
    Check if island conditions are met
    '''
    return (island_rule_1_ok(seq) and island_rule_2_ok(seq))


def find_islands(nucleotide_seq):
    '''
    Find CpG islands from the given nucleotide sequence
    
    '''    
    cpg_islands_list = []
    SAMPLE_LENGTH = 20
    start_index = 0
    end_index = SAMPLE_LENGTH
    sample_seq = nucleotide_seq[start_index:end_index]
     
    while(len(sample_seq) >= SAMPLE_LENGTH):
        if(island_conditions_ok(sample_seq)):
            end_index = end_index + 1
            sample_seq = nucleotide_seq[start_index:end_index]
        elif(not island_conditions_ok(sample_seq)):
            if(len(sample_seq) == SAMPLE_LENGTH):
                start_index = start_index + 1
                end_index = start_index + SAMPLE_LENGTH
                sample_seq = nucleotide_seq[start_index:end_index]
            else:
                cpg_islands_list.append(sample_seq)
                start_index = end_index
                end_index = start_index + SAMPLE_LENGTH
                sample_seq = nucleotide_seq[start_index:end_index]
        else:
            raise Exception(TEXTS_DICT['invalid_condition'])
            
    return cpg_islands_list
    
    
def start():
    '''
    The main thing.
    '''
    test_fasta_file = 'SH1_genome.fasta'
    
    organism1 = input(TEXTS_DICT['organism_input_1'])
    organism2 = input(TEXTS_DICT['organism_input_2'])
    email = input(TEXTS_DICT['email_input'])
    islands_org1 = []
    islands_org2 = []
    islands_test = []
    
    if(organism1 == 'testing'):
        islands_test = find_islands(parse_genome_from_fasta_file(
                read_fasta_file_from_filesystem(test_fasta_file)))
    elif(organism1 is not '' and email is not ''):
        islands_org1 = find_islands(parse_genome_from_fasta_file(
                read_fasta_file_from_the_internet(organism1, email)))
        if(organism2 is not ''):
            islands_org2 = find_islands(parse_genome_from_fasta_file(
                    read_fasta_file_from_the_internet(organism2, email)))
    else:
        if(organism1 is ''):
            print(TEXTS_DICT['empty_organism_input'])
        if(email is ''):
            print(TEXTS_DICT['empty_email_input'])
             
    print('islands found:', islands_org1) 

    # TODO: compare CpG-islands if two genomes were given 
    
    # TODO: visualize CpG-islands in genome


# guard to only execute code when a file is invoked as a script
if __name__ == '__main__':
    start()
