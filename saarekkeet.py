#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Bio.Entrez
import Bio.SeqIO



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
    print('read_fasta_from_filesystem')
    

def parse_genome_from_fasta_file(file):
    '''
    Read genome from fasta file to an immutable list (tuple) of nucleotides.
    '''
    nucleotides_str = '' # 
    
    # TODO: parse file and populate the nucleotides_str
    
    dummy_test_seq = 'CTGGACACCAGCGTAGACCTGCGGTTCAAGTGACCATGCCGGGAATCGTCTCACAGTACGTGCTCCCCGT'
    
    # return nucleotides_str
    return dummy_test_seq


def island_rule_1_ok(nucleotide_seq_str):
    '''
    Check if the CpG sites rule #1 is satisfied
    A GC percentage greater than 50% (https://en.wikipedia.org/wiki/CpG_site)
    '''
    verdict = False
    
    sample_seq_length = len(nucleotide_seq_str)
    gc_pairs_count = nucleotide_seq_str.count('GC')
    gc_nucleotides_count = 2 * gc_pairs_count
    
    gc_percentage = (gc_nucleotides_count/sample_seq_length) * 100
    
    if(gc_percentage >= 50):
        verdict = True
    
    return verdict


def island_rule_2_ok(nucleotide_seq_str):
    '''
    Check if the CpG sites rule #2 is satisfied
    An observed-to-expected CpG ratio greater than 60 % (https://en.wikipedia.org/wiki/CpG_site)
    '''
    print('island_rule_2_ok')
    
    return false


# TODO: what should this return?
def find_islands(nucleotide_seq):
    '''
    Find CpG islands from the given nucleotide sequence
    
    '''
    print('find_islands')
    
    first_200 = nucleotide_seq[:40] 
    
    print(first_200)
    island_rule_1_ok(first_200)
    
    
    island_min_length = 200
    index = 0
    
    for nucleotide in nucleotide_seq:
        index = index + 1
        #print(nucleotide)

    

def start():
    '''
    The main thing.
    '''
    test_fasta_file = 'SH1_genome.fasta'
    
    organism1 = input('Anna eliö 1: ')
    organism2 = input('Anna eliö 2: ')
    email = input('Anna email: ')
    
    if(organism1 == 'testing'):
        fasta_for_org1  = read_fasta_file_from_filesystem(test_fasta_file)
    elif(organism1 is not '' and email is not ''):
        find_islands(parse_genome_from_fasta_file(
                read_fasta_file_from_the_internet(organism1, email)))
        if(organism2 is not ''):
            find_islands(parse_genome_from_fasta_file(
                    read_fasta_file_from_the_internet(organism2, email)))
    else:
        if(organism1 is ''):
            print('Anna jokin eliö!')
        if(email is ''):
            print('Sähköpostiosoite tarvitaan internet-hakuun!')


# guard to only execute code when a file is invoked as a script
if __name__ == '__main__':
    start()
