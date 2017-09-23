#!/usr/bin/env python3
# -*- coding: utf-8 -*-



def read_fasta_file_from_the_internet(organism_name):
    '''
    Try to find genome for given organism from the internet in a fasta format.
    Throws exception in case of an unknown organism.
    '''
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
    nucleotides_tuple = () # lets use immutable list, tuple (we don't need mutations at this point!)
    
    # TODO: parse file and populate the tuple
    
    return nucleotides_tuple


def island_rule_1_ok(nucleotide_list):
    '''
    Check if the CpG sites rule #1 is satisfied
    A region with at least 200 bp (https://en.wikipedia.org/wiki/CpG_site)
    '''
    print('island_rule_1_ok')
    
    return false


def island_rule_2_ok(nucleotide_list):
    '''
    Check if the CpG sites rule #2 is satisfied
    A GC percentage greater than 50% (https://en.wikipedia.org/wiki/CpG_site)
    '''
    print('island_rule_2_ok')
    
    return false


def island_rule_3_ok(nucleotide_list):
    '''
    Check if the CpG sites rule #3 is satisfied
    An observed-to-expected CpG ratio greater than 60 % (https://en.wikipedia.org/wiki/CpG_site)
    '''
    print('island_rule_3_ok')
    
    return false


def start():
    '''
    The main thing.
    '''
    
    print('started!')



# guard to only execute code when a file is invoked as a script
if __name__ == '__main__':
    start()
