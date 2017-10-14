'''
Search for CpG islands
'''

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Bio.Entrez
import Bio.SeqIO


class Genome:
    '''
    Genome under research
    '''
    name = ''
    genome_sequence = ''
    possible_cpg_islands = []
    cpg_ratio_threshold = 0
    gc_percentage_threshold = 0
    window_size = 0
    researcher_email = ''

    def populate_subject_from_input(self):
        '''
        Populate genome from input
        '''
        self.name = input('  Tutkittavan eliön nimi (tai Q niinkuin quit): ')

        if self.name == 'Q':
            return
        if self.name == '' or self.name is None:
            raise Exception('Eliön nimi tarvitaan!')

        self.researcher_email = input('  Anna sähköposti: ')

        if self.researcher_email is None or self.researcher_email == '':
            raise Exception('Email tarvitaan!')

        # get the genomes
        if self.name == 'test':
            test_fasta_file = 'SH1_genome.fasta'
            self.genome_sequence = read_fasta_file_from_filesystem(test_fasta_file)
            self.name = 'SH1' # this is the name of the test subject
        else:
            self.genome_sequence = get_sequence_from_the_internet(
                self.name, self.researcher_email)

    def display(self, verbose):
        ''' Print information of the genome including island search results '''
        print('')
        print('Genome Information:')
        print('- name:', self.name)
        print('- sequence length:', self.get_genome_length())
        print('- number of possible CpG islands:', len(self.possible_cpg_islands))
        print('Search parameters for CpG islands:')
        print('- window_size:', self.window_size)
        print('- cpg_ratio_threshold:', self.cpg_ratio_threshold)
        print('- gc_percentage_threshold:', self.gc_percentage_threshold)
        print('')
        if verbose:
            print('verbose mode, show detailed island info')

    def get_name(self):
        ''' Return genome name '''
        return self.name

    def get_genome(self):
        ''' Return genome sequence string '''
        return self.genome_sequence

    def set_window_size(self, size):
        ''' Set the sample window size '''
        self.window_size = size

    def get_window_size(self):
        ''' Return the sample window size '''
        return self.window_size

    def set_cpg_ratio_threshold(self, cpg_ratio):
        ''' Set the CpG ratio '''
        self.cpg_ratio_threshold = cpg_ratio

    def get_cpg_ratio_threshold(self):
        ''' Return Cpg ratio '''
        return self.cpg_ratio_threshold

    def set_gc_percentage_threshold(self, gc_percentage):
        ''' Set the GC percentage '''
        self.gc_percentage_threshold = gc_percentage

    def get_gc_percentage_threshold(self):
        ''' Return GC percentage value '''
        return self.gc_percentage_threshold

    def get_genome_length(self):
        ''' Return the length of the genome sequence '''
        return len(self.genome_sequence)

    def set_cpg_island(self, island):
        ''' Add CpG island to island list '''
        self.possible_cpg_islands.append(island)

    def set_cpg_islands(self, islands):
        ''' Set CpG islands '''
        self.possible_cpg_islands = islands

    def get_cpg_islands(self):
        ''' Return Cpg islands list '''
        return self.possible_cpg_islands


def read_fasta_file_from_filesystem(file_name):
    '''
    Mainly for testing purposes.
    '''
    dummy_test_seq = 'CTGGACACCAGCGTAGACCTGCGGTTCAAGTGACCATGCCGGGAATCGTCTCACAGTACGTGCTCCCCGT'
    genome_seq = []

    with open(file_name) as fasta:
        #header_line = f.readline()
        for line in fasta:
            genome_seq.append(line.strip())

    return ''.join(genome_seq)


def get_sequence_from_the_internet(organism_name, email):
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


def island_rule_1_ok(subject):
    '''
    Check if the CpG sites rule #1 is satisfied
    A GC percentage greater than 50% (https://en.wikipedia.org/wiki/CpG_site)
    '''
    nucleotide_seq_str = subject.get_genome()
    gc_percentage_threshold = subject.get_gc_percentage_threshold()

    sample_seq_length = len(nucleotide_seq_str)
    g_nucleotides_count = nucleotide_seq_str.count('G')
    c_nucleotides_count = nucleotide_seq_str.count('C')

    g_percentage = (g_nucleotides_count/sample_seq_length) * 100
    c_percentage = (c_nucleotides_count/sample_seq_length) * 100
    gc_percentage = g_percentage + c_percentage

    if int(gc_percentage) >= int(gc_percentage_threshold):
        return True

    return False


def island_rule_2_ok(subject):
    '''
    Check if the CpG sites rule #2 is satisfied
    An observed-to-expected CpG ratio greater than 60 %
    (https://en.wikipedia.org/wiki/CpG_site)
    '''
    nucleotide_seq_str = subject.get_genome()

    cpg_ratio_threshold = subject.get_cpg_ratio_threshold()

    gc_pairs_obs = nucleotide_seq_str.count('CG')
    gc_pairs_exp = (nucleotide_seq_str.count('C') *
                    nucleotide_seq_str.count('G')) / len(nucleotide_seq_str)
    if gc_pairs_exp == 0:
        return False

    cpg_ratio = (gc_pairs_obs/gc_pairs_exp) * 100

    if int(cpg_ratio) > int(cpg_ratio_threshold):
        return True

    return False


def island_conditions_ok(subject):
    '''
    Check if island conditions are met
    '''
    return island_rule_1_ok(subject) and island_rule_2_ok(subject)


def find_islands(subject):
    '''
    Find CpG islands from the given nucleotide sequence

    '''
    nucleotide_seq = subject.get_genome()
    #print("DADADADAAADAAA 1", nucleotide_seq)
    sample_seq_length = int(subject.get_window_size())
    #print("DADADADAAADAAA 2", sample_seq_length)
    start_index = 0
    end_index = sample_seq_length
    sample_seq = nucleotide_seq[start_index:end_index]
    cpg_islands_list = []
    genome_length = subject.get_genome_length()
    # let's put the original genome sequence first
    cpg_islands_list.append(nucleotide_seq)

    #print('genome_length: ', genome_length)

    while True:
        # this could take the subject
        if island_conditions_ok(subject):
            end_index = end_index + 1
            if end_index > genome_length:
                cpg_islands_list.append((start_index, end_index-1))
                print('break 1:', end_index)
                break

            sample_seq = nucleotide_seq[start_index:end_index]
        else:
            if len(sample_seq) == sample_seq_length:
                # no island in default sample, lets move on
                start_index = start_index + 1
                end_index = start_index + sample_seq_length
                if end_index > genome_length:
                    print('break 2:', end_index)
                    break

                sample_seq = nucleotide_seq[start_index:end_index]
            else:
                cpg_islands_list.append((start_index, end_index))
                start_index = start_index + 1 # move start index by one
                end_index = start_index + sample_seq_length
                if end_index > genome_length:
                    print('break 3:', end_index)
                    break
                sample_seq = nucleotide_seq[start_index:end_index]

    subject.set_cpg_islands(cpg_islands_list)


def create_study_subjects(window_size, gc_percentage, cgp_ratio):
    ''' Create genomes  '''
    subjects = []

    #print('Anna tutkittavan genomin nimi, syötä nimeksi "test" jos haluat testata SH1 genomilla.')

    while 1:
        subject = Genome()
        subject.populate_subject_from_input()

        if subject.name == 'Q':
            break

        subject.set_cpg_ratio_threshold(cgp_ratio)
        subject.set_gc_percentage_threshold(gc_percentage)
        subject.set_window_size(window_size)

        subjects.append(subject)

    return subjects


def search_islands(subjects):
    '''
    Gets list of organisms whos genomes are studied for islands

    '''
    for subject in subjects:
        print('')
        print('')
        print('*** Etsitään saarekkeita ', subject.get_name(), ' genomille')
        print('')
        find_islands(subject)


def visualize_results(subject_list):
    ''' Visualize results  '''
    for subject in subject_list:
        print('visualize_results NOT IMPLEMENTED YET')


def compare_results(subject_list):
    ''' Compare results, displays information of each genome '''
    for subject in subject_list:
        subject.display(True)


def start():
    '''
    The main thing. Asks user input.

    '''
    subjects = []

    print('')
    print('--------------------------------------------------------------------')
    print('Anna saarekkeiden etsintäparametrit, enterillä suoraan default arvot')
    print('')
    sample_window_size = input('  Koeikkunan aloituskoko (default 200):') or 200
    gc_percentage_threshold = input('  GC-pitoisuuden raja-arvo s(default 50):') or 50
    cpg_ratio_threshold = input('  CpG-suhteen raja-arvo (default 60):') or 60
    print('')
    print('--------------------------------------------------------------------')
    print('Seuraavaksi tutkittavat genomit.')
    print('Nimellä "test" testataan SH1_genome.fasta tiedostolla')
    print('')
    subjects = create_study_subjects(
        sample_window_size, gc_percentage_threshold, cpg_ratio_threshold)

    search_islands(subjects)

    compare_results(subjects)

    visualize_results(subjects)


# guard to only execute code when a file is invoked as a script
if __name__ == '__main__':
    start()
