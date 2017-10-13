'''
Search for CpG islands
'''

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Bio.Entrez
import Bio.SeqIO


class ResearchSubject:
    '''
    Genome under research
    '''
    name = ''
    genome_sequence = ''
    cpg_islands = [] # list holding indexes of found islands and maybe the calculated values of
                     # cpg_ratio and gc_percentage per island
    cpg_ratio_threshold = 0
    gc_percentage_threshold = 0
    window_size = 0
    researcher_email = ''

    def populate_subject_from_input(self):
        '''
        Populate genome from input
        '''
        self.name = input('Anna tutkittavan eliön nimi (tai Q niinkuin quit): ')

        if self.name == 'Q':
            return
        if self.name == '' or self.name is None:
            raise Exception('Eliön nimi tarvitaan!')

        self.researcher_email = input('Anna sähköposti: ')

        if self.researcher_email is None or self.researcher_email == '':
            raise Exception('Email tarvitaan!')

        # get the genomes
        if self.name == 'test':
            self.genome_sequence = self.__read_test_fasta_file_from_filesystem()
            self.name = 'SH1' # this is the name of the test subject
        else:
            self.genome_sequence = self.__read_fasta_file_from_the_internet(
                self.name, self.researcher_email)


    def get_genome(self):
        ''' Return genome sequence string '''
        return self.genome_sequence


    def set_window_size(self, size):
        ''' Set the sample window size '''
        self.window_size = size


    def set_cpg_ratio_threshold(self, cpg_ratio):
        ''' Set the CpG ratio '''
        self.cpg_ratio_threshold = cpg_ratio


    def set_gc_percentage_threshold(self, gc_percentage):
        ''' Set the GC percentage '''
        self.gc_percentage_threshold = gc_percentage


    def get_cpg_ratio_threshold(self):
        ''' Return Cpg ratio '''
        return self.cpg_ratio_threshold


    def get_gc_percentage_threshold(self):
        ''' Return GC percentage value '''
        return self.gc_percentage_threshold

    def get_window_size(self):
        ''' Return the sample window size '''
        return self.window_size

    def get_genome_length(self):
        ''' Return the length of the genome sequence '''
        return len(self.genome_sequence)

    def set_cpg_islands(self, islands):
        ''' Set CpG islands '''
        self.cpg_islands = islands

    def get_cpg_islands(self):
        ''' Return Cpg islands list '''
        return self.cpg_islands


    def __read_test_fasta_file_from_filesystem(self):
        '''
        Mainly for testing purposes.
        '''
        dummy_test_seq = 'CTGGACACCAGCGTAGACCTGCGGTTCAAGTGACCATGCCGGGAATCGTCTCACAGTACGTGCTCCCCGT'

        test_fasta_file = 'SH1_genome.fasta'
        genome_seq = []

        with open(test_fasta_file) as fasta:
            #header_line = f.readline()
            for line in fasta:
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
    nucleotide_seq = subject.get_genome()
    print("DADADADAAADAAA 1", nucleotide_seq)
    sample_seq_length = int(subject.get_window_size())
    print("DADADADAAADAAA 2", sample_seq_length)
    start_index = 0
    end_index = sample_seq_length
    sample_seq = nucleotide_seq[start_index:end_index]
    cpg_islands_list = []
    genome_length = subject.get_genome_length()
    # let's put the original genome sequence first
    cpg_islands_list.append(nucleotide_seq)

    print('genome_length: ', genome_length)

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

    #return cpg_islands_list


def create_study_subjects(window_size, gc, cgp):
    ''' Create genomes  '''
    subjects = []

    print('Anna tutkittavan genomin nimi, syötä nimeksi "test" jos haluat testata SH1 genomilla.')

    while 1:
        subject = ResearchSubject()
        subject.populate_subject_from_input()

        if subject.name == 'Q':
            break

        subject.set_cpg_ratio_threshold(cgp)
        subject.set_gc_percentage_threshold(gc)
        subject.set_window_size(window_size)

        subjects.append(subject)

    return subjects


def search_islands(subjects):
    '''
    Gets list of organisms whos genomes are studied for islands

    '''
    for subject in subjects:
        find_islands(subject)


def visualize_results(subject_list):
    ''' Visualize results  '''
    for subject in subject_list:
        print('visualize_results NOT IMPLEMENTED YET', len(subject.get_cpg_islands()))


def compare_results(subject_list):
    ''' Compare results '''
    for subject in subject_list:
        print('compare_results NOT IMPLEMENTED YET', len(subject.get_cpg_islands()))


def start():
    '''
    The main thing. Asks user input.

    '''
    subjects = []

    sample_window_size = input('Anna koeikkunan aloituskoko (default 200):') or 200
    gc_percentage_threshold = input('GC-pitoisuuden raja-arvo (default 50):') or 50
    cpg_ratio_threshold = input('CpG-suhteen raja-arvo (default 60):') or 60
    print('Anna nimeksi "test" jos haluat testata SH1_genome.fasta tiedostolla')
    subjects = create_study_subjects(
        sample_window_size, gc_percentage_threshold, cpg_ratio_threshold)

    search_islands(subjects)

    compare_results(subjects)

    visualize_results(subjects)


# guard to only execute code when a file is invoked as a script
if __name__ == '__main__':
    start()
