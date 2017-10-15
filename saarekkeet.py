'''
Python -ohjelmointia biotieteilijöille harjoitustyö 2017.


Etsitään mahdollisia CpG-saarekkeita annetuista genomeista.

CpG-saarekkeet ovat genomissa sellaisia alueita, joissa CG-nukleotidiparien 
osuus on huomattavan suuri verrattuna C- ja G-nukleotidien kokonaismäärään.
CpG-saareke ei ole kuitenkaan yksikäsitteisesti määritelty.
(https://en.wikipedia.org/wiki/CpG_site).


@author Toni Räsänen

'''

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Bio.Entrez
import Bio.SeqIO
import numpy as np
import matplotlib.pyplot as plt



class ResearchSubject:
    '''
    Luokka pyrkii kokoamaan yhteen tietoja tutkittavasta genomista. 
    
    Luokka sisältää genomin nimen, sekvenssin ja listan mahdollisista löydetyistä 
    CpG-saarekkeista. Luokassa on sähköpostiosoite jota käytetään genomin sekvenssin
    hakemisessa internet-palvelusta. Luokka sisältää myös tiedon saarekkeiden 
    haussa käytetyistä raja-arvoista; näyteikkunan aloituskoko (window_size),
    CG-pitoisuus (gc_percentage_threshold) ja (“Obs/Exp”) CpG-parien osuus
    (gc_percentage_threshold).
    
    Luokalla on metodi display jolla se tulostaa genomin tiedot. Tätä käytetään
    vertailtaessa tuloksia saareke-etsinnän lopuksi. 
    
    Luokalla on metodi populate_subject_from_input, jota käytetään 
    genomi-instanssien luomiseen.
    
    '''    
    
    genome_name = ''
    genome_sequence = ''
    researcher_email = ''
    window_size = 0
    cpg_ratio_threshold = 0
    cg_percentage_threshold = 0
    possible_cpg_islands = []
    

    def populate_subject_from_input(self):
        '''
        Käytetään genomi-instanssien luomiseen. 
        
        Kysytään käyttäjältä tietoja tutkittavasta genomista. 
        Haetaan genomin tunnisteen perusteella genomille sekvenssi internet-
        palvelusta. Jos genominen nimeksi annetaan "test" luetaan genomille 
        sekvenssi testitiedostosta (SH1_genome.fasta).
        
        '''
        self.genome_name = input('  Tutkittavan genomin nimi (tai Q niinkuin quit): ')

        if self.genome_name == 'Q':
            return
        if self.genome_name == '' or self.genome_name is None:
            raise Exception('Genomin nimi tarvitaan!')

        # get the genomes
        if self.genome_name == 'test':
            test_fasta_file = 'SH1_genome.fasta'
            self.genome_sequence = read_fasta_file_from_filesystem(test_fasta_file)
            self.genome_name = 'SH1'
        else:
            self.researcher_email = input('  Anna sähköposti: ')
            
            if self.researcher_email is None or self.researcher_email == '':
                raise Exception('Email tarvitaan!')            
            
            self.genome_sequence = get_sequence_from_the_internet(
                self.genome_name, self.researcher_email)


    def display(self, verbose):
        ''' 
        Tulostetaan genomin tiedot konsoliin. Käytetään vertailtaessa tuloksia 
        saareke-etsinnän lopuksi. 
        
        Tulostaa genomin tunnisteen (nimen), genomin sekvenssin pituuden,
        mahdollisten löydettyjen saarekkeiden lukumäärän. 
        
        Listaa lisäksi kaikki löydetyt saarekkeet ja tulostaa niiden pituuden,
        saarekkeen alku- ja loppukohdan (indeksin) genomin sekvenssissä sekä
        saarekkeelle lasketut arvot CG-pitoisuudesta ja CpG-parien osuudesta.
        
        '''
        print('')
        print('Genomi:')
        print('- nimi:', self.genome_name)
        print('- sekvenssin pituus:', self.get_genome_length())
        print('- mahdollisten CpG-saarekkeiden lkm:', len(self.possible_cpg_islands))
        if verbose:
            index = 1
            for island in self.possible_cpg_islands:
                start_idx, end_idx = island
                print('   saareke', index, ': pituus', (end_idx - start_idx), 
                      'start_index:', start_idx, 'end_index:', end_idx, 
                      (calculate_details(self.genome_sequence, start_idx, end_idx)))
                index = index + 1
        print('')

    def get_name(self):
        ''' Palauttaa genomine nimen (tai tunnuksen) '''
        return self.genome_name

    def get_genome(self):
        ''' Palauttaa genomin nukleotidisekvenssin '''
        return self.genome_sequence

    def set_window_size(self, size):
        ''' Asettaa saarekkeiden etsinnässä käytetyn näyteikkunan aloituskoon '''
        self.window_size = size

    def get_window_size(self):
        ''' Palauttaa näyteikkunan koon '''
        return self.window_size

    def set_cpg_ratio_threshold(self, cpg_ratio):
        ''' Asetetaan CpG-parien suhdeluku '''
        self.cpg_ratio_threshold = cpg_ratio

    def get_cpg_ratio_threshold(self):
        ''' Palautetaan CpG-parien suhdeluku '''
        return self.cpg_ratio_threshold

    def set_cg_percentage_threshold(self, cg_percentage):
        ''' Asetetaan käytetty CG-pitoisuuden raja-arvo '''
        self.cg_percentage_threshold = cg_percentage

    def get_cg_percentage_threshold(self):
        ''' Asetetaan käytetty CG-pitoisuuden raja-arvo '''
        return self.cg_percentage_threshold

    def get_genome_length(self):
        ''' Palauttaa genomin nukleotidisekvenssin pituuden '''
        return len(self.genome_sequence)

    def set_cpg_island(self, island):
        ''' Lisää saarekkeen saarekelistaan '''
        self.possible_cpg_islands.append(island)

    def set_cpg_islands(self, islands):
        ''' Alustaa genomin saarekelistan '''
        self.possible_cpg_islands = islands

    def get_cpg_islands(self):
        ''' Palauttaa genomin saarekelistan '''
        return self.possible_cpg_islands


def calculate_details(seq, region_start, region_end):
    '''
    Lasketaan löydetylle saarekkeelle CG-pitoisuus (“CG%”) ja CpG-parien osuus
    (“Obs/Exp”).
    Palautetaan arvot merkkijonona.
    
    Tietoja käytetään ja ne näytetään löydettyjen saarekkeiden vertaiuosiossa
    jossa löydetyt saarekkeet listatataan.
    
    '''
    
    def obs(nucleotide_seq_str):
        cg_pairs_obs = nucleotide_seq_str.count('CG')
        cg_pairs_exp = (nucleotide_seq_str.count('C') * 
                        nucleotide_seq_str.count('G')) / len(nucleotide_seq_str)

        cpg_ratio = (cg_pairs_obs/cg_pairs_exp) * 100
        return cpg_ratio
    
    def cg(nucleotide_seq_str):
        sample_seq_length = len(nucleotide_seq_str)
        g_nucleotides_count = nucleotide_seq_str.count('G')
        c_nucleotides_count = nucleotide_seq_str.count('C')
    
        g_percentage = (g_nucleotides_count/sample_seq_length) * 100
        c_percentage = (c_nucleotides_count/sample_seq_length) * 100
        cg_percentage = g_percentage + c_percentage
        return cg_percentage
    
    region = seq[region_start:region_end] # island sequence
    obs = round(obs(region), 2)
    cg = round(cg(region), 2)
    
    return ('(Obs/Exp = ' + str(obs) + ' ja %CG = ' + str(cg) + ')')


def read_fasta_file_from_filesystem(file_name):
    '''
    Luetaan genomin sekvenssi tiedostosta ja palautetaan se.
    
    Tätä funktiota käytetään testisekvenssin lukemiseen tiedostosta.
    
    '''
    genome_seq = []

    with open(file_name) as fasta:
        for line in fasta:
            genome_seq.append(line.strip())

    return ''.join(genome_seq)


def get_sequence_from_the_internet(organism_name, email):
    '''
    Haetaan annetulla genomin nimellä genomille sekvenssi internetistä
    Bio.Entrez palvelusta. Palvelun käyttöön tarvitaan sähköpostiosoite.
    
    Palautetaan löytynyt sekvenssi. Jos palvelu palauttaa virheen sitä ei
    käsitellä vaan se menee eteenpäin käyttäjälle.
    
    (http://biopython.org/DIST/docs/api/Bio.Entrez-module.html)

    '''
    Bio.Entrez.email = email
    handle = Bio.Entrez.efetch(db="nucleotide", rettype="gb", retmode="text",
                               id=organism_name)
    seq_record = Bio.SeqIO.read(handle, "gb")
    handle.close()

    return seq_record.seq


def island_rule_1_ok(subject, seq):
    '''
    Tarkistetaan ylittyykö CG-pitoisuus annetulla sekvenssialueella (seq).
    
    CG-pitoisuus (“CG%”) on vähintään 50%. CG-pitoisuus lasketaan G- ja 
    C-nukleotidien yhteisenä prosenttiosuutena (C% + G%) ko. alueella.
    
    Raja-arvoa on mahdollista muuttaa ohjelman käynnistyksen yhteydessä.
    
    '''
    nucleotide_seq_str = seq
    cg_percentage_threshold = subject.get_cg_percentage_threshold()

    sample_seq_length = len(nucleotide_seq_str)
    g_nucleotides_count = nucleotide_seq_str.count('G')
    c_nucleotides_count = nucleotide_seq_str.count('C')

    g_percentage = (g_nucleotides_count/sample_seq_length) * 100
    c_percentage = (c_nucleotides_count/sample_seq_length) * 100
    cg_percentage = g_percentage + c_percentage

    if int(cg_percentage) >= int(cg_percentage_threshold):
        return True

    return False


def island_rule_2_ok(subject, seq):
    '''
    Tarkistetaan onko havannoitujen-vs-oletusarvoisten (“Obs/Exp”) CpG-parien 
    osuus yli 0.6 (60%) annetulla sekvenssialueella (seq)
    
    
    Havainnoitujen (Obs) CpG-parien määrä on yksinkertaisesti välittömästi 
    toisiaan seuraavien C- ja G-nukleotidien määrä ko. alueella.

    Oletusarvo CpG-pareille (Exp) voidaan arvioida kaikkien C- ja G-nukleotidien 
    määrästä ko. alueella kertomalla löydettyjen C-nukleotidien määrä 
    G-nukleotidien määrällä ja normalisoimalla tulo alueen pituudella: 
    Exp = (C:n määrä * G:n määrä) / alueen_pituus.
    
    Raja-arvoa on mahdollista muuttaa ohjelman käynnistyksen yhteydessä.
    
    '''
    nucleotide_seq_str = seq

    cpg_ratio_threshold = subject.get_cpg_ratio_threshold()

    cg_pairs_obs = nucleotide_seq_str.count('CG')
    cg_pairs_exp = (nucleotide_seq_str.count('C') *
                    nucleotide_seq_str.count('G')) / len(nucleotide_seq_str)
    if cg_pairs_exp == 0:
        return False

    cpg_ratio = (cg_pairs_obs/cg_pairs_exp) * 100

    if int(cpg_ratio) > int(cpg_ratio_threshold):
        return True

    return False


def island_conditions_ok(subject, seq):
    '''
    Tarkistetaan toteutuuko saarekkeille asetetut vaatimukset annetulla
    sekvenssialuella (seq)
    '''
    return island_rule_1_ok(subject, seq) and island_rule_2_ok(subject, seq)


def find_islands(subject):
    '''
    Etsitään CpG-saarekkeita genomin sekvenssistä käymällä sitä läpi näyteikkunan
    kokoisina paloina kasvattaen näyteikkunaa kunnes saarekkeen ehdot eivät
    enää täyty. Siirretään näyteikkunan aloitussekvenssiä yhdellä kunnes koko
    genomi on käyty läpi.

    '''
    # genomin nukleotidisekvenssi
    nucleotide_seq = subject.get_genome()
    # näytealueen aloituskoko
    seq_sample_length = int(subject.get_window_size())
    # näytteen alkuindeksin alkuarvo
    start_index = 0
    # näytteen loppuindeksin alkuarvo
    end_index = seq_sample_length
    # sekvenssin pala
    seq_sample = nucleotide_seq[start_index:end_index]
    # saarekelista
    cpg_islands_list = []

    while True:
        
        # otetaan näyte
        seq_sample = nucleotide_seq[start_index:end_index]
        
        # tarkistetaan että näyte on minimipituinen
        if not len(seq_sample) >= seq_sample_length:
            break
        
        # tarkistetaan täyttääkö näyte saarekkeen ehdot
        if island_conditions_ok(subject, seq_sample):
            # jos täyttää siirretään loppuindeksiä yhdellä eteenpäin
            end_index = end_index + 1
            # tarkistetaan ollaanko genomin sekvenssin lopussa 
            if end_index > subject.get_genome_length():
                # lisätään saareke listaan
                cpg_islands_list.append((start_index, end_index-1))
                # siirretään alkuindeksiä yhdellä, samoin loppuindeksiä
                start_index = start_index + 1
                end_index = start_index + seq_sample_length
        else:
            # saarekkeen ehdot eivät täyttyneet
            # jos näyte on suurempi kuin näytealueen oletuskoko, otetaan se talteen
            if len(seq_sample) > seq_sample_length:
                cpg_islands_list.append((start_index, end_index-1))
            # kasvatetaan alku- ja loppuindeksejä
            start_index = start_index + 1
            end_index = start_index + seq_sample_length
            # tarkistetaan onko genomin sekvenssiä vielä jäljellä
            if end_index > subject.get_genome_length():
                break

    subject.set_cpg_islands(cpg_islands_list)


def create_study_subjects(window_size, gc_percentage, cgp_ratio):
    ''' 
    Luodaan, instantioidaan, tutkittavat genomit
    
    '''
    subjects = []

    while 1:
        subject = ResearchSubject()
        subject.populate_subject_from_input()

        if subject.genome_name == 'Q':
            break

        subject.set_cpg_ratio_threshold(cgp_ratio)
        subject.set_cg_percentage_threshold(gc_percentage)
        subject.set_window_size(window_size)

        subjects.append(subject)

    return subjects


def search_islands(subjects):
    '''
    Käydään läpi tutkittavat genomit ja etsitään saarekkeita kustakin.

    '''
    print('')
    print('*** Etsitään saarekkeita ***')
    print('')
    for subject in subjects:
        find_islands(subject)


def visualize_results(subject_list):
    ''' 
    Näytetään kuvaajat tutkittujen genomien löydöksistä.
     
    '''
    for subject in subject_list:
        bar_chart(subject)


def bar_chart(subject):
    ''' 
    Näytetään genomin mahdolliset saarekkeet pylväsdiagrammiesityksenä.

    Saarekkeet ovat rivissä x-akselilla ja pylväiden pituus kuvaa saarekkeen
    kokoa, saarekkeen pituus näkyy siis y-akselilla.
    
    Kuvaajan heikkoutena on mm. se että kuvaaja menee nopeasti tukkoon jos 
    genomista löytyy paljon saarekkeita.
    
    '''
    N = len(subject.get_cpg_islands())
    
    island_lengths = []
    island_labels = []
    for island in subject.get_cpg_islands():
        start_idx, end_idx = island
        island_lengths.append((end_idx - start_idx))
        island_labels.append(str(start_idx) + '-' + str(end_idx))
    
    ind = np.arange(N) 
    width = 0.35       
    
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, island_lengths, width, color='y')
    
    ax.set_ylabel('Saarekkeen pituus')
    ax.set_xlabel('Saarekkeet, ' + str(N) + (' kappale' if N == 1 else ' kappaletta' ))
    ax.set_title( subject.genome_name + ' genomin mahdolliset CpG-saarekkeet')
    ax.set_xticks(ind)
    ax.set_xticklabels(island_labels)
    
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 1.05 * height,
                    '%d' % int(height),
                    ha='center', va='bottom')
    
    autolabel(rects1)
    
    plt.show()


def compare_results(subject_list):
    ''' Compare results, displays information of each genome '''
    print('')
    print('*** Saarekkeiden etsinnän tulokset ***')
    for subject in subject_list:
        subject.display(True)
    print('')   


def start():
    '''
    Ohjelman aloitus. 
    Kysytään käyttäjältä saarekkeiden etsinnässä käytetyt raja-arvot, 
    käyttäjä voi enterillä hyväksyä oletusarvot.  
    
    Kysytään käyttäjältä tutkittavien genomien nimet, tunnisteet ja luodaan 
    tutkittavat genomit (create_study_subjects).
    
    Etsitään saarekkeita (search_islands).
    
    Tulostetaan etsinnäin tulokset kultakin genomilta (compare_results).
    
    Näytetään kuvaaja kullekin genomille saareke-etsinnäin tuloksista 
    (visualize_results).

    '''
    subjects = []

    print('')
    sample_window_size = input('  Koeikkunan aloituskoko (default 200):') or 200
    cg_percentage_threshold = input('  CG-pitoisuuden raja-arvo s(default 50):') or 50
    cpg_ratio_threshold = input('  CpG-suhteen raja-arvo (default 60):') or 60
    print('')
        
    print('--------------------------------------------------------------------')
    print('Seuraavaksi tutkittavat genomit.')
    print('Nimellä "test" käytetään testitiedostoa SH1_genome.fasta')
    print('')
    subjects = create_study_subjects(
        sample_window_size, cg_percentage_threshold, cpg_ratio_threshold)

    search_islands(subjects)

    compare_results(subjects)
    print('')
    print('Saarekkeiden hakuun käytetyt arvot:')
    print('- hakuikkunan aloituskoko:', sample_window_size)
    print('- CpG-suhteen raja-arvo:', cpg_ratio_threshold)
    print('- CG-pitoisuuden raja-arvo:', cg_percentage_threshold)
    print('')

    visualize_results(subjects)


# guard to only execute code when a file is invoked as a script
if __name__ == '__main__':
    start()
