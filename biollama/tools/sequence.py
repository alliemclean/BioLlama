import random
import requests
import xml.etree.ElementTree as ET
from collections import Counter


class Sequence(object):
    def __init__(self, length=50, alphabet=['A', 'C', 'T', 'G']):
        self.length = length
        self.alphabet = alphabet
        self.sequence = self._generate_random()

    def _generate_random(self):
        seq = ''
        for i in range(self.length):
            base = random.choice(self.alphabet)
            seq = seq + base
        return seq

    def region(self, region):
        """ get dna sequence from genomic location from UCSC """
        url = 'http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment='
        if '-' in region:
            region.replace('-', ',')
        response = requests.get("{}{}".format(url, region))
        try:
            et = ET.fromstring(response.content)
            dna = et.getchildren()[0].getchildren()[0].text  # DNA
        except ET.ParseError:
            print('Sequence not found.  {}\n{}'.format(region, response.content))
            dna = ''
        self.sequence = dna.replace('\n', '')
        self.length = len(self.sequence)
        self.alphabet = set([i for i in self.sequence])
        return self.sequence

    def gc_content(self, bases='CGcg'):
        """ GC content of sequence """
        if not len(self.sequence):
            return 0
        ct = Counter([i for i in self.sequence])
        total = 0
        for i in bases:
            total += ct[i]
        return float(total)/sum(ct.values())

    def calculate_tm(self, bases='CGcg'):
        """ Melting temperature for sequence
           https://primerdigital.com/fastpcr/m7.html
           Marmur and Doty equation for long duplexes
        """
        if not len(self.sequence):
            return 0
        ct = Counter([i for i in self.sequence])
        total = 0
        for i in bases:
            total += ct[i]
        return (41 * float(total) - 650)/sum(ct.values()) + 69.3

    def __str__(self):
        return self.sequence
