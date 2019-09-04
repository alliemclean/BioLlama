import random
import requests
import xml.etree.ElementTree as ET


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
        print(response)
        et = ET.fromstring(response.content)
        dna = et.getchildren()[0].getchildren()[0].text  # DNA
        return dna.replace('\n', '')

    def __str__(self):
        return self.sequence
