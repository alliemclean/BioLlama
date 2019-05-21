import random


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

    def __str__(self):
        return self.sequence
