from random import randint
from DNA_data import Data


class Tokenizer:

    def __init__(self, filename, window_size, lower=3, upper=8):
        self.lower = lower
        self.upper = upper
        self.window_size = window_size
        self.counter = window_size

        self.f = open(filename, "r")
        self.f.readline()

    def generate_kmers(self):
        seq = ""

        for line in self.f:
            seq += line[:-1]

        return self.string_to_kmers(seq)

    def string_to_kmers(self, sequence):
        if not sequence:
            return []

        start = 0
        end = 0
        length = len(sequence)
        k_mers = []

        while end < length:
            r = randint(self.lower, self.upper)
            end = start+r

            if end > length:
                end = length

            dna_code = Data.index_of_kmer[sequence[start:end]]
            k_mers.append(dna_code)

            if dna_code in Data.kmer_count:
                Data.kmer_count[dna_code] += 1
            else:
                Data.kmer_count[dna_code] = 1

            start += 1

        return k_mers
