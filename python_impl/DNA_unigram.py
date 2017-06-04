from random import randint
from DNA_data import Data


class Unigram:

    def __init__(self):
        self.table = None
        self.size = 0

        self.count_to_probability()
        self.create_unigram_table()

    def create_unigram_table(self):
        size = 10000000
        self.table = [0] * size

        cur_pos = 0
        end_pos = 0

        for dna_code, probability in Data.kmer_distribution.items():
            end_pos = int(cur_pos + probability * size)

            for i in range(cur_pos, end_pos):
                self.table[i] = dna_code
            cur_pos = end_pos

        for i in range(0, size-end_pos):
            self.table.pop()
        self.size = len(self.table)

    def count_to_probability(self):
        denominator = 0.0

        for dna_code, count in Data.kmer_count.items():
            denominator += count**0.75

        for dna_code, count in Data.kmer_count.items():
            Data.kmer_distribution[dna_code] = count ** 0.75 / denominator

    def choose_n(self, n, forbidden_value):
        values = set()

        while len(values) is not n:
            r = randint(0, self.size-1)
            if self.table[r] is not forbidden_value:
                values.add(self.table[r])

        return list(values)
