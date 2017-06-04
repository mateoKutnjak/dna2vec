import numpy as np
from Bio import pairwise2


def create_vocabulary():
    print "Creating vocabularies..."

    vocabulary = dict()
    reversed_vocabulary = dict()
    counter = 0

    for len in range(3, 8+1):
        for i in range(0, 4**len):
            kmer = int_to_dna(i, len)

            vocabulary[kmer] = counter
            reversed_vocabulary[counter] = kmer

            counter += 1
    return vocabulary, reversed_vocabulary


def get_dimension():
    return 4**3 + 4**4 + 4**5 + 4**6 + 4**7 + 4**8


def int_to_dna(value, length):
    numerals = []

    for i in range(0, length):
        value, r = divmod(value, 4)
        numerals.append(str(r))

    numstr = ''.join(reversed(numerals))
    return numstr.\
        replace("0", "A").\
        replace("1", "C").\
        replace("2", "G").\
        replace("3", "T")


def softmax_list(l):
    for i in range(0, len(l)):
        l[i] = softmax(l[i])
    return l


def softmax(vector):
    e_x = np.exp(vector - np.max(vector))
    return e_x / e_x.sum()


def sigmoid(x):
    if x > 6:
        return 1.0
    elif x < -6:
        return 0.0
    else:
        return 1 / (1 + np.exp(-x))


def random_matrix(rows, columns):
    return np.matrix(np.random.rand(rows, columns))


def range_for_kmer_length(length):
    if length < 3 or length > 8:
        raise Exception("Invalid kmer length")

    cur_len = 0
    for i in range(3, 8+1):
        if i == length:
            return cur_len, cur_len + 4**i
        cur_len += 4**i


def sim(v, w):
    return (np.dot(v.transpose(), w) / np.linalg.norm(v) / np.linalg.norm(w)).item(0, 0)


def Needleman_Wunsch_score(kmer1, kmer2):
    alignments = pairwise2.align.globalxx(kmer1, kmer2)
    result = pairwise2.format_alignment(*alignments[0])

    return int(result.split("=")[-1][:-1])