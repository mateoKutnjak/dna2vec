import DNA_util


class Data:

    index_of_kmer, kmer_of_index = DNA_util.create_vocabulary()
    kmer_count = dict()
    kmer_distribution = dict()