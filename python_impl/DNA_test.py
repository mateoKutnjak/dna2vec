import pickle
from DNA_kmerinfo import Kmer_info
from DNA_data import Data
import DNA_util
from random import randint


def main():
    network = pickle.load(open("final_network", "rb"))

    dimension = DNA_util.get_dimension()
    print dimension
    exit()
    from_3, to_3 = DNA_util.range_for_kmer_length(3)
    from_4, to_4 = DNA_util.range_for_kmer_length(4)
    from_5, to_5 = DNA_util.range_for_kmer_length(5)
    from_6, to_6 = DNA_util.range_for_kmer_length(6)
    from_7, to_7 = DNA_util.range_for_kmer_length(7)
    from_8, to_8 = DNA_util.range_for_kmer_length(8)

    data = dict()



    while True:
        index1 = randint(from_8, to_8 - 1)
        index2 = randint(from_8, to_8 - 1)

        kmer1 = Data.kmer_of_index[index1]
        kmer2 = Data.kmer_of_index[index2]

        score = DNA_util.Needleman_Wunsch_score(kmer1, kmer2)

        if score in data and data[score][0] >= 100:
            continue
        if score > 2 or score < 1:
            continue

        vec1 = network.result(index1)
        vec2 = network.result(index2)

        sim = DNA_util.sim(vec1, vec2)

        if score not in data:
            data[score] = (1, sim)
        else:
            t = data[score]
            data[score] = (t[0] + 1, t[1] + sim)

        stop = True
        print "(",
        for key, value in data.items():
            print value[0],
            if value[0] < 100:
                stop = False

        print ")"

        if stop:
            break

    for key, value in data.items():
        print key, "->", value[1] / value[0]

if __name__ == "__main__":
    main()