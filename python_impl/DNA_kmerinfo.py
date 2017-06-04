from DNA_data import Data
import DNA_util


class Kmer_info:

    def __init__(self, network):
        self.network = network

    def nearest_neighbor(self, v, neighbor_len):
        if neighbor_len < 3 or neighbor_len > 8:
            raise Exception("k must be in range [3,8]", neighbor_len)

        if type(v) is str:
            v = self.network.forward_propagation(Data.index_of_kmer[v])

        _from, _to = DNA_util.range_for_kmer_length(neighbor_len)
        max_sim, max_index = 0.0, -1

        for i in range(_from, _to):
            output = self.network.forward_propagation(i)
            cur_sim = DNA_util.sim(v, output)

            if cur_sim > max_sim:
                max_sim = cur_sim
                max_index = i

        print max_sim
        return Data.kmer_of_index[max_index]

    def n_nearest_neighbors(self, v, neighbor_len, n):
        if neighbor_len < 3 or neighbor_len > 8:
            raise Exception("k must be in range [3,8]", neighbor_len)

        if n == 1:
            return self.nearest_neighbor(v, neighbor_len)

        if type(v) is str:
            v = self.network.forward_propagation(Data.index_of_kmer[v])

        _from, _to = DNA_util.range_for_kmer_length(neighbor_len)
        l = []

        print _from, _to

        for i in range(_from, _to):
            print i
            output_vec = self.network.forward_propagation(i)
            cur_sim = DNA_util.sim(v, output_vec)

            l.append((cur_sim, Data.kmer_of_index[i]))

        l.sort()
        l = l[::-1]

        print l

        return [elem[1] for elem in l][:n]