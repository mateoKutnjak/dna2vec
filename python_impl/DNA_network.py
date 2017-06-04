import numpy as np
import DNA_util


class Network:

    def __init__(self, hidden_size, context_size, unigram, ns_size=5, alpha=0.05):
        self.learning_rate = alpha

        self.V = DNA_util.get_dimension()
        self.N = hidden_size
        self.C = context_size

        self.theta1 = DNA_util.random_matrix(self.V, self.N)
        self.theta2 = [None] * context_size
        for c in range(0, context_size):
            self.theta2[c] = DNA_util.random_matrix(self.N, self.V)

        self.h = None
        self.u = None
        self.y = None

        self.unigram = unigram

        self.neg_samples = None
        self.pos_sample = 0
        self.ns_size = ns_size

        self.loss = None

    def result(self, input_index):
        self.h = np.matrix(self.theta1[input_index, :]).copy().transpose()
        self.u = [num.transpose() * self.h for num in self.theta2]
        self.y = DNA_util.softmax_list(self.u)
        return np.concatenate(self.y)

    def train(self, input_index, output_indexes):
        self.loss = [0] * self.C

        for c in range(0, self.C):
            output_index = output_indexes[c]
            self.neg_samples = self.unigram.choose_n(self.ns_size, output_index)

            vWI = self.theta1[input_index, :].transpose()
            #vwO = self.theta1[output_index, :].transpose()
            #v_wI = self.theta2[c][:, input_index]
            v_wO = self.theta2[c][:, output_index]

            sigmoid_out = DNA_util.sigmoid((v_wO.transpose() * vWI).item(0, 0))

            gwI = (sigmoid_out-1) * v_wO
            gwO = sigmoid_out-1
            gwj = []

            self.loss[c] -= np.log(sigmoid_out)

            for neg_sample in self.neg_samples:
                v_wj = self.theta2[c][:, neg_sample]

                data_for_sigmoid = (v_wj.transpose() * vWI).item(0, 0)
                sigmoid_wj = DNA_util.sigmoid(data_for_sigmoid)

                gwj.append(sigmoid_wj)

                gwI += sigmoid_wj * v_wj

                self.loss -= np.log(DNA_util.sigmoid(-data_for_sigmoid))

            for i in range(0, len(self.neg_samples)):
                self.theta2[c][:, self.neg_samples[i]] -= self.learning_rate * gwj[i] * vWI

            self.theta2[c][:, output_index] -= self.learning_rate * gwO * vWI
            self.theta1[input_index, :] -= self.learning_rate * gwI.transpose()
