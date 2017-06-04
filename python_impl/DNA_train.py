from DNA_network import Network
from DNA_tokenizer import Tokenizer
from DNA_unigram import Unigram
from DNA_iterator import Iterator
import pickle
import os.path

WINDOW_SIZE = 5
HIDDEN_LAYER_SIZE = 5

LAST_LOADED_COUNTER = 0

SAVE_SWITCH = 0

SWITCH_FILENAME = "save_switch"
ITERATOR_FILENAME = "save_iterator"
NETWORK_FILENAME = "save_network"

FINAL_NETWORK = "final_network"


def load_iterator():
    if not os.path.exists(ITERATOR_FILENAME + str(SAVE_SWITCH)):
        print "Tokenizing and creating iterator..."
        kmers = Tokenizer("eColi.txt", WINDOW_SIZE).generate_kmers()
        return Iterator(WINDOW_SIZE, kmers)

    print "Loading iterator (switch =", SAVE_SWITCH, ")..."
    return pickle.load(open(ITERATOR_FILENAME + str(SAVE_SWITCH), "rb"))


def load_network():
    if not os.path.exists(NETWORK_FILENAME + str(SAVE_SWITCH)):
        print "Making unigram table and SGNS network..."
        unigram = Unigram()
        return Network(HIDDEN_LAYER_SIZE, WINDOW_SIZE*2, unigram, alpha=0.2)

    print "Loading network (switch =", SAVE_SWITCH, ")..."
    return pickle.load(open(NETWORK_FILENAME + str(SAVE_SWITCH), "rb"))


def load_switch():
    global SAVE_SWITCH

    if not os.path.exists(SWITCH_FILENAME):
        SAVE_SWITCH = 0
    else:
        SAVE_SWITCH = pickle.load(open(SWITCH_FILENAME, "rb"))


def load_data():
    global SAVE_SWITCH
    global LAST_LOADED_COUNTER

    load_switch()

    iterator = load_iterator()
    network = load_network()

    SAVE_SWITCH = (SAVE_SWITCH+1) % 2
    LAST_LOADED_COUNTER = iterator.counter

    return iterator, network


def save_data(iterator, network):
    global SAVE_SWITCH

    print "Saving iterator (switch =", SAVE_SWITCH, ")..."
    pickle.dump(iterator, open(ITERATOR_FILENAME + str(SAVE_SWITCH), "wb"))

    print "Saving network (switch =", SAVE_SWITCH, ")..."
    pickle.dump(network, open(NETWORK_FILENAME + str(SAVE_SWITCH), "wb"))

    pickle.dump(SAVE_SWITCH, open(SWITCH_FILENAME, "wb"))
    SAVE_SWITCH = (SAVE_SWITCH+1) % 2


def main():
    iterator, network, = load_data()

    print "Training network..."

    while iterator.has_next():

        input_index, output_indexes = iterator.next_kmer()
        network.train(input_index, output_indexes)

        if iterator.counter % 100 == 0:
            print "iterations:", iterator.counter, "/", len(iterator.kmers)
            print "loss: ", network.loss
        if iterator.counter % 100000 == 0 and iterator.counter != LAST_LOADED_COUNTER:
            save_data(iterator, network)

    network.neg_samples = None
    network.dEdW = None
    network.unigram = None

    pickle.dump(network, open(FINAL_NETWORK, "wb"))
if __name__ == "__main__":
    main()
