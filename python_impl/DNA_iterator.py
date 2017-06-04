class Iterator:

    def __init__(self, window_size, kmers):
        self.window_size = window_size
        self.counter = window_size

        self.kmers = kmers

    def has_next(self):
        if self.counter < self.window_size:
            return False
        elif self.counter >= len(self.kmers)-self.window_size:
            return False
        return True

    def next_kmer(self):
        window = self.kmers[self.counter - self.window_size:self.counter + self.window_size + 1]
        middle_value = window.pop(self.window_size)

        self.counter += 1

        return middle_value, window

