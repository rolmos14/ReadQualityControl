import gzip


class FASTQRead:

    def __init__(self, identifier: str, features: str, sequence: str, repetition: str, quality: str):
        self.identifier = identifier
        self.features = features
        self.sequence = sequence
        self.repetition = repetition
        self.quality = quality

    def __str__(self):
        read = self.identifier + " " + self.features + '\n'
        read += self.sequence + '\n'
        read += self.repetition + '\n'
        read += self.quality
        return read

    def get_length(self):
        return len(self.sequence)

    def get_average_gc(self):
        return round((self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence) * 100, 2)

    def get_num_n(self):
        return self.sequence.count('N')

    def get_n_ratio(self):
        return self.get_num_n() / self.get_length() * 100

    def __eq__(self, other):
        return self.sequence == other.sequence

    def __hash__(self):
        return hash(self.sequence)


class FASTQFile:

    def __init__(self, file_lines):
        self.reads = []
        for line in range(0, len(file_lines), 4):
            ident = file_lines[line].split()[0]
            # Features field can be empty
            feat = ""
            try:
                feat = file_lines[line].split()[1] + " " + file_lines[line].split()[2]
            except IndexError:
                pass
            seq = file_lines[line + 1]
            rep = file_lines[line + 2]
            qual = file_lines[line + 3]
            new_read = FASTQRead(ident, feat, seq, rep, qual)
            self.reads.append(new_read)

    def get_read(self, number):
        return self.reads[number]

    def get_num_reads(self):
        return len(self.reads)

    def get_read_avg_length(self):
        return round(sum([read.get_length() for read in self.reads]) / self.get_num_reads())

    def get_read_avg_gc(self):
        return round(sum([read.get_average_gc() for read in self.reads]) / self.get_num_reads(), 2)

    def get_total_repeats(self):
        reads_set = set(self.reads)
        return self.get_num_reads() - len(reads_set)

    def get_num_reads_with_n(self):
        return len([read for read in self.reads if read.get_num_n() > 0])

    def get_read_n_ratio(self):
        return round(sum([read.get_n_ratio() for read in self.reads]) / self.get_num_reads(), 2)

    def show_stats(self):
        # read_lengths = [read.get_length() for read in self.reads]
        # Dictionary made of {length: quantity} pairs
        # length_quantity_dict = {read_length: read_lengths.count(read_length) for read_length in read_lengths}
        # sorted_lengths = sorted(length_quantity_dict.keys())
        # sorted_len_qty_dict = {length: length_quantity_dict[length] for length in sorted_lengths}
        stats = f"Reads in the file = {self.get_num_reads()}\n"
        # Stage 2 statistics:
        # for length, quantity in sorted_len_qty_dict.items():
        #     stats += f"\twith length {length} = {quantity}\n"
        stats += f"Reads sequence average length = {self.get_read_avg_length()}"
        stats += f"\n\nRepeats = {self.get_total_repeats()}"
        stats += f"\nReads with Ns = {self.get_num_reads_with_n()}"
        stats += f"\n\nGC content average = {self.get_read_avg_gc()}%"
        stats += f"\nNs per read sequence = {self.get_read_n_ratio()}%"
        print(stats)


class ReadQualityControl:

    def __init__(self, paths):
        self.files = []
        for path in paths:
            self.files.append(FASTQFile(self.read_file_lines(path)))

    @staticmethod
    def read_file_lines(path):
        with gzip.open(path, 'rt') as file:
            file_lines = file.read().splitlines()
        return file_lines

    def choose_best_file(self):
        best_file = self.files[0]
        for file in self.files[1:]:
            if file.get_total_repeats() + file.get_num_reads_with_n() < \
                    best_file.get_total_repeats() + best_file.get_num_reads_with_n():
                best_file = file
        best_file.show_stats()


fastq_file_paths = [input() for _ in range(3)]
read_quality_control = ReadQualityControl(fastq_file_paths)
read_quality_control.choose_best_file()
