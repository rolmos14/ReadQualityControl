class FASTQRead:

    def __init__(self, identifier, features, sequence, repetition, quality):
        self.identifier = identifier;
        self.features = features;
        self.sequence = sequence;
        self.repetition = repetition;
        self.quality = quality;

    def __str__(self):
        read = self.identifier + " " + self.features + '\n'
        read += self.sequence + '\n'
        read += self.repetition + '\n'
        read += self.quality
        return read


class FASTQFile:

    def __init__(self, file_lines):
        self.reads = []
        for line in range(0, len(file_lines), 4):
            ident = file_lines[line].split()[0]
            feat = file_lines[line].split()[1] + " " + file_lines[line].split()[2]
            seq = file_lines[line + 1]
            rep = file_lines[line + 2]
            qual = file_lines[line + 3]
            new_read = FASTQRead(ident, feat, seq, rep, qual)
            self.reads.append(new_read)

    def get_read(self, number):
        return self.reads[number]


class ReadQualityControl:

    def __init__(self, path):
        self.file = FASTQFile(self.read_file_lines(path))

    @staticmethod
    def read_file_lines(path):
        with open(path) as file:
            file_lines = file.readlines()
        return file_lines


fastq_file_path = input()
read_quality_control = ReadQualityControl(fastq_file_path)
print(read_quality_control.file.get_read(0))
