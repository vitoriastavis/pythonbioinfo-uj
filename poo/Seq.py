import re

class Seq:

    FASTA_REGEX = '\>([^\ \n]+)\ ?(.*)\n([^\>]+)'

    # constructor
    def __init__(self, seqid, title, seq):
        self.seqid = seqid
        self.title = title
        self.seq = seq

    # split string into lines
    @staticmethod
    def lines(seq, lw = 60):
        sequence = '\n'.join(
            seq[i:i+lw] for i in range (0, len(seq), lw)
        )     

        return sequence


    # return a sequence in fasta format
    def fasta(self, lw = 60):
        seqid = self.seqid
        title = self.title
        seq = self.seq
    
        if title != '': title = f'{title}'
        seq = Seq.lines(seq)

        fasta_seq = f'{seqid}{title}\n{seq}'

        return fasta_seq

    # deserialize sequences from a file
    @classmethod    
    def from_file(cls, filename):
        with open(filename) as f: buf = f.read()

        seqs = {}

        for match in re.finditer(cls.FASTA_REGEX, buf):

            print(match)

            seqid, title, seq = match.groups()

            if seqid in seqs:
                raise Exception(f'Non-unique: {seqid}')

            seq = re.sub('[\ \n\t]+', '', seq)

            seqs[seqid] = cls(seqid, title, seq)

        return seqs

    ## check lines, re methods, 