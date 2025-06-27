from Bio import SeqIO
from Bio.Seq import Seq
import sys

class FastaAnalyzer:
    """
    Analyze fasta file using biopython
    """
    def __init__(self, fasta_file: str):
        """
        Opens the fasta file in read mode
        :param fasta_file:
        """
        self.handle = None
        self.seqs = {}
        name = None
        if fasta_file.endswith(".gz"):
            try:
                self.handle = open(fasta_file, 'rt')
            except IOError as e:
                sys.stderr.write(f"Error opening {fasta_file}: {e}")

        for line in self.handle:
            line = line.strip()
            if line.startswith(">"):
                words=line[1:].split()
                name=words [0][1:]
                self.seqs[name]=''
            else:
                self.seqs[name]= self.seqs[name] + line.lower()

    def __del__(self):
        self.handle.close()
        pass

    def how_many_records(self):
        records = len(self.seqs)
        return records

    def sequence_lengths(self):
        lengths = {}
        longest = {}
        for name, seq in self.seqs.items():
            lengths[name] = len(seq)
            if longest[name] < len(seq):
                longest[name] = len(seq)
        print(longest)
        return lengths

    # noinspection SpellCheckingInspection
    def reading_frame(self):
        reading_frames = {}
        for name,seq in self.seqs.items():
           bioseq = []
           for i in range(3):
                seqt = seq[i:]
                bioseq.append(seqt)
           reading_frames[name] = bioseq
        return reading_frames

    def orf(self, reading_frames :dict = None):
        if reading_frames:
            for names, rfs in reading_frames.items():
                if len(reading_frames[rfs]) != 3:
                    print("Reading frame not valid")
                    return None
        else:
            reading_frames = self.reading_frame()
        open_reading_frames = {}
        longest = {}
        stop_codons = ['tga', 'tag', 'tga']
        start_codon = 'atg'
        for name, rfs in reading_frames.items():
            for i in range(3):
                orfs = {}
                seq = rfs[i]
                rfs[i] = ''
                orf = ''
                start = len(seq) + 1
                for j in range(len(seq),3):
                    codon = seq[j:j+3]
                    if codon == start_codon:
                        start = j+1
                    if start < j:
                        orf = orf + codon
                    if codon in stop_codons:
                        orfs[start] = orf + codon
                        if longest[rfs] < len(orf):
                            longest[rfs] = len(orf)
                        orf = ''
                        start = len(seq)+1
                        rfs[i] = orfs
                        open_reading_frames[name] = rfs[i]
        print(longest)
        return open_reading_frames

    def longest_orf_of_id(self, id : str, open_reading_frames:dict = None):
        longest = {}
        if open_reading_frames:
            for names, rfs in open_reading_frames.items():
                for i in range(3):
                    for j in range(len(rfs[i])):
                        if rfs[i][j][0:3] != 'atg':
                            print('Open reading frame not valid')
                            return None
        else:
            open_reading_frames = self.orf()
        # noinspection PyInconsistentReturns
        for names, rfs in open_reading_frames.items():
            if names == id:
                for i in range(3):
                    for j in range(len(rfs[i])):
                        if longest[names][rfs[i][j]] < len(rfs[i][j]):
                            longest[names][rfs[i]] = len(rfs[i][j])
                return longest
            else:
                print('id not found')
                return None






