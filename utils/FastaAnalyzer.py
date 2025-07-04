from Bio import SeqIO
from Bio.Seq import Seq
import sys
import collections

def phred33_to_q(phred33):
   return ord(phred33) - 33



class FastaAnalyzer:
    """
    Analyze fasta file using biopython
    """
    def __init__(self, fasta_file: str, fastq_file: bool=False):
        """
        Opens the fasta file in read mode
        :param fasta_file:
        """
        self.handle = None
        self.seqs = {}
        self.seq_with_quality = {}
        self.fastq_sequences = []
        self.fastq_qualities = []
        name = None
        try:
            self.handle = open(fasta_file, 'r')
        except IOError as e:
            sys.stderr.write(f"Error opening {fasta_file}: {e}")
        if not fastq_file:
            for line in self.handle:
                line = line.strip()
                if line.startswith(">"):
                    name=line[1:].split()[0]
                    self.seqs[name]=''
                else:
                    self.seqs[name]= self.seqs[name] + line.lower()
        else:
            while True:
                self.handle.readline()
                seq = self.handle.readline().rstrip().lower()
                self.handle.readline()
                qual = self.handle.readline().rstrip().lower()
                if len(seq) ==0:
                    break
                self.fastq_sequences.append(seq)
                self.fastq_qualities.append(qual)

    def get_sequences(self):
        return self.seqs

    def naive(self, p):
        occurrences = {}
        for name,seq in self.seqs.items():
            matches = []
            for i in range( len(seq) - len(p) + 1):
                match = True
                for j in range(len(p)):
                    if seq[i + j] != p[j]:
                        match = False
                        break
                if match:
                    matches.append(i)
                    occurrences[name] = matches
        return occurrences

    def reverse_complement(self):
        complement = { 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n' }
        reverse = {}
        for name,seq in self.seqs.items():
            t = ''
            for base in seq:
                t += complement[base]
            reverse[name] = t
        return reverse

    def fastq_create_hist(self):
        hist = [0] * 50
        for qual in self.fastq_qualities:
            for phred in qual:
                q = phred33_to_q(phred)
                hist[q] += 1
        return hist

    def nucleotide_counts(self):
        count = collections.Counter()
        if self.fastq_sequences:
           for seq in self.fastq_sequences:
               count.update(seq)
           return count
        else:
            for name, seq in self.seqs.items():
                count.update(seq)
            return count

    def fastq_find_gc_by_pos(self):
        gc = [0] * 100
        totals = [0] * 100
        for read in self.fastq_sequences:
            for i in range(len(read)):
                if read[i] == 'c' or read[i] == 'g':
                    gc[i] += 1
                totals[i] += 1
        for i in range(len(gc)):
            if totals[i] > 0:
                gc[i] /= float(totals[i])
        return gc
    def __del__(self):
        if self.handle:
            self.handle.close()

    def how_many_records(self):
        records = len(self.seqs)
        return records

    def sequence_lengths(self):
        lengths = {name: len(seq) for name, seq in self.seqs.items()}
        max_len = max(lengths.values())
        longest_ids = [n for n, L in lengths.items() if L == max_len]
        print(f"Longest length is {max_len}, found in: {longest_ids}")
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
                if len(rfs) != 3:
                    print("Reading frame not valid")
                    return None
        else:
            reading_frames = self.reading_frame()
        open_reading_frames = {}
        longest = {name: 0 for name in reading_frames}
        stop_codons = ['tga', 'tag', 'taa']
        start_codon = 'atg'
        for name, rfs in reading_frames.items():
            open_reading_frames[name] = []
            for i in range(3):
                orfs = {}
                seq = rfs[i]
                orf = ''
                start = None
                in_orf = False
                for j in range(0,len(seq),3):
                    if seq[j:j+3] != start_codon:
                        continue
                    for k in range(j+3, len(seq), 3):
                        codon = seq[k:k + 3]
                        if codon in stop_codons:
                            orf = seq[j:k + 3]
                            orfs[j] = orf
                            break
                open_reading_frames[name].append(orfs)
        return open_reading_frames

    def longest_orf_of_id(self, seq_id : str, open_reading_frames:dict = None):
        longest = {}
        if open_reading_frames is None:
            open_reading_frames = self.orf()
        if seq_id not in open_reading_frames:
            print(f"ID {seq_id!r} not found")
            return None
        best_start, best_orf = None, ""
        for frame_dict in open_reading_frames[seq_id]:
            for start, orf in frame_dict.items():
                if len(orf) > len(best_orf):
                    best_orf = orf
                    best_start = start

        return {seq_id: {"start": best_start, "orf": best_orf}}


    def find_repeats(self, n: int):
        counts = {}
        for name, seq in self.seqs.items():
            L = len(seq)
            for i in range(L - n+1):
                kmer =seq[i:i+n]
                counts[kmer] = counts.get(kmer, 0) + 1
        repeats = {k: v for k, v in counts.items() if v > 1}
        most_frequent = max(repeats.items(), key=lambda kv: kv[1]) if repeats else None
        return repeats, most_frequent

    def count_bases(self, seq_id : str=None):
        counts = {'g':0, 't':0, 'a':0, 'c':0}
        if seq_id is None:
            for name, seq in self.seqs.items():
                for base in seq:
                    counts[base] += 1
            return counts
        elif seq_id not in self.seqs:
            print(f"ID {seq_id!r} not found")
            return None
        else:
            for base in self.seqs[seq_id]:
                counts[base] += 1
            return counts


