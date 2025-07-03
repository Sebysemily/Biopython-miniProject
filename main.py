import matplotlib

from utils.FastaAnalyzer import FastaAnalyzer
import matplotlib.pyplot as plt

reference_fa = FastaAnalyzer("data/phix.fa")
reads_fa = FastaAnalyzer("data/")

for read_name, read_seq in reads_fa.get_sequences().items():
    matches = reference_fa.naive(read_seq)


'''
gc = fa.fastq_find_gc_by_pos()
plt.plot(range(len(gc)), gc)
plt.savefig("figures/gc%_of_SRR835775_1.png",
            dpi=300,
            bbox_inches='tight')
plt.show()

h = fa.fastq_create_hist()
plt.bar(range(len(h)), h)
plt.savefig("figures/hqual_of_SRR835775_1.png",
            dpi=300,
            bbox_inches='tight')
plt.show()


base_counts = fa.count_bases()
print(base_counts)
length = fa.sequence_lengths()
print(length)

records = fa.how_many_records()
print(f"there is {records}")

lengths = fa.sequence_lengths()
min_len = min(lengths.values())
shortest_ids = [n for n, L in lengths.items() if L == min_len]
print(f"Shortest length is {min_len}, found in: {shortest_ids}")
rfs = fa.reading_frame()
orfs = fa.orf(rfs)
for seq_id in orfs:
    summary = fa.longest_orf_of_id(seq_id, orfs)
    start = summary[seq_id]["start"]
    orf_seq = summary[seq_id]["orf"]
    print(f"   Seq {seq_id}: longest ORF starts at {start}, length {len(orf_seq)}")

repeats_dict, top = fa.find_repeats(12)
if top:
    max_count = top[1]
else:
    max_count = 0

num_max = sum(1 for cnt in repeats_dict.values() if cnt == max_count)
print(f"{num_max} different 12-base sequences occur {max_count} times")
frame2_orfs = {seq_id: frames[2]
               for seq_id, frames in orfs.items()}


for seq_id, orf_dict in frame2_orfs.items():
    print(f"\n— seq {seq_id}, frame 2 —")
    if not orf_dict:
        print("  (no ORFs found)")
    for start, orf_seq in orf_dict.items():
        print(f"  at {start}: length={len(orf_seq)}  seq={orf_seq}")

max_len = 0
winner = None   # will hold (seq_id, start, orf_seq)
for seq_id, orf_dict in frame2_orfs.items():
    for start, orf_seq in orf_dict.items():
        if len(orf_seq) > max_len:
            max_len = len(orf_seq)
            winner = (seq_id, start, orf_seq)
'''