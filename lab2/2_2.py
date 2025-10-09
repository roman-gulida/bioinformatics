S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"
S2= 'ABAA'

def count_kmers(seq, nucleotides):
    counts = {}
    for i in range(len(seq) - nucleotides + 1):
        nucleotide_fragment = seq[i:i+nucleotides]
        counts[nucleotide_fragment] = counts.get(nucleotide_fragment, 0) + 1
    total_positions = len(seq) - nucleotides + 1
    percentages = {nucleotide_fragment: (count / total_positions) * 100 for nucleotide_fragment, count in counts.items()}
    return percentages

dinucleotide_percentages = count_kmers(S2, 2)
trinucleotide_percentages = count_kmers(S2, 3)

print("Dinucleotide Percentages:")
for nucleotide_fragment, pct in dinucleotide_percentages.items():
    print(f"{nucleotide_fragment}: {pct:.2f}%")

print("\nTrinucleotide Percentages:")
for nucleotide_fragment, pct in trinucleotide_percentages.items():
    print(f"{nucleotide_fragment}: {pct:.2f}%")
