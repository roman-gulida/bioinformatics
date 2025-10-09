S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"
nucleotides = ['A', 'T', 'C', 'G']

dinucleotides = []
for n1 in nucleotides:
    for n2 in nucleotides:
        dinucleotides.append(n1 + n2)

trinucleotides = []
for n1 in nucleotides:
    for n2 in nucleotides:
        for n3 in nucleotides:
            trinucleotides.append(n1 + n2 + n3)

def calculate_percentage(seq, combinations):
    total_positions = len(seq) - len(combinations[0]) + 1
    return {
        combo: (seq.count(combo) / total_positions) * 100
        for combo in combinations
    }

dinucleotide_percentages =calculate_percentage(S, dinucleotides)
trinucleotide_percentages = calculate_percentage(S, trinucleotides)

print("Dinucleotide Percentages:")
for k, v in dinucleotide_percentages.items():
    print(f"{k}: {v:.2f}%")

print("\nTrinucleotide Percentages:")
for k, v in trinucleotide_percentages.items():
    print(f"{k}: {v:.2f}%")
