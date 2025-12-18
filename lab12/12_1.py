import math

motifs = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTGAGT",
    "AAGGTAAGT",
]

S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

bases = ["A", "C", "G", "T"]
motif_len = len(motifs[0])
num_seqs = len(motifs)

# COUNT MATRIX
count_matrix = {b: [0] * motif_len for b in bases}

for seq in motifs:
    for i, nt in enumerate(seq):
        count_matrix[nt][i] += 1

print("COUNT MATRIX")
for b in bases:
    print(b, count_matrix[b])
print()

# WEIGHT & RELATIVE FREQUENCY MATRIX
freq_matrix = {
    b: [count_matrix[b][i] / num_seqs for i in range(motif_len)] for b in bases
}

print("RELATIVE FREQUENCY MATRIX")
for b in bases:
    print(b, freq_matrix[b])
print()

# LOG-LIKELIHOOD MATRIX
null_prob = 0.25
ll_matrix = {b: [] for b in bases}

for b in bases:
    for p in freq_matrix[b]:
        if p == 0:
            ll_matrix[b].append(float("-inf"))
        else:
            ll_matrix[b].append(math.log(p / null_prob))

print("LOG-LIKELIHOOD MATRIX")
for b in bases:
    print(b, ll_matrix[b])
print()


# SLIDING WINDOW SCORING
def score_window(window):
    score = 0
    for i, nt in enumerate(window):
        score += ll_matrix[nt][i]
    return score


scores = []

print("SLIDING WINDOW SCORES")
for i in range(len(S) - motif_len + 1):
    window = S[i : i + motif_len]
    score = score_window(window)
    scores.append((i, window, score))
    print(f"Pos {i:2d}  {window}  Score = {score:.3f}")

# SIGNAL DETECTION
best = max(scores, key=lambda x: x[2])

print("\nBEST MATCH")
print("Position:", best[0])
print("Sequence:", best[1])
print("Score:", best[2])

if best[2] > 0:
    print("\nConclusion: Exon–intron border signal detected.")
else:
    print("\nConclusion: No significant exon–intron border detected.")


print(
    "\nMost sliding windows obtain a score of −∞ because the log-likelihood model is built directly from observed frequencies without pseudocounts. "
    "Any nucleotide that never appeared at a given position in the training motifs has probability zero, causing the log-likelihood to be −∞."
    " This behavior is expected and highlights highly conserved positions in the exon–intron boundary motif."
)
