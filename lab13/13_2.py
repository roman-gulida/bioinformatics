import numpy as np
import json


def compute_transition_matrix(sequence):
    nucleotides = ["A", "C", "G", "T"]
    n = len(nucleotides)

    # Count transitions
    counts = np.zeros((n, n))
    nuc_to_idx = {nuc: i for i, nuc in enumerate(nucleotides)}

    for i in range(len(sequence) - 1):
        from_nuc = sequence[i]
        to_nuc = sequence[i + 1]
        if from_nuc in nuc_to_idx and to_nuc in nuc_to_idx:
            counts[nuc_to_idx[from_nuc], nuc_to_idx[to_nuc]] += 1

    # Convert counts to probabilities
    transition_matrix = np.zeros((n, n))
    for i in range(n):
        row_sum = counts[i].sum()
        if row_sum > 0:
            transition_matrix[i] = counts[i] / row_sum

    return transition_matrix, nucleotides


def save_matrix_to_json(matrix, nucleotides, filename="transition_matrix_t2.json"):
    matrix_dict = {"nucleotides": nucleotides, "matrix": matrix.tolist()}

    with open(filename, "w") as f:
        json.dump(matrix_dict, f, indent=2)

    print(f"Transition matrix saved to {filename}")


def main():
    dna_sequence = "TCTTCGTAAGAACGGTAATAAAGGAGCTGGTGGCCATAGTTACGGCGCCGAT"

    print(f"DNA Sequence: {dna_sequence}")
    print(f"Length: {len(dna_sequence)}\n")

    transition_matrix, nucleotides = compute_transition_matrix(dna_sequence)

    print("Transition Matrix:")
    print("     ", "  ".join(nucleotides))
    for i, nuc in enumerate(nucleotides):
        print(
            f"{nuc}:   ",
            "  ".join(
                f"{transition_matrix[i, j]:.2f}" for j in range(len(nucleotides))
            ),
        )

    save_matrix_to_json(transition_matrix, nucleotides)


if __name__ == "__main__":
    main()
