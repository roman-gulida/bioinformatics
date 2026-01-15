import numpy as np
import json
import random


def load_transition_matrix(filename):
    with open(filename, "r") as f:
        data = json.load(f)
    return data


def synthesize_sequence(states, matrix, length, initial_state=None):
    if initial_state is None:
        current_idx = random.randint(0, len(states) - 1)
    else:
        current_idx = states.index(initial_state)

    sequence = [states[current_idx]]

    for _ in range(length - 1):
        probabilities = matrix[current_idx]
        current_idx = np.random.choice(len(states), p=probabilities)
        sequence.append(states[current_idx])

    return sequence


def main():
    # DNA sequence synthesis
    print("DNA SEQUENCE SYNTHESIS")

    dna_data = load_transition_matrix("lab13/transition_matrix_t2.json")
    dna_states = dna_data["nucleotides"]
    dna_matrix = np.array(dna_data["matrix"])

    dna_seq = synthesize_sequence(dna_states, dna_matrix, 50)
    print("Synthesized DNA sequence (50 nucleotides):")
    print("".join(dna_seq))

    print("Original sequence:")
    print("TCTTCGTAAGAACGGTAATAAAGGAGCTGGTGGCCATAGTTACGGCGCCGAT")

    # Word sequence synthesis
    print("\nWORD SEQUENCE SYNTHESIS")

    word_data = load_transition_matrix("lab13/word_transition_matrix_t3.json")
    word_states = word_data["words"]
    word_matrix = np.array(word_data["matrix"])

    word_seq = synthesize_sequence(word_states, word_matrix, 50, initial_state="the")
    print("Synthesized text (50 words):")
    print(" ".join(word_seq))

    print("Original text:")
    print("""The quick brown fox jumps over the lazy dog. The dog was sleeping under a tree. 
    A bird flew over the tree and landed on a branch. The fox saw the bird and tried to catch it. 
    The bird quickly flew away. The dog woke up and barked at the fox. The fox ran into the forest. 
    The forest was dark and quiet. The moon shone brightly in the sky.""")


if __name__ == "__main__":
    main()
