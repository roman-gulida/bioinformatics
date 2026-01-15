import numpy as np
import json
import re


def compute_word_transition_matrix(text):
    # Extract words
    words = re.findall(r"\b\w+\b", text.lower())

    # Get unique words and create mapping to ASCII symbols
    unique_words = []
    seen = set()
    for word in words:
        if word not in seen:
            unique_words.append(word)
            seen.add(word)

    # Map each word to an ASCII symbol (starting from '!')
    word_to_symbol = {word: chr(33 + i) for i, word in enumerate(unique_words)}
    word_to_idx = {word: i for i, word in enumerate(unique_words)}

    n = len(unique_words)
    counts = np.zeros((n, n))

    # Count transitions between consecutive words
    for i in range(len(words) - 1):
        from_word = words[i]
        to_word = words[i + 1]
        counts[word_to_idx[from_word], word_to_idx[to_word]] += 1

    # Convert counts to probabilities
    transition_matrix = np.zeros((n, n))
    for i in range(n):
        row_sum = counts[i].sum()
        if row_sum > 0:
            transition_matrix[i] = counts[i] / row_sum

    return transition_matrix, unique_words, word_to_symbol


def save_matrix_to_json(
    matrix, words, word_to_symbol, filename="word_transition_matrix_t3.json"
):
    matrix_dict = {
        "words": words,
        "word_to_symbol": word_to_symbol,
        "matrix": matrix.tolist(),
    }

    with open(filename, "w") as f:
        json.dump(matrix_dict, f, indent=2)

    print(f"Transition matrix saved to {filename}")


def main():
    text = """The quick brown fox jumps over the lazy dog. The dog was sleeping under a tree. 
    A bird flew over the tree and landed on a branch. The fox saw the bird and tried to catch it. 
    The bird quickly flew away. The dog woke up and barked at the fox. The fox ran into the forest. 
    The forest was dark and quiet. The moon shone brightly in the sky."""

    print(f"Text: {text}")
    print(f"Length: {len(text)} characters\n")

    transition_matrix, words, word_to_symbol = compute_word_transition_matrix(text)

    print(f"Number of unique words: {len(words)}")
    print("\nWord to Symbol Mapping:")
    for word, symbol in word_to_symbol.items():
        print(f"{word}: {symbol}")

    save_matrix_to_json(transition_matrix, words, word_to_symbol)


if __name__ == "__main__":
    main()
