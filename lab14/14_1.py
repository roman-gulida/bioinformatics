import math

# Known sequences
S1 = "ATCGATTCGATATCATACACGTAT"  # CpG island (+)
S2 = "CTCGACTAGTATGAAGTCCACGCTTG"  # Non-CpG island (-)
test_sequence = "CAGGTTGGAAACGTAA"

bases = ["A", "C", "G", "T"]


def count_transitions(sequence):
    """Count transition frequencies from a sequence"""
    # Initialize count matrix
    counts = {base: {b: 0 for b in bases} for base in bases}

    # Count transitions
    for i in range(len(sequence) - 1):
        from_base = sequence[i]
        to_base = sequence[i + 1]
        if from_base in bases and to_base in bases:
            counts[from_base][to_base] += 1

    return counts


def calculate_frequencies(counts):
    """Convert counts to frequencies"""
    frequencies = {base: {b: 0.0 for b in bases} for base in bases}

    for from_base in bases:
        total = sum(counts[from_base].values())
        if total > 0:
            for to_base in bases:
                frequencies[from_base][to_base] = counts[from_base][to_base] / total

    return frequencies


def print_matrix(matrix, title):
    """Print a matrix in a readable format"""
    print(f"\n{title}")
    print(f"{'':>5}", end="")
    for base in bases:
        print(f"{base:>8}", end="")
    print()

    for from_base in bases:
        print(f"{from_base:>5}", end="")
        for to_base in bases:
            print(f"{matrix[from_base][to_base]:>8.3f}", end="")
        print()


def calculate_log_likelihood_matrix(freq_plus, freq_minus):
    """Calculate log2 likelihood ratio matrix"""
    log_likelihood = {base: {b: 0.0 for b in bases} for base in bases}

    for from_base in bases:
        for to_base in bases:
            if freq_plus[from_base][to_base] > 0 and freq_minus[from_base][to_base] > 0:
                ratio = freq_plus[from_base][to_base] / freq_minus[from_base][to_base]
                # Using change of base formula: log2(x) = log(x) / log(2)
                log_likelihood[from_base][to_base] = math.log(ratio) / math.log(2)
            elif (
                freq_plus[from_base][to_base] > 0
                and freq_minus[from_base][to_base] == 0
            ):
                # Handle division by zero - assign a large positive value
                log_likelihood[from_base][to_base] = 10.0
            elif (
                freq_plus[from_base][to_base] == 0
                and freq_minus[from_base][to_base] > 0
            ):
                # Handle zero numerator - assign a large negative value
                log_likelihood[from_base][to_base] = -10.0
            else:
                # Both are zero
                log_likelihood[from_base][to_base] = 0.0

    return log_likelihood


def score_sequence(sequence, log_likelihood):
    """Calculate the log-likelihood score for a sequence"""
    score = 0.0
    transitions = []

    for i in range(len(sequence) - 1):
        from_base = sequence[i]
        to_base = sequence[i + 1]

        if from_base in bases and to_base in bases:
            transition_score = log_likelihood[from_base][to_base]
            score += transition_score
            transitions.append((from_base, to_base, transition_score))

    return score, transitions


# Step 1: Count transitions for CpG+ model (S1)
print("STEP 1: Count transitions for CpG+ model (S1)")
print(f"S1 sequence: {S1}")
counts_plus = count_transitions(S1)
print_matrix(counts_plus, "Transition Counts (+)")

# Step 2: Count transitions for CpG- model (S2)
print("\nSTEP 2: Count transitions for CpG- model (S2)")
print(f"S2 sequence: {S2}")
counts_minus = count_transitions(S2)
print_matrix(counts_minus, "Transition Counts (-)")

# Calculate frequencies
freq_plus = calculate_frequencies(counts_plus)
freq_minus = calculate_frequencies(counts_minus)

print_matrix(freq_plus, "Transition Frequencies (+)")
print_matrix(freq_minus, "Transition Frequencies (-)")

# Step 3: Calculate log-likelihood matrix
print("\nSTEP 3: Calculate Log-Likelihood Matrix")
log_likelihood = calculate_log_likelihood_matrix(freq_plus, freq_minus)
print_matrix(log_likelihood, "Log-Likelihood Matrix β = log2(Tr+/Tr-)")

# Step 4: Score the test sequence
print("\nSTEP 4: Score Test Sequence")
print(f"Test sequence: {test_sequence}")

score, transitions = score_sequence(test_sequence, log_likelihood)

print("\nDetailed scoring:")
print(f"{'Position':>10} {'Transition':>12} {'Score':>10}")
print("-" * 35)
for i, (from_base, to_base, trans_score) in enumerate(transitions):
    print(f"{i:>10} {from_base}->{to_base:>10} {trans_score:>10.3f}")

print(f"\n{'Total Log-Likelihood Score:':<30} {score:.3f}")
print(f"{'RESULT:':<20}", end="")

if score > 0:
    print("CpG ISLAND (+)")
    print("The sequence is MORE LIKELY to belong to a CpG island.")
    print("Positive score indicates higher similarity to the (+) model.")
elif score < 0:
    print("NON-CpG ISLAND (-)")
    print("The sequence is MORE LIKELY to belong to non-CpG regions.")
    print("Negative score indicates higher similarity to the (-) model.")
else:
    print("UNCERTAIN")
    print("The sequence shows no preference for either model.")
