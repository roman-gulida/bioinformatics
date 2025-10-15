import matplotlib.pyplot as plt
import math

def read_fasta(filename):
    sequence = ""
    reading = False 
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if reading:
                    break
                reading = True
                continue
            if reading:
                sequence += line
    return sequence.upper()


def calculate_tm(seq, method = 2, Na = 0.0001):
    seq = seq.upper()
  
    A = seq.count("A")
    T = seq.count("T")
    G = seq.count("G")
    C = seq.count("C")
    length = len(seq)

    if method == 1:
        tm = 4 * (G + C) + 2 * (A + T)
    elif method == 2:
        gc_percent = ((G + C) / length) * 100
        tm = -(81.5 + 16.6 * math.log10(Na) + 0.41 * (gc_percent / 100) - (600 / length))
    else:
        raise ValueError("Method must be 1 or 2.")
    return tm


def sliding_window_tm(sequence, window_size = 8, Na = 0.001):
    results = []
    for i in range(len(sequence) - window_size + 1):
        window_seq = sequence[i:i+window_size]
        tm_values = calculate_tm(window_seq)
        results.append({
            "start": i + 1,
            "end": i + window_size,
            "window": window_seq,
            "tm": tm_values
        })
    return results


def plot_tm(results):
    positions = [r["start"] for r in results]
    tm_values = [r["tm"] for r in results]

    plt.figure(figsize=(10, 5))
    plt.plot(positions, tm_values, marker='.', linestyle='-', linewidth=1.5)
    plt.title("Sliding Window Tm (window size = 8)")
    plt.xlabel("Window start position along sequence")
    plt.ylabel("Melting Temperature (C)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def main():
    fasta_file = "mitochondrion.1.1.genomic.fna"
    try:
        seq = read_fasta(fasta_file)
    except FileNotFoundError:
        print("File not found")
        return

    results = sliding_window_tm(seq)
    
    print("\nSliding Window tm Results:")
    for r in results:
        print(f"{r['start']:>4}-{r['end']:>4} | {r['window']} | tm={r['tm']:.5f} C")

    plot_tm(results)

if __name__ == "__main__":
    main()
