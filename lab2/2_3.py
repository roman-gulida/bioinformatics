import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt

# --- Function to read FASTA file ---
def read_fasta(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
    # Skip header lines starting with '>'
    sequence = "".join([line.strip() for line in lines if not line.startswith(">")])
    return sequence.upper()  # ensure uppercase

# --- Function to calculate relative frequencies with sliding window ---
def sliding_window_frequencies(sequence, window_size=30):
    alphabet = sorted(list(set(sequence)))  # all unique symbols
    frequencies = {symbol: [] for symbol in alphabet}
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]
        for symbol in alphabet:
            freq = window.count(symbol) / window_size
            frequencies[symbol].append(freq)
    return frequencies

# --- Function to plot results ---
def plot_frequencies(frequencies):
    plt.figure(figsize=(12, 6))
    for symbol, values in frequencies.items():
        plt.plot(values, label=symbol)
    plt.xlabel("Window Position")
    plt.ylabel("Relative Frequency")
    plt.title("Sliding Window Relative Frequencies")
    plt.legend()
    plt.show()

# --- GUI function to select file and run analysis ---
def select_file():
    file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fna")])
    if not file_path:
        return
    try:
        sequence = read_fasta(file_path)
        if len(sequence) < 30:
            messagebox.showerror("Error", "Sequence is too short for a sliding window of 30.")
            return
        frequencies = sliding_window_frequencies(sequence, 30)
        plot_frequencies(frequencies)
    except Exception as e:
        messagebox.showerror("Error", str(e))

# --- GUI Setup ---
root = tk.Tk()
root.title("DNA Sliding Window Analyzer")
root.geometry("400x150")

label = tk.Label(root, text="Select a FASTA file to analyze:", font=("Arial", 12))
label.pack(pady=20)

btn = tk.Button(root, text="Browse FASTA File", command=select_file)
btn.pack(pady=10)

root.mainloop()
