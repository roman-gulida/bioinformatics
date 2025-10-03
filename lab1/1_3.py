import tkinter as tk
from tkinter import filedialog, ttk, messagebox

def parse_fasta(filepath):
    """Read FASTA file and return the sequence as a generator (memory-efficient)."""
    with open(filepath, "r") as f:
        next(f)  # skip header
        for line in f:
            line = line.strip()
            if line.startswith(">"):  # start of a new sequence → stop
                break
            yield line

def calculate_frequencies(sequence_generator):
    """Calculate relative frequencies of symbols in a sequence."""
    counts = {}
    total = 0
    for seq_part in sequence_generator:
        for char in seq_part:
            counts[char] = counts.get(char, 0) + 1
            total += 1
    
    frequencies = {char: (count / total) * 100 for char, count in counts.items()}
    return frequencies

def open_file():
    filepath = filedialog.askopenfilename(
        filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
    )
    if not filepath:
        return

    try:
        sequence_generator = parse_fasta(filepath)
        freqs = calculate_frequencies(sequence_generator)

        # Clear old table
        for row in tree.get_children():
            tree.delete(row)

        # Insert new data
        for symbol, percent in sorted(freqs.items()):
            tree.insert("", "end", values=(symbol, f"{percent:.4f}%"))

        messagebox.showinfo("Success", "Frequencies calculated successfully!")

    except Exception as e:
        messagebox.showerror("Error", str(e))

# GUI Setup
root = tk.Tk()
root.title("FASTA Sequence Analyzer")
root.geometry("400x300")

btn = tk.Button(root, text="Open FASTA File", command=open_file)
btn.pack(pady=10)

tree = ttk.Treeview(root, columns=("Symbol", "Percentage"), show="headings")
tree.heading("Symbol", text="Symbol")
tree.heading("Percentage", text="Percentage")
tree.pack(expand=True, fill="both")

root.mainloop()
