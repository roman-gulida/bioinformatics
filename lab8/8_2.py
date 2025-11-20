import sys
from pathlib import Path

def reverse_complement(seq):
    """Return reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def read_fasta(filename):
    """Read FASTA file and return sequence as single string"""
    seq = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                seq.append(line.upper())
    return ''.join(seq)

def find_inverted_repeats(sequence, min_size=4, max_size=6):
    """Find all inverted repeats in sequence"""
    repeats = []
    seq_len = len(sequence)
    
    for size in range(min_size, max_size + 1):
        for i in range(seq_len - size + 1):
            left_ir = sequence[i:i + size]
            left_rc = reverse_complement(left_ir)
            
            # Search for matching inverted repeat downstream
            for j in range(i + size, seq_len - size + 1):
                right_ir = sequence[j:j + size]
                
                if left_rc == right_ir:
                    repeats.append({
                        'left_start': i,
                        'left_end': i + size - 1,
                        'right_start': j,
                        'right_end': j + size - 1,
                        'ir_size': size,
                        'spacer_size': j - (i + size),
                        'total_length': (j + size) - i
                    })
    
    return repeats

def merge_overlapping_tes(repeats):
    """Merge overlapping or nested transposable elements"""
    if not repeats:
        return []
    
    # Sort by left_start, then by total_length (descending)
    sorted_repeats = sorted(repeats, key=lambda x: (x['left_start'], -x['total_length']))
    
    merged = []
    current = sorted_repeats[0]
    
    for next_te in sorted_repeats[1:]:
        # Check if next TE overlaps or is nested within current
        if next_te['left_start'] <= current['right_end']:
            # Keep the longer one or extend if they overlap
            if next_te['right_end'] > current['right_end']:
                # Extend current if next extends beyond
                if next_te['total_length'] > current['total_length']:
                    current = next_te
        else:
            # No overlap, save current and move to next
            merged.append(current)
            current = next_te
    
    merged.append(current)
    return merged

def detect_transposable_elements(fasta_files):
    """Detect transposable elements in multiple genome files"""
    
    for fasta_file in fasta_files:
        print(f"Processing: {fasta_file}")
        
        sequence = read_fasta(fasta_file)
        print(f"Genome length: {len(sequence)} bp")
        
        print("\nSearching for inverted repeats (4-6 bp)...")
        repeats = find_inverted_repeats(sequence, min_size=4, max_size=6)
        
        print(f"Found {len(repeats)} potential transposable elements")
        
        print("\nMerging overlapping elements...")
        merged = merge_overlapping_tes(repeats)
        
        print(f"\nFinal count: {merged.__len__()} transposable elements\n")
        print(f"{'Position':<20} {'Length (bp)':<15} {'IR Size':<10} {'Spacer'}")
        print("-" * 80)
        
        for te in merged:
            position = f"{te['left_start']}-{te['right_end']}"
            print(f"{position:<20} {te['total_length']:<15} {te['ir_size']:<10} {te['spacer_size']}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python te_detector.py <fasta_file1> <fasta_file2> <fasta_file3>")
        print("Example: python te_detector.py genome1.fasta genome2.fasta genome3.fasta")
        sys.exit(1)
    
    fasta_files = sys.argv[1:]
    
    # Verify files exist
    for f in fasta_files:
        if not Path(f).exists():
            print(f"Error: File '{f}' not found")
            sys.exit(1)
    
    detect_transposable_elements(fasta_files)