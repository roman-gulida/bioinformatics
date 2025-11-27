seq = """ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAA
CGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAAC
TAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTG
TTGCAGCCGATCATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTC
CCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTAC
GTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGG
CTTAGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTTGAACAGCCCTATGTGTTCATCAAACGTTCGGAT
GCTCGAACTGCACCTCATGGTCATGTTATGGTTGAGCTGGTAGCAGAACTCGAAGGCATTCAGTACGGTC
GTAGTGGTGAGACACTTGGTGTCCTTGTCCCTCATGTGGGCGAAATACCAGTGGCTTACCGCAAGGTTCT
TCTTCGTAAGAACGGTAATAAAGGAGCTGGTGGCCATAGTTACGGCGCCGATCTAAAGTCATTTGACTTA
GGCGACGAGCTTGGCACTGATCCTTATGAAGATTTTCAAGAAAACTGGAACACTAAACATAGCAGTGGTG
TTACCCGTGAACTCATGCGTGAGCTTAACGGAGGGGCATACACTCGCTATGTCGATAACAACTTCTGTGG
CCCTGATGGCTACCCTCTTGAGTGCATTAAAGACCTTCTAGCACGTGCTGGTAAAGCTTCATGCACTTTG
TCCGAACAACTGGACTTTATTGACACTAAGAGGGGTGTATACTGCTGCCGTGAACATGAGCATGAAATTG
CTTGGTACACGGAACGTTCTGAAAAGAGCTATGAATTGCAGACACCTTTTGAAATTAAATTGGCAAAGAA
ATTTGACACCTTCAATGGGGAATGTCCAAATTTTGTATTTCCCTTAAATTCCATAATCAAGACTATTCAA
CCAAGGGTTGAAAAGAAAAAGCTTGATGGCTTTATGGGTAGAATTCGATCTGTCTATCCAGTTGCGTCAC
CAAATGAATGCAACCAAATGTGCCTTTCAACTCTCATGAAGTGTGATCATTGTGGTGAAACTTCATGGCA
GACGGGCGATTTTGTTAAAGCCACTTGCGAATTTTGTGGCACTGAGAATTTGACTAAAGAAGGTGCCACT
ACTTGTGGTTACTTACCCCAAAATGCTGTTGTTAAAATTTATTGTCCAGCATGTCACAATTCAGAAGTAG
GACCTGAGCATAGTCTTGCCGAATACCATAATGAATCTGGCTTGAAAACCATTCTTCGTAAGGGTGGTCG
CACTATTGCCTTTGGAGGCTGTGTGTTCTCTTATGTTGGTTGCCATAACAAGTGTGCCTATTGGGTTCCA
CGTGCTAGCGCTAACATAGGTTGTAACCATACAGGTGTTGTTGGAGAAGGTTCCGAAGGTCTTAATGACA
ACCTTCTTGAAATACTCCAAAAAGAGAAAGTCAACATCAATATTGTTGGTGACTTTAAACTTAATGAAGA
GATCGCCATTATTTTGGCATCTTTTTCTGCTTCCACAAGTGCTTTTGTGGAAACTGTGAAAGGTTTGGAT
TATAAAGCATTCAAACAAATTGTTGAATCCTGTGGTAATTTTAAAGTTACAAAAGGAAAAGCTAAAAAAG""".replace('\n', '')

# Restriction enzymes with recognition sequences
enzymes = {
    'EcoRI': 'GAATTC',
    'BamHI': 'GGATCC',
    'HindIII': 'AAGCTT',
    'TaqI': 'TCGA',
    'HaeIII': 'GGCC'
}

# Cleavage positions (after which nucleotide the cut occurs)
cleavage_offsets = {
    'EcoRI': 1,      # G^AATTC
    'BamHI': 1,      # G^GATCC
    'HindIII': 1,    # A^AGCTT
    'TaqI': 1,       # T^CGA
    'HaeIII': 2      # GG^CC
}

def find_restriction_sites(sequence, recognition_seq, offset):
    positions = []
    start = 0
    while True:
        pos = sequence.find(recognition_seq, start)
        if pos == -1:
            break
        cleavage_pos = pos + offset
        positions.append(cleavage_pos)
        start = pos + 1
    return positions

def calculate_fragments(sequence_length, cleavage_positions):
    if not cleavage_positions:
        return [sequence_length]
    
    sorted_positions = sorted(cleavage_positions)
    fragments = []
    
    fragments.append(sorted_positions[0])
    
    for i in range(1, len(sorted_positions)):
        fragments.append(sorted_positions[i] - sorted_positions[i-1])
    
    fragments.append(sequence_length - sorted_positions[-1])
    
    return fragments

print("DNA SEQUENCE LENGTH:", len(seq), "nucleotides\n")

results = {}

for enzyme_name, recognition_seq in enzymes.items():
    offset = cleavage_offsets[enzyme_name]
    cleavage_positions = find_restriction_sites(seq, recognition_seq, offset)
    num_cleavages = len(cleavage_positions)
    fragments = calculate_fragments(len(seq), cleavage_positions)
    
    results[enzyme_name] = {
        'positions': cleavage_positions,
        'fragments': fragments,
        'num_cleavages': num_cleavages
    }
    
    print(f"\n{enzyme_name} (Recognition: {recognition_seq})")
    print(f"Number of cleavages: {num_cleavages}")
    
    if num_cleavages > 0:
        print(f"Cleavage positions: {cleavage_positions}")
        print(f"Fragment lengths: {fragments}")
        print(f"Number of fragments: {len(fragments)}")
    else:
        print("No cleavage sites found")
        print(f"Fragment lengths: [{len(seq)}]")

print("\n" + "="*40)
print("GEL ELECTROPHORESIS SIMULATION")
print("="*40)

max_length = len(seq)
gel_height = 30

for enzyme_name in enzymes.keys():
    fragments = results[enzyme_name]['fragments']
    sorted_fragments = sorted(fragments, reverse=True)
    
    print(f"\n{enzyme_name}:")
    print("-" * 50)
    
    if len(fragments) == 1:
        print("No digestion - single band at origin")
        continue
    
    for i, frag_length in enumerate(sorted_fragments):
        migration_distance = int((1 - frag_length / max_length) * gel_height)
        
        spaces = " " * migration_distance
        band = "█" * 10
        
        print(f"{spaces}{band}  {frag_length} bp")

