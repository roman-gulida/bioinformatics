import matplotlib.pyplot as plt
from collections import Counter

GENETIC_CODE = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

def parse_fasta(filename):
    sequence = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence += line.strip().upper()
    return sequence.replace('T', 'U')

def extract_codons_from_aug(seq):
    codons = []
    for i in range(len(seq) - 2):
        if seq[i:i+3] == 'AUG':
            for j in range(i, len(seq) - 2, 3):
                codon = seq[j:j+3]
                if len(codon) == 3:
                    codons.append(codon)
                    if codon in ['UAA', 'UAG', 'UGA']:
                        break
    return codons

def analyze_codons(codons):
    codon_freq = Counter(codons)
    top10_codons = codon_freq.most_common(10)
    
    aa_list = [GENETIC_CODE.get(codon, 'Unknown') for codon in codons]
    aa_freq = Counter(aa_list)
    top3_aa = aa_freq.most_common(3)
    
    return top10_codons, top3_aa

def plot_codons(top10, title, color):
    codons = [f"{c[0]}\n({GENETIC_CODE.get(c[0], '?')})" for c in top10]
    counts = [c[1] for c in top10]
    
    plt.figure(figsize=(10, 5))
    plt.bar(codons, counts, color=color, alpha=0.7, edgecolor='black')
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xlabel('Codon (Amino Acid)', fontsize=11)
    plt.ylabel('Frequency', fontsize=11)
    plt.xticks(fontsize=9)
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.show()

def main():
    covid_file = 'covid19_genome.fna'
    flu_file = 'influenza_genome.fna'
    
    covid_rna = parse_fasta(covid_file)
    flu_rna = parse_fasta(flu_file)
    
    print("Extracting codons from AUG...")
    covid_codons = extract_codons_from_aug(covid_rna)
    flu_codons = extract_codons_from_aug(flu_rna)
    
    print(f"COVID-19: {len(covid_codons)} codons analyzed")
    print(f"Influenza: {len(flu_codons)} codons analyzed\n")
    
    covid_top10, covid_top3_aa = analyze_codons(covid_codons)
    flu_top10, flu_top3_aa = analyze_codons(flu_codons)
    
    # a) COVID-19 top 10 chart
    print("=" * 50)
    print("a) COVID-19 Top 10 Most Frequent Codons:")
    for i, (codon, count) in enumerate(covid_top10, 1):
        aa = GENETIC_CODE.get(codon, 'Unknown')
        print(f"{i:2}. {codon} ({aa}): {count}")
    plot_codons(covid_top10, "a) COVID-19 Top 10 Most Frequent Codons", '#3B82F6')
    
    # b) Influenza top 10 chart
    print("\n" + "=" * 50)
    print("b) Influenza Top 10 Most Frequent Codons:")
    for i, (codon, count) in enumerate(flu_top10, 1):
        aa = GENETIC_CODE.get(codon, 'Unknown')
        print(f"{i:2}. {codon} ({aa}): {count}")
    plot_codons(flu_top10, "b) Influenza Top 10 Most Frequent Codons", '#9333EA')
    
    # c) Compare common codons
    print("\n" + "=" * 50)
    print("c) Common Codons in Top 10 of Both Genomes:")
    covid_set = {c[0] for c in covid_top10}
    flu_set = {c[0] for c in flu_top10}
    common = covid_set & flu_set
    
    if common:
        for codon in sorted(common):
            covid_count = next(c[1] for c in covid_top10 if c[0] == codon)
            flu_count = next(c[1] for c in flu_top10 if c[0] == codon)
            aa = GENETIC_CODE.get(codon, 'Unknown')
            print(f"  {codon} ({aa}) - COVID: {covid_count}, Flu: {flu_count}")
    else:
        print("  No common codons found")
    
    # d) Top 3 amino acids
    print("\n" + "=" * 50)
    print("d) Top 3 Amino Acids for Each Genome:")
    print("\nCOVID-19:")
    for i, (aa, count) in enumerate(covid_top3_aa, 1):
        print(f"  {i}. {aa}: {count} occurrences")
    
    print("\nInfluenza:")
    for i, (aa, count) in enumerate(flu_top3_aa, 1):
        print(f"  {i}. {aa}: {count} occurrences")
    print("=" * 50)

if __name__ == "__main__":
    main()