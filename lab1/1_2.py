def frequencies(seq):
    total = len(seq)
    alphabet = set(seq)
    
    frequencies = {}
    for symbol in alphabet:
        count = seq.count(symbol)
        frequencies[symbol] = count / total   
    
    return frequencies

S = "ATTTCGCCGATA"
print(frequencies(S))
