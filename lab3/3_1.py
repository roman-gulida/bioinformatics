import math

def calculate_tm(seq, method = 1, Na = 0.001):
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


if __name__ == "__main__":
    seq = "ATTGCCC"
    tm_value1 = calculate_tm(seq, method=1)
    tm_value2 = calculate_tm(seq, method=2)
    
    print(f"\n tm with method 1: {tm_value1:.2f} C")
    print(f"\n tm with method 2: {tm_value2:.2f} C")

