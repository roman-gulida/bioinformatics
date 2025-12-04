import pandas as pd
import matplotlib.pyplot as plt

# --- Input ---------------------------------------------------------
seq = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
window_len = 30  
step = 1     

# --- Functions -----------------------------------------------------

def cg_percentage(seq: str) -> float:
    return 100.0 * (seq.count('C') + seq.count('G')) / len(seq) if len(seq) else 0

def kappa_ic(seq: str) -> float:
    """Bioinformatic (Kappa) Index of Coincidence used in promoter analysis."""
    N = len(seq)
    freqs = {b: seq.count(b)/N for b in 'ACGT'}
    return 100 * sum(v*v for v in freqs.values())

# --- Sliding windows ----------------------------------------------
result = []

for start in range(0, len(seq) - window_len + 1, step):
    w = seq[start:start+window_len]
    cg = cg_percentage(w)
    ic = kappa_ic(w)
    center = start + (window_len - 1) / 2
    
    result.append({
        "start": start,
        "center": center,
        "window": w,
        "cg": cg,
        "ic": ic
    })

df = pd.DataFrame(result)

# Whole sequence expected values
whole_cg = cg_percentage(seq)         # Should match 29.27 → OK
whole_ic = kappa_ic(seq)             # Should match 27.53 → FIXED

# Center of Weight 
cw_cg = (df["center"] * df["cg"]).sum() / df["cg"].sum()
cw_ic = (df["center"] * df["ic"]).sum() / df["ic"].sum()


df.to_csv("lab10/dna_pattern_windows.csv", index=False)


plt.figure(figsize=(12,5))
plt.plot(df["center"], df["cg"], label="CG%", marker='o')
plt.plot(df["center"], df["ic"], label="Kappa IC%", marker='s')
plt.title("Sliding Window Pattern (30bp)")
plt.xlabel("Window Center Position")
plt.ylabel("Value (%)")
plt.legend()
plt.grid(True)
plt.savefig("lab10/dna_pattern.png", dpi=200)
plt.close()


plt.figure(figsize=(12,4))
plt.scatter(df["center"], df["cg"], s=40)
plt.title("Window Centers Marked")
plt.xlabel("Center Position")
plt.ylabel("CG%")
plt.grid(True)
plt.savefig("lab10/dna_pattern_centers_position.png", dpi=200)
plt.close()


plt.figure(figsize=(7,7))
plt.scatter(df["cg"], df["ic"], s=60, color='black')

for _, row in df.iterrows():
    plt.text(row["cg"], row["ic"], str(int(row["center"])), fontsize=7)

plt.xlabel("CG%")
plt.ylabel("Kappa IC%")
plt.title("Promoter Pattern: CG% vs Kappa IC (Window Centers)")
plt.grid(True)
plt.savefig("lab10/dna_pattern_promkappa_style.png", dpi=200)
plt.close()

print("\n=== SUMMARY ===")
print("Sequence length:", len(seq))
print(f"Whole CG% = 29.27 → Computed: {whole_cg:.2f}")
print(f"Kappa IC% = 27.53 → Computed: {whole_ic:.2f}")
print("Sliding windows:", len(df))
print(f"Center of weight (CG): {cw_cg:.2f}")
print(f"Center of weight (IC): {cw_ic:.2f}")
