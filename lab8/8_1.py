import random
import re


def random_dna(length):
    return "".join(random.choice("ACGT") for _ in range(length))

host_dna = random_dna(random.randint(200, 400))

transposons = [
    ("TIRL1", "ATGCGA", "TCGCAT"),   
    ("TIRL2", "GGATCC", "GGATCC"),
    ("TIRL3", "TTAACC", "GGT TAA".replace(" ", "")),
    ("TIRL4", "CCCGGG", "CCCGGG")
]

cores = [
    "ACGTACGTACGT",
    "TTGACCTTGACC",
    "GATCGATCGATC",
    "TAACTAACTAAC"
]


dna_list = list(host_dna)

insert_positions = sorted(random.sample(range(50, len(host_dna)-50), 4))
transposon_locations = []

for i, (name, IRL, IRR) in enumerate(transposons):
    pos = insert_positions[i]
    seq = IRL + cores[i] + IRR

    # Overwrite into DNA (allows overlap)
    dna_list[pos:pos+len(seq)] = seq
    transposon_locations.append((name, pos, pos+len(seq)))

host_dna = "".join(dna_list)

print("Generated DNA: ", host_dna)
print("Generated DNA length:", len(host_dna))
print("Actual inserted transposons:")
for t in transposon_locations:
    print(" ", t)


def detect_transposons(dna, transposons):
    detections = []

    for (name, IRL, IRR) in transposons:
        for match in re.finditer(IRL, dna):
            start_IRL = match.start()

            search_region = dna[start_IRL + len(IRL): start_IRL + 100]
            irr_match = search_region.find(IRR)

            if irr_match != -1:
                start = start_IRL
                end = start_IRL + len(IRL) + irr_match + len(IRR)
                detections.append((name, start, end))

    return detections


print("\nDetected transposons:")
for d in detect_transposons(host_dna, transposons):
    print(" ", d)
