import matplotlib.pyplot as plt


def smith_waterman(seq1, seq2, match=3, mismatch=-3, gap=-2):
    """
    Smith-Waterman local alignment algorithm.
    Returns aligned sequences and the alignment score.
    """
    n, m = len(seq1), len(seq2)

    # Initialize scoring matrix
    matrix = [[0] * (m + 1) for _ in range(n + 1)]
    max_score = 0
    max_pos = (0, 0)

    # Fill the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Calculate scores for each direction
            match_score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diagonal = matrix[i - 1][j - 1] + match_score
            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap

            # Take maximum (or 0 for local alignment)
            matrix[i][j] = max(0, diagonal, up, left)

            # Track maximum score position
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)

    # Traceback from maximum score
    aligned1, aligned2 = "", ""
    i, j = max_pos

    while i > 0 and j > 0 and matrix[i][j] > 0:
        current = matrix[i][j]
        diagonal = matrix[i - 1][j - 1]
        up = matrix[i - 1][j]
        left = matrix[i][j - 1]

        match_score = match if seq1[i - 1] == seq2[j - 1] else mismatch

        if current == diagonal + match_score:
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = seq2[j - 1] + aligned2
            i -= 1
            j -= 1
        elif current == up + gap:
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1
        elif current == left + gap:
            aligned1 = "-" + aligned1
            aligned2 = seq2[j - 1] + aligned2
            j -= 1
        else:
            break

    return aligned1, aligned2, max_score


def scoring_equation_1(aligned1, aligned2):
    """
    Scoring Equation 1: Percent Identity
    Calculates the percentage of exact matches (excluding gaps) over alignment length.
    Formula: (matches / alignment_length) * 100
    """
    if len(aligned1) == 0:
        return 0.0
    matches = sum(1 for c1, c2 in zip(aligned1, aligned2) if c1 == c2 and c1 != "-")
    return (matches / len(aligned1)) * 100


def scoring_equation_2(aligned1, aligned2):
    """
    Scoring Equation 2: Similarity Index with Gap Penalty
    Accounts for matches, mismatches, and gaps with weighted scoring.
    Formula: (3*matches - 1*mismatches - 2*gaps) / alignment_length
    """
    if len(aligned1) == 0:
        return 0.0

    matches = 0
    mismatches = 0
    gaps = 0

    for c1, c2 in zip(aligned1, aligned2):
        if c1 == "-" or c2 == "-":
            gaps += 1
        elif c1 == c2:
            matches += 1
        else:
            mismatches += 1

    score = (3 * matches - 1 * mismatches - 2 * gaps) / len(aligned1)
    return score


def scoring_equation_3(aligned1, aligned2):
    """
    Scoring Equation 3: Normalized Alignment Score
    Combines match quality and alignment coverage.
    Formula: (matches - 0.5*mismatches - 0.3*gaps) / max(len(seq1_ungapped), len(seq2_ungapped))
    """
    matches = 0
    mismatches = 0
    gaps = 0

    for c1, c2 in zip(aligned1, aligned2):
        if c1 == "-" or c2 == "-":
            gaps += 1
        elif c1 == c2:
            matches += 1
        else:
            mismatches += 1

    seq1_len = len(aligned1) - aligned1.count("-")
    seq2_len = len(aligned2) - aligned2.count("-")
    max_len = max(seq1_len, seq2_len)

    if max_len == 0:
        return 0.0

    score = (matches - 0.5 * mismatches - 0.3 * gaps) / max_len
    return score


def chunked_alignment(seq1, seq2, chunk_size=200, overlap=50, threshold=30):
    """
    Performs alignment on overlapping chunks to handle large sequences.
    Returns list of alignment regions with their scores and positions.
    """
    alignments = []
    step = chunk_size - overlap

    print(f"Analyzing sequences: {len(seq1)}bp vs {len(seq2)}bp")
    print(f"Chunk size: {chunk_size}, Overlap: {overlap}")

    chunk_count = 0
    for i in range(0, len(seq1) - chunk_size + 1, step):
        chunk1 = seq1[i : i + chunk_size]

        for j in range(0, len(seq2) - chunk_size + 1, step):
            chunk2 = seq2[j : j + chunk_size]

            a1, a2, score = smith_waterman(chunk1, chunk2)

            chunk_count += 1

            if score > threshold:
                # Calculate all three scoring equations
                score1 = scoring_equation_1(a1, a2)
                score2 = scoring_equation_2(a1, a2)
                score3 = scoring_equation_3(a1, a2)

                alignments.append(
                    {
                        "seq1_start": i,
                        "seq1_end": i + chunk_size,
                        "seq2_start": j,
                        "seq2_end": j + chunk_size,
                        "score": score,
                        "score1_percent_identity": score1,
                        "score2_similarity_index": score2,
                        "score3_normalized_score": score3,
                        "aligned1": a1,
                        "aligned2": a2,
                        "chunk1": chunk1,
                        "chunk2": chunk2,
                    }
                )

    print(f"Processed {chunk_count} chunk pairs")
    print(f"Found {len(alignments)} significant alignments (score > {threshold})")

    alignments.sort(key=lambda x: x["score"], reverse=True)

    return alignments[:15]


def show_alignment_viewer(aln):
    """
    Interactive alignment viewer with scrolling capability and scoring metrics.
    """
    a1 = aln["aligned1"]
    a2 = aln["aligned2"]
    score = aln["score"]

    fig = plt.figure(figsize=(14, 8))
    fig.canvas.manager.set_window_title("Alignment Viewer")

    # Info section with all scores
    ax_info = fig.add_axes([0.1, 0.85, 0.8, 0.1])
    ax_info.axis("off")
    info_text = f"SW Score: {score}  |  Score1 (% Identity): {aln['score1_percent_identity']:.2f}%  |  Score2 (Similarity): {aln['score2_similarity_index']:.2f}  |  Score3 (Normalized): {aln['score3_normalized_score']:.2f}"
    ax_info.text(
        0.5,
        0.5,
        info_text,
        ha="center",
        fontsize=11,
        weight="bold",
        bbox=dict(boxstyle="round", fc="#e6f2ff"),
    )

    ax_align = fig.add_axes([0.05, 0.3, 0.9, 0.45])
    ax_align.axis("off")

    window_size = 80

    txt_seq1 = ax_align.text(
        0, 0.75, "", fontfamily="monospace", fontsize=11, color="#0066cc"
    )
    txt_match = ax_align.text(
        0, 0.5, "", fontfamily="monospace", fontsize=11, color="#009900"
    )
    txt_seq2 = ax_align.text(
        0, 0.25, "", fontfamily="monospace", fontsize=11, color="#cc0000"
    )

    ax_align.text(
        -0.01, 0.75, "Flu:", ha="right", fontsize=11, weight="bold", color="#0066cc"
    )
    ax_align.text(
        -0.01, 0.25, "Cov:", ha="right", fontsize=11, weight="bold", color="#cc0000"
    )

    def update_display(pos):
        start = int(pos)
        end = min(start + window_size, len(a1))

        sub1 = a1[start:end]
        sub2 = a2[start:end]

        match_line = ""
        for c1, c2 in zip(sub1, sub2):
            if c1 == c2 and c1 != "-":
                match_line += "|"
            elif c1 == "-" or c2 == "-":
                match_line += " "
            else:
                match_line += " "

        txt_seq1.set_text(sub1)
        txt_match.set_text(match_line)
        txt_seq2.set_text(sub2)
        fig.canvas.draw_idle()

    update_display(0)
    plt.show()


def visualize_genome_comparison(alignments, len1, len2):
    """
    Visualizes genome alignments as connecting lines between two sequences.
    Click on a line to view detailed alignment.
    """
    fig, ax = plt.subplots(figsize=(14, 7))
    plt.subplots_adjust(bottom=0.15)

    ax.plot(
        [0, len1],
        [1.5, 1.5],
        linewidth=8,
        color="#0066cc",
        label="Influenza Genome",
        solid_capstyle="round",
    )
    ax.plot(
        [0, len2],
        [0.5, 0.5],
        linewidth=8,
        color="#cc0000",
        label="COVID-19 Genome",
        solid_capstyle="round",
    )

    line_data = {}

    for idx, aln in enumerate(alignments):
        mid1 = (aln["seq1_start"] + aln["seq1_end"]) / 2
        mid2 = (aln["seq2_start"] + aln["seq2_end"]) / 2

        alpha = min(0.7, aln["score1_percent_identity"] / 100)
        color = plt.cm.viridis(aln["score1_percent_identity"] / 100)

        (line,) = ax.plot(
            [mid1, mid2], [1.5, 0.5], color=color, alpha=alpha, linewidth=2, picker=5
        )
        line_data[line] = aln

    ax.set_ylim(0, 2)
    ax.set_xlim(-50, max(len1, len2) + 50)
    ax.set_yticks([0.5, 1.5])
    ax.set_yticklabels(["COVID-19", "Influenza"])
    ax.set_xlabel("Genome Position (bp)", fontsize=12)
    ax.set_title(
        "Genome Alignment Visualization\nClick on connection lines to view detailed alignment",
        fontsize=14,
        weight="bold",
    )
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    # Interaction handler
    def on_click(event):
        if event.artist in line_data:
            aln = line_data[event.artist]
            print(
                f"\nOpening alignment: SW Score={aln['score']}, Score1={aln['score1_percent_identity']:.2f}%, Score2={aln['score2_similarity_index']:.2f}, Score3={aln['score3_normalized_score']:.2f}"
            )
            show_alignment_viewer(aln)

    fig.canvas.mpl_connect("pick_event", on_click)

    sm = plt.cm.ScalarMappable(
        cmap=plt.cm.viridis, norm=plt.Normalize(vmin=0, vmax=100)
    )
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation="horizontal", pad=0.08, aspect=30)
    cbar.set_label("Sequence Identity (%)", fontsize=11)

    plt.show()


if __name__ == "__main__":
    seqFlu = """GCATATCATGTGGGTAAGCCGAAACTCGTGTTCATCATTTTTCGCCTCTTGTTACATAA
TCAAGCAGCTTTGTTAATAATTGGGCCATGCCGTTAAAAGGAGAAATCTCGTATACTATCCCTAAACATC
ATTATAATCCAATGGTAGAGCATCGTTTTAGAAGACCCCTGAGTACTTTCTCCATTAACCGTTCGTGCGT
CCTGGACAACGTTTGTGTTGGCAATGTCTGTACATCGTGTGGCGAGTTTGTTTGACGGTAATACGAACAT
CAGAAGTATGTGACCCAGTTGAGTAGTCTCTTGTATATTCCGACGGTAGTCTGAGAAGGTTTGTTTAAGT
AACCTAAATCGTAAGTGGTTGTTTCCATTTGATCCTTTGTATGAGAACTGTTCATCATAGAATACAACAT
TGCATCCAATACAATGTCTCTTATTACGTCCCCTGTCAAGAGCATCTGCTTGTCAGCATTGATTCTCTCG
CTAAGCTTCTGTTGGCATATGCACCGTATACAAGATAATCGCTCTTGAAGTTGATCCTGGTCGGACTAGT
GAGGTAAGGATCTTCCGTGAGTATATAGTTATGATCTGTCTAATAAAATCGTAACACTTACGTCTTCAAA
TTGGAGATACACAATTATTGCTTGTGTGACTTATTGCGGCAATTGTAATATTATGGACGGGAGTACTTTC
TTTCTTTTATTATGATTAACACAAGAATGAATTCACAACTGATGTTCTAGAGTTTGTGACTGATTTAAAA
CGTAAGGGCATAGACCACTTTATCACCATTATTCTGATTATGTCGCATTACTTACTCCGATTAATACCGC
TTTTATGAACTTCGCGTTATAAAGTGTGCCACAATCCCCAGTTACAGTATGCTATCGCGCGTATTTTACT
TTGAACAGAAAATTTTTGGTTTCCTAGCATCATTTGGCCGTTAAAAAAGAGATCCGATATCGTGAAGTAC
TGCCTCCTTCGAGGACTTGGGCTCTTCCGGGTGTTTGGTTCATT""".replace("\n", "")

    seqCov = """ATGGATGTCAATCCGACTTTACTTTTCTTGAAAGTGCCAGTGCAAAATGCTATAAGTACCACTTTCCCTT
    ATACTGGAGACCCTCCATACAGCCATGGAACAGGAACAGGATACACCATGGACACAGTCAACAGAACACA
    TAAATACTCAGAAAAAGGAAAGTGGACAACGAACACAGAGACTGGAGCACCCCAACTCAATCCAATTGAT
    GGACCATTACCTGAGGACAACGAGCCGAGTGGGTATGCACAAACGGATTGTGTATTGGAAGCAATGGCTT
    TCCTTGAAGAATCTCACCCAGGGATCTTTGAAAACTCGTGTCTCGAAACGATGGAAATTGTTCAGCAAAC
    AAGAGTGGATAAACTGACCCAAGGCCGCCAGACCTATGACTGGACGTTGAATAGAAATCAGCCGGCTGCT
    ACCGCATTGGCCAACACTATAGAGGTATTCAGATCGAATGGCCTGACAGCCAATGAATCAGGAAGGTTGA
    TCGATTTCCTCAAGGACGTGATGGATTCAATGGATAAGGAAGAAATGGAGATTACAACACATTTCCAGAG
    GAAGAGGAGAGTGAGGGACAACATGACCAAGAAAATGGTCACACAGAGAACAATAGGAAAGAAAAAACAA
    AGACTGAACAAAAGGAGCTACCTAATAAGAGCACTTACATTGAACACAATGACAAAGGATGCTGAAAGAG
    GCAAGCTGAAAAGGAGGGCAATCGCAACACCCGGGATGCAAATCAGAGGATTCGTGTATTTTGTAGAAGC
    ACTAGCGAGGAGCATCTGTGAGAAACTTGAGCAATCTGGCCTCCCTGTCGGAGGGAATGAGAAGAAAGCT
    AAATTGGCAAATGTTGTGAGGAAGATGATGACTAATTCACAAGATACAGAGCTCTCCTTCACAATTACTG
    GGGACAACACCAAATGGAATGAGAATCAAAACCCCCGGATGTTTCTAGCAATGATAACATACATCACAAG
    AAACCAGCCAGAATGGTTTAGAAATGTCTTAAGCATTGCTCCTATAATGTTCTCAAACAAGATGGCGAGA
    TTAGGAAAAGGGTACATGTTCGAAAGTAAGAGTATGAAGTTACGGACACAAGTACCAGCGGAAATGCTCG
    CAAATATTGACCTGAAATACTTCAACAAATCAACAAGAGAGAAAATCGAGAAAATAAGACCTCTACTGAT
    AGATGGCACAGCCTCATTGAGTCCTGGAATGATGATGGGCATGTTCAACATGTTGAGTACAGTCTTAGGA
    GTTTCAATTCTGAATCTCGGGCAGAAGAAGTACACCAAAACCACATATTGGTGGGACGGACTCCAATCCT
    CAGATGACTTCGCCCTCATAGTGAATGCACCGAATCATGAGGGAATACAGGCAGGAGTAGATAGGTTCTA
    TAGAACCTGCAAATTAGTTGGGATAAACATGAGCAAGAAGAAATCCTACATAAATCGGACAGGAACATTC
    GAATTCACAAGCTTTTTCTACCGCTATGGATTCGTAGCTAACTTCAGTATGGAGTTGCCCAGTTTTGGAG
    TGTCCGGGATTAATGAGTCAGCTGACATGAGCGTTGGTGTTACAGTAATAAAGAACAATATGATAAACAA
    CGATCTTGGACCAGCAACAGCCCAAATGGCCCTTCAGCTATTTATCAAAGACTACAGATACACATACCGA
    TGTCACAGGGGTGATACGCAAATTCAAACGAGGAGAGCATTCGAGCTGAAGAAGCTGTGGGAGCAGACCC
    GTTCGAAGGCAGGACTGTTGGTTTCAGATGGAGGGCCAAACCTGTACAATATCCGGAACCTCCACATTCC
    AGAGGTCTGCTTGAAATGGGAATTGATGGATGAAGACTACCAAGGCAGGTTGTGTAATCCTATGAACCCG
    TTTGTCAGTCATAAGGAAATTGATTCAGTCAACAATGCTGTGGTGATGCCAGCTCATGGCCCAGCCAAAA
    GCATGGAGTATGATGCCGTTGCAACCACACATTCATGGATTCCTAAGAGGAATCGCTCCATTCTCAACAC
    CAGCCAAAGGGGGATTCTTGAGGACGAACAGATGTACCAGAAGTGCTGCAACCTATTCGAAAAGTTCTTC
    CCCAGCAGTTCGTACAGGAGGCCAGTTGGAATTTCCAGCATGGTGGAGGCCATGGTGTCTAGGGCCCGAA
    TTGATGCACGAATTGACTTCGAATCTGGAAGGATTAAGAAAGAAGAGTTTGCTGAGATCATGAAGATCTG
    TTCCACCATTGAAGAGCTCAGACGGCAAAAATAGTGAATTTAGCTTGTCCTTCATGA""".replace("\n", "")

    print("GENOME LOCAL ALIGNMENT ANALYSIS")

    alignments = chunked_alignment(
        seqFlu, seqCov, chunk_size=50, overlap=30, threshold=40
    )

    if alignments:
        print("\nSCORING METRICS FOR TOP ALIGNMENT:")
        print(f"\nSmith-Waterman Score: {alignments[0]['score']}")
        print(
            f"Score 1 (Percent Identity): {alignments[0]['score1_percent_identity']:.2f}%"
        )
        print(
            f"Score 2 (Similarity Index): {alignments[0]['score2_similarity_index']:.2f}"
        )
        print(
            f"Score 3 (Normalized Score): {alignments[0]['score3_normalized_score']:.2f}"
        )
        visualize_genome_comparison(alignments, len(seqFlu), len(seqCov))
    else:
        print("\nNo significant alignments found. Try adjusting parameters.")
