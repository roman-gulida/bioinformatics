import math
import re
from collections import defaultdict

# Sample poetry by Mihai Eminescu (translated excerpt from "Luceafărul")
EMINESCU_TEXT = """
She was an emperor's daughter, with eyes of azure light,
Each evening at her window she watched the stars at night.
And one among them all she loved, the brightest in the sky,
The evening star that rose and set, and caught her tender eye.
She gazed upon that shining star with longing and desire,
And whispered to the evening star, her heart was all afire.
Come down to me, oh evening star, descend from heaven's height,
And walk with me through palace halls, beneath the pale moonlight.
"""

# Sample poetry by Nichita Stănescu (translated excerpt, similar style)
STANESCU_TEXT = """
The word flows like water through my hands and fingers,
Each syllable a drop of light that in the silence lingers.
I speak to you in whispers soft, in language made of dreams,
Where meaning bends and transforms, nothing is what it seems.
The world dissolves in metaphor, in symbols and in sound,
Where poetry becomes the place where truth is always found.
I write with ink of starlight, on pages made of air,
And every verse I craft for you is like a whispered prayer.
"""

# Mihai's disputed text (combination of both styles)
MIHAI_TEXT = """
She gazed upon the shining heavens with longing in her eyes,
Each evening at her window watching stars across the skies.
The word flows through the palace halls like water pure and bright,
Where meaning bends and transforms in the pale enchanted night.
Come down to me, she whispered soft, in language made of dreams,
The evening star that rose and set, nothing is what it seems.
Her heart was all afire with symbols and with sound,
I write with ink of starlight where truth is always found.
The brightest star in heaven's height descends through silver air,
And every verse she speaks to him is like a whispered prayer.
"""


def preprocess_text(text):
    """Convert text to lowercase and extract words"""
    # Remove punctuation and convert to lowercase
    text = text.lower()
    # Split into words
    words = re.findall(r"\b[a-z]+\b", text)
    return words


def build_transition_matrix(words):
    """Build word transition count matrix"""
    transitions = defaultdict(lambda: defaultdict(int))
    vocabulary = set(words)

    # Count transitions
    for i in range(len(words) - 1):
        word1 = words[i]
        word2 = words[i + 1]
        transitions[word1][word2] += 1

    return transitions, vocabulary


def calculate_probabilities(transitions, vocabulary):
    """Convert transition counts to probabilities"""
    probabilities = defaultdict(lambda: defaultdict(float))

    for word1 in transitions:
        total = sum(transitions[word1].values())
        if total > 0:
            for word2 in transitions[word1]:
                probabilities[word1][word2] = transitions[word1][word2] / total

    return probabilities


def build_log_likelihood_matrix(prob_eminescu, prob_stanescu, vocab):
    """Create log-likelihood matrix: log2(P_eminescu / P_stanescu)"""
    log_likelihood = defaultdict(lambda: defaultdict(float))

    for word1 in vocab:
        for word2 in vocab:
            p_eminescu = prob_eminescu[word1][word2]
            p_stanescu = prob_stanescu[word1][word2]

            # Handle edge cases
            if p_eminescu > 0 and p_stanescu > 0:
                ratio = p_eminescu / p_stanescu
                log_likelihood[word1][word2] = math.log(ratio) / math.log(2)
            elif p_eminescu > 0 and p_stanescu == 0:
                log_likelihood[word1][word2] = 5.0  # Strong Eminescu signal
            elif p_eminescu == 0 and p_stanescu > 0:
                log_likelihood[word1][word2] = -5.0  # Strong Stănescu signal
            else:
                log_likelihood[word1][word2] = 0.0  # Neither

    return log_likelihood


def analyze_with_sliding_window(words, log_likelihood, window_size=5):
    """Analyze text using sliding window approach"""
    results = []

    for i in range(len(words) - window_size + 1):
        window = words[i : i + window_size]
        score = 0.0
        transitions_in_window = []

        # Calculate score for this window
        for j in range(len(window) - 1):
            word1 = window[j]
            word2 = window[j + 1]
            transition_score = log_likelihood[word1][word2]
            score += transition_score
            transitions_in_window.append((word1, word2, transition_score))

        # Determine attribution
        if score > 1.0:
            attribution = "EMINESCU"
        elif score < -1.0:
            attribution = "STĂNESCU"
        else:
            attribution = "UNCLEAR/ORIGINAL"

        results.append(
            {
                "position": i,
                "window": " ".join(window),
                "score": score,
                "attribution": attribution,
                "transitions": transitions_in_window,
            }
        )

    return results


def generate_report(results, words):
    """Generate detailed plagiarism report"""
    print("PLAGIARISM ANALYSIS REPORT")

    print("\nEXECUTIVE SUMMARY")
    print("-" * 50)

    eminescu_count = sum(1 for r in results if r["attribution"] == "EMINESCU")
    stanescu_count = sum(1 for r in results if r["attribution"] == "STĂNESCU")
    original_count = sum(1 for r in results if r["attribution"] == "UNCLEAR/ORIGINAL")

    total = len(results)

    print(f"Total windows analyzed: {total}")
    print(
        f"Windows matching Eminescu style: {eminescu_count} ({100 * eminescu_count / total:.1f}%)"
    )
    print(
        f"Windows matching Stănescu style: {stanescu_count} ({100 * stanescu_count / total:.1f}%)"
    )
    print(
        f"Original/Unclear sections: {original_count} ({100 * original_count / total:.1f}%)"
    )

    print("\nDETAILED ANALYSIS (Sliding Window)")
    print("-" * 80)
    print(f"{'Pos':<5} {'Score':<8} {'Attribution':<20} {'Text Window'}")
    print("-" * 80)

    for result in results:
        score_str = f"{result['score']:+.2f}"
        print(
            f"{result['position']:<5} {score_str:<8} {result['attribution']:<20} {result['window'][:50]}..."
        )

    print("\nDETAILED TRANSITION ANALYSIS (First 10 windows)")

    for idx, result in enumerate(results[:10]):
        print(
            f'\nWindow {idx + 1} (Position {result["position"]}): "{result["window"]}"'
        )
        print(f"Overall Score: {result['score']:+.3f} → {result['attribution']}")
        print("Transitions:")
        for word1, word2, score in result["transitions"]:
            print(f"  {word1} → {word2}: {score:+.3f}")

    print("CONCLUSION")

    # Identify continuous plagiarized sections
    plagiarized_sections = []
    current_section = None

    for result in results:
        if result["attribution"] in ["EMINESCU", "STĂNESCU"]:
            if (
                current_section is None
                or current_section["author"] != result["attribution"]
            ):
                if current_section:
                    plagiarized_sections.append(current_section)
                current_section = {
                    "author": result["attribution"],
                    "start": result["position"],
                    "end": result["position"] + 5,
                    "text": result["window"],
                }
            else:
                current_section["end"] = result["position"] + 5
        else:
            if current_section:
                plagiarized_sections.append(current_section)
                current_section = None

    if current_section:
        plagiarized_sections.append(current_section)

    if plagiarized_sections:
        print("\nPLAGIARISM DETECTED!")
        print(
            f"\nThe defendant Mihai has incorporated {len(plagiarized_sections)} section(s)"
        )
        print("from the works of established poets:\n")

        for i, section in enumerate(plagiarized_sections, 1):
            word_range = f"words {section['start'] + 1}-{section['end']}"
            print(f"{i}. Source: {section['author']} ({word_range})")
            print(f'   Sample: "{section["text"]}..."')

        print("\nVERDICT:")
        print("GUILTY of plagiarism")
    else:
        print("\nNO SIGNIFICANT PLAGIARISM DETECTED")
        print("\nThe text appears to be original or sufficiently transformative.")
        print("NOT GUILTY")


def main():
    words_eminescu = preprocess_text(EMINESCU_TEXT)
    transitions_eminescu, vocab_eminescu = build_transition_matrix(words_eminescu)
    prob_eminescu = calculate_probabilities(transitions_eminescu, vocab_eminescu)

    words_stanescu = preprocess_text(STANESCU_TEXT)
    transitions_stanescu, vocab_stanescu = build_transition_matrix(words_stanescu)
    prob_stanescu = calculate_probabilities(transitions_stanescu, vocab_stanescu)

    # Combined vocabulary
    combined_vocab = vocab_eminescu.union(vocab_stanescu)

    log_likelihood = build_log_likelihood_matrix(
        prob_eminescu, prob_stanescu, combined_vocab
    )

    words_mihai = preprocess_text(MIHAI_TEXT)

    # Analyze with sliding window
    results = analyze_with_sliding_window(words_mihai, log_likelihood, window_size=5)

    # Generate report
    generate_report(results, words_mihai)


if __name__ == "__main__":
    main()
