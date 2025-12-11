Task 1:
Implement a software application that aligns two DNA sequences by using the Needleman–Wunsch algorithm:

S1="ACCGTGAAGCCAATAC"
S2="AGCGTGCAGCCAATAC"

Screenshot_1

Task 2:
Download the influenza genome and the COVID-19 genome from NCBI and align them using a local alignment method.
 Note that you need to implement an intermediate layer solution in order to align the two genome files.
 The alignment should be performed step by step on large regions, with connections between the partial results.
Standard sequence alignment algorithms do not allow the alignment of very large sequences because the scoring matrix
 grows exponentially with sequence length. The main result should be a visualization of the similarities between the two genomes.

Note: Do not use shortcuts proposed by AI; use native code without any hidden modules.

Screenshot_2, Screenshot_3

Task 3:
Formulate 3 scoring equations that are able to show the level of similarity between the 2 seq. 
alligned in the prev assigment. implement each of this pre scoring equations in your current implementation.