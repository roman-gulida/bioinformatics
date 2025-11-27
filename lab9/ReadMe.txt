Task 1:
1. Take an arbitrary DNA sequence from the NCBI (National Center for Biotechnology), between 1000 and 3000 nucleotides (letters).
2. Use 5 restriction enzymes (enzyme name, recognized sequence, cleavage site):

EcoRI     	
5'GAATTC
3'CTTAAG
     
5'---G     AATTC---3'
3'---CTTAA     G---5'

BamHI     	
5'GGATCC
3'CCTAGG
     
5'---G     GATCC---3'
3'---CCTAG     G---5'

HindIII     	
5'AAGCTT
3'TTCGAA
     
5'---A     AGCTT---3'
3'---TTCGA     A---5'

TaqI     	
5'TCGA
3'AGCT
     
5'---T   CGA---3'
3'---AGC   T---5'

HaeIII*     	
5'GGCC
3'CCGG
     
5'---GG  CC---3'
3'---CC  GG---5'

Input of the implementation:
1. The recognized sequence for each restriction enzyme.
2. A DNA sequence to be digested.

Output of the implementation:
1. Number of cleavages for a each restriction enzyme.
2. Cleavage positions and length of fragments.
3. A simulation of the electrophoresis gel based on the number of restriction enzymes used.

Task 2:
Download 10 influenza viruses varianths from NCBI and analyze their genome by using the application from 
assignment 1. 
A) Make an electrophoresis gel for each genome.
B) Eliminate all lines that are in common between the gel simulation, such that differences will be shown.
C) Merge all electrophoresis gel simulations(that shown only in the differences) in one genome electrophoresis gel.