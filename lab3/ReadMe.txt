Task 1:
The melting temp (tm) is the temp at which one half of a particular DNA will disociate and become a single stand of DNA
 primeland and seq are critical importnace in designing the parameters of a successful amplification. The tm of a nucleic acid 
 duplex increases both with its length and with increasing G C content.
 tm = 4(G + C) + 2(A + T)
 tm = 81.5 + 16.6(log10([Na+])) + 41*(%GC) - 600/length
 Implement an app that computes the tm of a DNA seq by using one of those formulas or both of them. 
 Input = a string of DNA, Output = temp in C

 Task 2:
 Design an app that uses a sliding window method in order to read the tm over the seq S use a sliding winow of 8 positions and chose
 a FASTA file as input.