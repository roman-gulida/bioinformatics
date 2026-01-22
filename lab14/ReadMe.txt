Task 1:
The CG island is a short stretch of DNA in which the frequency of the CG sequence is higher than other regions.
 It is also called the CpG island, where "p" simply indicates that "C" and "G" are connected by a phosphodiester bond.
CpG islands are often located around the promoters of housekeeping genes (which are essential for general cell functions) or
 other genes frequently expressed in a cell. At these locations, the CG sequence is not methylated. By contrast, the CG sequences
  in inactive genes are usually methylated to suppress their expression.

Two models are given in wich the observation (+) is tested against the null hypothesis (-):

M1 the (+) model (ISLAND)  
M2 the (-) model (NON ISLAND)

For M1 we have: S1 = “ATCGATTCGATATCATACACGTAT”, known to belong to a CpG island.
For M2 we have: S2 = “CTCGACTAGTATGAAGTCCACGCTTG”, known to belong to other regions in the genome.

Follow the steps below to implement a software application:
1. for the CpG+ model: Count the transition frequencies from the known sequence "S1" which does belong to a CpG island.
2. for the CpG- model: Also count the transition frequencies from the sequence "S2" which does not belong to a CpG island.
3. Make the log likelihood matrix

Note: on how to take log of any base if the log function of your programming language causes problems:

Given a new sequence S=“CAGGTTGGAAACGTAA”, use the log-likelihood matrix to find if sequence "S" belongs to a CpG island.

Task 2:
We are inside a courtroom, Mihai is accused of plagiarism. The lawyer must prove his innocence and asks you for help.
You take 2 poetries: one of M. Eminescu and the other is N. Stanescu. You use these 2 texts to create 2 transition martices, which is able 
to capture the probability of transition between words. You combine the matricse into a log likelihood matrix. We use the LLM to scan,
by using a sliding window, the text of Mihai which has been accused of plagiarism.

Negative scores -> one model
Positive scores -> other model
Zero scores -> neither

Simulate the text of Mihai by using chatgpt and ask it to make it a combination between 2 poetries of Eminescu and Stanescu.
At the end of this process you should be able to tell which parts of the text is written by Eminescu or by Stanescu or by neither.