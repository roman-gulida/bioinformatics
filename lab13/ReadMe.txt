Task 1:
A square matrix of an arbitrary size and the corresponding initial vector are given. Implement a software application that makes
a prediction on a total of 5 discrete steps, using this matrix and the corresponding vector.

Task 2:
Use a random DNA sequence of about 50 letters, use the sequence to compute the transition probabilities between letters. 
Output should be the transition matrix stored as a json file.

Task 3:
Use a random english text of about 300 letters(that implies spaces and punctuation), and compute the transition probabilities between words.
Store the transition matrix as a json file, for ease of implementation you could represent each new word by using 1 symbol of your choice(ASCII).

Task 3 output:
Text: The quick brown fox jumps over the lazy dog. The dog was sleeping under a tree. 
    A bird flew over the tree and landed on a branch. The fox saw the bird and tried to catch it. 
    The bird quickly flew away. The dog woke up and barked at the fox. The fox ran into the forest. 
    The forest was dark and quiet. The moon shone brightly in the sky.
Length: 351 characters

Number of unique words: 40

Word to Symbol Mapping:
the: !
quick: "
brown: #
fox: $
jumps: %
over: &
lazy: '
dog: (
was: )
sleeping: *
under: +
a: ,
tree: -
bird: .
flew: /
and: 0
landed: 1
on: 2
branch: 3
saw: 4
tried: 5
to: 6
catch: 7
it: 8
quickly: 9
away: :
woke: ;
up: <
barked: =
at: >
ran: ?
into: @
forest: A
dark: B
quiet: C
moon: D
shone: E
brightly: F
in: G
sky: H

Task 4:
Use the transition matrix from the json output in order to synthesize new sequences of text based on the transition matrix