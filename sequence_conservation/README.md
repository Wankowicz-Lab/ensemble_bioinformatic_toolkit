This is the initial analysis to identify similar sequences between highly related structures. The goal here is to build a script that will do the following"

1) Read in metadata
2) Create FASTA file
3) Run MAFFT alignment of the sequences (https://mafft.cbrc.jp/alignment/server/index.html)
4) Calculate a shannon entropy for each position (https://en.wikipedia.org/wiki/Entropy_(information_theory)) in the sequence and plot this (THINK ABOUT THIS).
5) Map the shannon entropy for each residue back to the residue in the PDB. 
