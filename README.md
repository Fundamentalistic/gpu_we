# V-Fold
Source code of program for protein distance map prediction

Training dataset was taken from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5860114/

We used PyTorch as neural network framework

## Tutorial

Dataset was developed with agreement to language processing models

Input is sequence of amino acid and first amino acid is target to we try to find distance map

Input example:
- ADNCBRJCBANDJFHCBJA

Input words volume is $`20`$ and we build sentences of any length

Output sequence words is sequence number of amino acid which to contact with target and distance $`d \in \mathbb{N}`$ as angstroms

5,6 7,5

All output words count is sequence maximum sequence length and count of natural numbers between 0 and cutoff distance
