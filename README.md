# alignerviz - Visualize many sequences in three dimensions.

## About

Large genomic sequences are often difficult for the human eye to
visualize. In the paper *Visual Exploration of Genomic Data*, Michail
Vlachos, Bahar Taneri, Eamonn Keogh, and Philip S. Yu (2007) suggest
converting sequences to trajectories. Each nucleotide is a vector,
based on:

    A = (0, 1)
    T = (1, 0)
    C = (0, -1)
    G = (-1, 0)

and these vectors are connected beginning to end from an initial
starting point. A large sequence difference will quickly push a
trajectory in a different direction, while more similar sequences will
have similar trajectories.

## Requirements

 - Python (>= 2.5)
 - matplotlib
 - numpy

## Usage

Run alignerviz with:

    python aviz.py sequences.fasta

If there is a specific FASTA entry that needs to be highlighted,
specify with its header:

    python aviz.py -H seq_header_name sequences.fasta

Specific graphic options can be set too. For example, if there are
many very similar sequences it may be useful to specify an alpha, and
print each one in black:

    python aviz.py -b -a 0.4 sequences.fasta

Line width can also be specified:

    python aviz.py -l 3 sequences.fastq


