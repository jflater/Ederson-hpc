#!/user/bin/python
"""This script subsets fasta files into a new file, population new file with the number specified"""

# Usage:
# python subset_fasta.py input1.fa <argument> > filename.fa

import sys
import screed

def main():
    """screed reading through a fasta file and print n sequences"""
    filename = sys.argv[1]
    num = int(sys.argv[2])
    with screed.open(filename) as seqfile:
        for i, read in enumerate(seqfile):
            if i < num:
                print '>%s\n%s' %(read.name, read.sequence)

if __name__ == '__main__':
    main()
