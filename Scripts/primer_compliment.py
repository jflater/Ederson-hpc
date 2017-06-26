"""Find reverse complement of reverse primer, print forward primer"""
# Input file should be in .fa format, with both forward and reverse primer
# the word "reverse" in record.name of each reverse primer

import sys
import screed

from Bio.Seq import Seq

def main():
    """this is the main, see above"""
    for record in screed.open(sys.argv[1]):
        seq = record.sequence
        my_seq = Seq(seq)
        if "reverse" in record.name:
            print ">" + record.name + "_complement"
            print my_seq.reverse_complement()
        if "forward" in record.name:        
            print ">" + record.name + "_complement"
            print my_seq.reverse_complement()

if __name__ == '__main__':
    main()


