"""use this to feed the body (150-(19*2))bp into a body.fa. We need to make sure no primers match any portion of the k-mer. In addtion, we also include this with all other .fa files we are comparing to.
"""
#python body.py

import sys
import screed

def rolling_window(seq, window_size):
    """Breaks down a sequence input into corresponding k-mers"""
    for i in xrange(len(seq) - window_size + 1):
        yield seq[i:i+window_size]


def main():
    """main for body.py"""
    for record in screed.open(sys.argv[1]):
        my_seq = record.sequence
        for i, seq in enumerate(rolling_window(my_seq, 112)):
            if i == 19:
                print ">" + str(i) + "_" + record.name
                print seq

if __name__ == '__main__':
    main()
