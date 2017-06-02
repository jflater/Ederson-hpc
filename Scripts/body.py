#python body.py
#use this to feed the body (150-(19*2))bp into a body.fa. We need to make sure no primers match any portion of the k-mer. In addtion, we also include this with all other .fa files we are comparing to.

import sys
import screed

def rolling_window(seq, window_size):
    """Breaks down a sequence input into corresponding k-mers"""
    for i in xrange(len(seq) - window_size + 1):
        yield seq[i:i+window_size]

#exclude first and last 19 bp (future primer) collect only center of 150 bp sequence
for record in screed.open(sys.argv[1]):
    my_seq = record.sequence
    for n, x in enumerate(rolling_window(my_seq, 112)):
        if n == 19:
            print ">" + str(n) + "_" + record.name
            print x
