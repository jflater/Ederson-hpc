"""Takes in kmers and removes the front and back 19 bp"""
#python pri-mer.py

import sys
import screed

#Breaks down a sequence input into corresponding k-mers
def rolling_window(seq, window_size):
    """rolling window slides along a sequence, in this case showing 19 bp at a time"""
    for i in xrange(len(seq) - window_size + 1):
        yield seq[i:i+window_size]

def main():
    """main for k-mer.py, find kmers in a .fa file."""
    #reading in your kmer file
    for record in screed.open(sys.argv[1]):
        my_seq = record.sequence
    #print my_seq
        for i, seq in enumerate(rolling_window(my_seq, 19)):
            if i == 0:
                print ">" + str(i) + "_" + record.name + "_forward"
                print seq
            if i == 131:
                print ">" + str(i) + "_" + record.name + "_reverse"
                print seq

if __name__ == '__main__':
    main()

