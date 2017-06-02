"""breaks up a sequence from a .fa file into kmerlength, and compares to kmers from comparision genomes
"""
#python k-mer.py [reference1.fa] [comparison1.fa] [comparison.fa] [comparison3.fa] > filename.fa

import sys
import screed

def consume_genome(fname):
    """reads through a .fa file, returning names and sequences"""
    genome = ''
    for i in fname:
        for record in screed.open(i):
            genome = genome + record.sequence
    return genome

def kmer_count(genomelookingin, kmerlength):
    """makes a dictionary of 150 bp kmers"""
    kmer150 = {}
    for i in range(len(genomelookingin)+1-kmerlength):
        kmer = genomelookingin[i:i+kmerlength]
        if kmer150.has_key(kmer):
            continue
        else:
            kmer150[kmer] = kmer150.get(kmer, 0)+1
    return kmer150

def rolling_window(seq, window_size):
    """rolling window slides along a sequence, showing only 150 bp at a time"""
    for i in xrange(len(seq) - window_size + 1):
        yield seq[i:i+window_size]


def main():
    """main for k-mer.py, find kmers in a .fa file."""
    fname = sys.argv[1:2]
    fname_compare = sys.argv[2:]


    ref_genome1 = consume_genome(fname)
    kmer_dict = kmer_count(ref_genome1, 150)
    ref_genome2 = consume_genome(fname_compare)

    for seq in rolling_window(ref_genome2, 150):
        if kmer_dict.has_key(seq):
            del kmer_dict[seq]

    for i, kmer in enumerate(kmer_dict.keys()):
        print ">" + str(i) + "_" + fname[0]
        print kmer

if __name__ == '__main__':
    main()



