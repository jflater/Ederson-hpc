"""similar to k_mer.py, but for primers from those kmers"""
#python primer_validate.py [primer.fa] [comparison1.fa] [comparison.fa] [comparison3.fa] > name.fa

#To take in primers forward and back
#And validate them by making sure only appear in genome of interest

import sys
import screed

#Not sure on this script, needs testing, was formated from k_mer.py
def consume_genome(fname):
    """reads through a .fa file, returning names and sequences"""
    genome = ''
    for i in fname:
        for record in screed.open(i):
            genome = genome + record.sequence
    return genome

def make_kmer_dict(fname):
    """making an object with names and sequences of a .fa"""
    for i in fname:
        for record in screed.open(i):
            print record.name
            print record.sequence

def kmer_count(kmerlookingin, primerlength):
    """makes dictionary of primers"""
    primer19 = {}
    print kmerlookingin
    for i in range(len(kmerlookingin)+1-primerlength):
        kmer = kmerlookingin[i:i+primerlength]
        if primer19.has_key(kmer):
            continue
        else:
            primer19[kmer] = 1
    return primer19

#breaks a sequence into its respective k-mers
def rolling_window(seq, window_size):
    """rolling window slides along a sequence, showing 19 bp at a time"""
    for i in xrange(len(seq) - window_size + 1):
        yield seq[i:i+window_size]

def make_primer_dict(fname):
    """primer dict function"""
    i = {}
    for filename in fname:
        for record in screed.open(filename):
            st = record.sequence
            i[st] = record.name
            
    return i

def main():
    """main for primer validation"""
    fname = sys.argv[1:2]
    fname_compare = sys.argv[2:]

    kmer_dict = make_primer_dict(fname)
    #kmer_dict is a dictionary of all 19-mers in fnames
    ref_genome2 = consume_genome(fname_compare)
    #breakdown fname_compare into 19-mers
    for seq in rolling_window(ref_genome2, 19):
        if kmer_dict.has_key(seq):
            del kmer_dict[seq]

    for i, kmer in enumerate(kmer_dict.keys()):
        print ">" + str(i) + "_" + kmer_dict[kmer]
        print kmer

if __name__ == '__main__':
    main()
