"""python k-mer.py [primer.fa] [/mnt/research/germs/jared/big_files/Genomes/compare_genomes/full_genomes.fa] [/mnt/research/germs/jared/big_files/Genomes/master_genomes.fa] [/mnt/research/germs/jared/big_files/Genomes/RefSoil_archaea_genomes.fa] [/mnt/research/germs/jared/big_files/Genomes/RefSoil_bacteria_genomes] [body.fa] > filename.fa"""

#Need to take in primers forward and back, and validate them by making sure they do not appear
#in know genomes of non-target organisms.

import sys
import screed

def consume_genome(fname):
    """takes multiple files and multiple genes (carrots) and returns the sequence as one variable"""
    genome = ''
    for i in fname:
        for record in screed.open(i):
            genome = genome + record.sequence
    return genome

def kmer_count(genomelookingin, primerlength):
    """makes dictionary of all k-mers in a sequence that are present"""
    kmer19 = {}
    for i in range(len(genomelookingin)+1-primerlength):
        kmer = genomelookingin[i:i+primerlength]
        if kmer19.has_key(kmer):
            continue
        else:
            kmer19[kmer] = kmer19.get(kmer, 0)+1
    return kmer19

def rolling_window(seq, window_size):
    """breaks a sequence into its respective k-mers"""
    for i in xrange(len(seq) - window_size + 1):
        yield seq[i:i+window_size]


def main():
    """this is the main, need it for python, google it"""
    # Remove hello, only useing for testing
    print "hello"
    fname = sys.argv[1:2]
    ref_genome1 = consume_genome(fname)
    kmer_dict = kmer_count(ref_genome1, 19)
    ##kmer_dict is a dictionary of all 19-mers in fnames
    fname_compare = sys.argv[2:]
    ref_genome2 = consume_genome(fname_compare)
    ##breakdown fname_compare into 19-mers
    for seq in rolling_window(ref_genome2, 19):
        if kmer_dict.has_key(seq):
            del kmer_dict[seq]

    #for record in screed.open(sys.argv[1]):
    #    y = record.name

    for i, kmer in enumerate(kmer_dict.keys()):
        print ">" + str(i) + "_" + fname[0]
        print kmer

if __name__ == '__main__':
    main()
