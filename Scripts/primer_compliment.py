"""working on python loop to print the compliment of a primer"""
import sys
import screed
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna


for record in screed.open(sys.argv[1]):
    seq = record.sequence

my_seq = Seq(seq)
print my_seq.complement()
