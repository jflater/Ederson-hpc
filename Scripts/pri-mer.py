#python pri-mer.py 
#use this to find unique ends of kmers 

import screed, sys

#Breaks down a sequence input into corresponding k-mers
def rolling_window(seq, window_size):
   for i in xrange(len(seq) - window_size + 1):
      yield seq[i:i+window_size]


#Make fname be first and last 19 bp or two fnames forward and reverse
for record in screed.open(sys.argv[1]): #reading in your kmer file
   my_seq = record.sequence
   #print my_seq
   for n, x in enumerate(rolling_window(my_seq, 19)):
      if n == 0:
         print ">" + str(n) + "_" + record.name + "_forward"
         print x
      if n == 131:
         print ">" + str(n) + "_" + record.name + "_reverse"
         print x


      

