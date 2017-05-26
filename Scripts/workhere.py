#python k-mer.py [primer.fa] [comparison1.fa] [comparison.fa] [comparison3.fa] > filename.fa

#Need to take in primers forward and back, and validate them by making sure they do not appear in know genomes of non-target organisms. 

import screed, sys

def make_kmer_dict(fname):
   d = {}
   for f in fname:
       for record in screed.open(f):
	   if d.has_key(record.sequence):
   	       d[record.sequence].append(record.name)
	   else:
               d[record.sequence] = [record.name]
   return d

#makes dictionary of all k-mers in a sequence that are present
def kmer_count(g, bp):
   f={}
   for x in range(len(g)+1-bp):
      kmer=g[x:x+bp]
      if f.has_key(kmer):
         continue
      else:
         f[kmer]=f.get(kmer,0)+1
   return(f)

#breaks a sequence into its respective k-mers
def rolling_window(seq, window_size):
   for i in xrange(len(seq) - window_size + 1):
      yield seq[i:i+window_size]

fname = sys.argv[1:2]
fname_compare = sys.argv[2:]

kmer_dict = make_kmer_dict(fname)
#put consume_genome
#kmer_dict is a dictionary of all 19-mers in fnames
ref_genome2 = consume_genome(fname_compare)
#breakdown fname_compare into 19-mers
for seq in rolling_window(ref_genome2, 19):
        if kmer_dict.has_key(seq):
        	continue
	else:
		print ">"+kmer_dict[seq][0]
		print seq
#for record in screed.open(sys.argv[1]):
#        y = record.name

#for n, kmer in enumerate(kmer_dict.keys()):
#        print ">" + str(n) + "_" + fname[0]
#        print kmer
