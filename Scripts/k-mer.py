#python k-mer.py [reference1.fa] [comparison1.fa] [comparison.fa] [comparison3.fa] > filename.fa
import screed, sys

def consume_genome(fname):
   genome = ''
   for f in fname:
      for record in screed.open(f):
         genome = genome + record.sequence
   return genome

def kmer_count(g, bp):
   f={}
   for x in range(len(g)+1-bp):
      kmer=g[x:x+bp]
      if f.has_key(kmer):
         continue
      else:
         f[kmer]=f.get(kmer,0)+1
   return(f)

def rolling_window(seq, window_size):
   for i in xrange(len(seq) - window_size + 1):
      yield seq[i:i+window_size]

fname = sys.argv[1:2]
fname_compare = sys.argv[2:]


ref_genome1 = consume_genome(fname)
kmer_dict = kmer_count(ref_genome1, 150)
ref_genome2 = consume_genome(fname_compare)

for seq in rolling_window(ref_genome2, 150):
	if kmer_dict.has_key(seq):
		del kmer_dict[seq]

for n, kmer in enumerate(kmer_dict.keys()):
	print ">" + str(n) + "_" + fname[0]
	print kmer 



