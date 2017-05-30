# Ederson
# Check size of .fa files in genome directories, some are less than 100b, indicating they lack actual sequences in the file and are most likely links to over view pages for that sequencing effort. 
# List of accesion numbers in Accession_numbers_for_Jared_and_Adina.xlsx
We are interested in IDing unique sequences in G1 and G2 that are not found in other rhizobium species as well as not found in RefSoil db or RefSeq db. Once we have those identified, we will id 19 by seq to target for primer design. 
![](Images/ideas.jpg)

----

# G1:Rhizobium tropici CIAT 899 (+chromosome) 
1. NC_020059.1+
2. NC_020060.1
3. NC_020061.1
4. NC_020062.1

----

# G2:Rhizobium freirei PRF 81
5. NZ_AQHN01000095.1
6. NZ_AQHN01000096.1
7. NZ_AQHN1000084.1 

----

The python script "k-mer.py" is a script that works through refrence and comparison genomes that we choose and finds unique sequnces of lenght k (k-mer) for this script k = 150 bp. 

Download RefSoil from https://figshare.com/articles/RefSoil_Database/4362812
----

First script in k-mer process:
```{python}
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
```
#Example shell commands and file names for first script:
----
Find k-mers in r.freire
```{bash}
python k-mer.py refrence_genomes/NZ_AQHN01000095.1.fa refrence_genomes/NZ_AQHN01000096.1.fa refrence_genomes/NC*.fa compare_genomes/*.fa > r.freirei.mers.fa
```
Find k-mers in r.tropici
```{bash}
phython k-mer.py refrence_genomes/NC_02006*.fa refrence_genomes/NC_020059.1.fa refrence_genomes/NZ*.fa compare_genomes/*.fa > r.tropici.mers.fa
```


