# Ederson
# Identification of primer targets for identifying biomarkers for specific bacteria

###  By Jared Flater in collaboration with Ederson Jesus and Adina Howe


# Inputs: List of accession numbers of targets in Accession_numbers_for_Jared_and_Adina.xlsx
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
6. NZ_AQHN0100009 6.1
7. NZ_AQHN1000084.1 

----

The python script "k_mer.py" is a script that works through repython k-mer.py frence and comparison genomes that we choose and finds unique sequnces of lenght k (k-mer) for this script k = 150 bp. 

#  Download RefSoilpython k-mer.py  from https://figshare.com/articles/RefSoil_Database
----

#  First script in k-mer process:
```python
"""breaks up a sequence from a .fa file into kmerlength, and compares to kmers from other genomes
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
```
```{bash}
bash-3.2$ python Scripts/k_mer.py test/test_NC.fa test/test_NZ.fa test/test_compare.fa > 1_NC_kmer_test.fa
```
# For each file output looks like:
```{bash}
bash-3.2$ head 1_NC_kmer_test.fa 
>0_test/test_NC.fa
ACTATATCGTGGTTTGCACTGTGACGGGTTCGACGCAAGCCGGCATGGTCGTTGGTTTCGCCAAGGATGGACGACAGCGCAATGTGATCGGTATCGATGCTTCGGCAACCCCTCTCCAAGCCCAGTCGCAGGTGCTTAACATTGCCCGGC
>1_test/test_NC.fa
CCATCAGAACACGCGGCGTGATCTGCCAGGATACGCCCCATCTGTCCTTGCACCAGCCGCATTCGCTTTCCTGGCCGCCATTGCCGACAATGGCGTTCCAATAACGATCCGTCTCTTCCTGGTCCTCGGTGTAGATCTGAAACGAAAAGG
>2_test/test_NC.fa
ATACCAATAAGACTTTGAAAAAGCTGGCGGCCCGAGGGTTGCTCCTTTGGCAGGAGAGGGGCTGCGAGGTTCTGGATGGGAAGGGGCTTGCAGAACTTGCGGGCTGGGACGGTTTTGATTCGGGCAAGCGCCCCTTCATCTAGGCACTGC
>3_test/test_NC.fa
AAGGCAGGTAAAACGGGTGAAAGATTCCTGTCGATACTCATTGTCCTGCATAAAATCGACATTGCTGTCGTAGATCTCATATTTGAGGCGGTCTTCCACATATTCTTGCGGTCGAGATCGCTGCTGCTTTTGCATAGCAGTGAGCGGCCG
>4_test/test_NC.fa
TGCCGCATCGCCTCGATCGCCATGCAGCAGATGTCATTGCCGCGGCCAATGGTGAAGGTCAGGCCGTGGCCAGCAAGGCCCGGCTTGTCGGTGTCGAGGATGACATAAGCTGCCGAATAATCCGGATCCGGATTCATCGCATCCGAACCA
```
This generated a file that was over 1 gb in size, I wan't really expecting a file this large, however, it may be possible. To test this, I made three files... one .fa file that contained the first 150bp kmer + 5bp of random. The sconed file is just a .fa with zzzzzzzzzzzzz filling it...the final is a .fa with the same first 150 bp as the first file. 

When we compare 155pb_test.fa to 155zbp_test.fa and 150bp_test.fa, we should get 5 kmers:
```{bash}
jflater-air:test jaredflater$ python ../Scripts/k_mer.py 155pb_test.fa 155zbp_test.fa 150bp_test.fa
>0_155pb_test.fa
ATATCGTGGTTTGCACTGTGACGGGTTCGACGCAAGCCGGCATGGTCGTTGGTTTCGCCAAGGATGGACGACAGCGCAATGTGATCGGTATCGATGCTTCGGCAACCCCTCTCCAAGCCCAGTCGCAGGTGCTTAACATTGCCCGGCGGT
>1_155pb_test.fa
TATCGTGGTTTGCACTGTGACGGGTTCGACGCAAGCCGGCATGGTCGTTGGTTTCGCCAAGGATGGACGACAGCGCAATGTGATCGGTATCGATGCTTCGGCAACCCCTCTCCAAGCCCAGTCGCAGGTGCTTAACATTGCCCGGCGGTT
>2_155pb_test.fa
CTATATCGTGGTTTGCACTGTGACGGGTTCGACGCAAGCCGGCATGGTCGTTGGTTTCGCCAAGGATGGACGACAGCGCAATGTGATCGGTATCGATGCTTCGGCAACCCCTCTCCAAGCCCAGTCGCAGGTGCTTAACATTGCCCGGCG
>3_155pb_test.fa
TATATCGTGGTTTGCACTGTGACGGGTTCGACGCAAGCCGGCATGGTCGTTGGTTTCGCCAAGGATGGACGACAGCGCAATGTGATCGGTATCGATGCTTCGGCAACCCCTCTCCAAGCCCAGTCGCAGGTGCTTAACATTGCCCGGCGG
>4_155pb_test.fa
ATCGTGGTTTGCACTGTGACGGGTTCGACGCAAGCCGGCATGGTCGTTGGTTTCGCCAAGGATGGACGACAGCGCAATGTGATCGGTATCGATGCTTCGGCAACCCCTCTCCAAGCCCAGTCGCAGGTGCTTAACATTGCCCGGCGGTTC
```
Now, if we make the same comparison, less the final 150bp_test.fa file, we should have 6 kmers returned:

```{bash}
jflater-air:test jaredflater$ python ../Scripts/k_mer.py 155pb_test.fa 155zbp_test.fa
>0_155pb_test.fa
ACTATATCGTGGTTTGCACTGTGACGGGTTCGACGCAAGCCGGCATGGTCGTTGGTTTCGCCAAGGATGGACGACAGCGCAATGTGATCGGTATCGATGCTTCGGCAACCCCTCTCCAAGCCCAGTCGCAGGTGCTTAACATTGCCCGGC
>1_155pb_test.fa
ATATCGTGGTTTGCACTGTGACGGGTTCGACGCAAGCCGGCATGGTCGTTGGTTTCGCCAAGGATGGACGACAGCGCAATGTGATCGGTATCGATGCTTCGGCAACCCCTCTCCAAGCCCAGTCGCAGGTGCTTAACATTGCCCGGCGGT
>2_155pb_test.fa
TATCGTGGTTTGCACTGTGACGGGTTCGACGCAAGCCGGCATGGTCGTTGGTTTCGCCAAGGATGGACGACAGCGCAATGTGATCGGTATCGATGCTTCGGCAACCCCTCTCCAAGCCCAGTCGCAGGTGCTTAACATTGCCCGGCGGTT
>3_155pb_test.fa
CTATATCGTGGTTTGCACTGTGACGGGTTCGACGCAAGCCGGCATGGTCGTTGGTTTCGCCAAGGATGGACGACAGCGCAATGTGATCGGTATCGATGCTTCGGCAACCCCTCTCCAAGCCCAGTCGCAGGTGCTTAACATTGCCCGGCG
>4_155pb_test.fa
TATATCGTGGTTTGCACTGTGACGGGTTCGACGCAAGCCGGCATGGTCGTTGGTTTCGCCAAGGATGGACGACAGCGCAATGTGATCGGTATCGATGCTTCGGCAACCCCTCTCCAAGCCCAGTCGCAGGTGCTTAACATTGCCCGGCGG
>5_155pb_test.fa
ATCGTGGTTTGCACTGTGACGGGTTCGACGCAAGCCGGCATGGTCGTTGGTTTCGCCAAGGATGGACGACAGCGCAATGTGATCGGTATCGATGCTTCGGCAACCCCTCTCCAAGCCCAGTCGCAGGTGCTTAACATTGCCCGGCGGTTC
```
# Example shell commands and file names for first script:
----
Find unique k-mers in r.freire
```bash
python k-mer.py refrence_genomes/NZ_AQHN01000095.1.fa refrence_genomes/NZ_AQHN01000096.1.fa refrence_genomes/NC*.fa compare_genomes/*.fa > r.freirei.mers.fa
```
Find unique k-mers in r.tropici
```bash
phython k-mer.py refrence_genomes/NC_02006*.fa refrence_genomes/NC_020059.1.fa refrence_genomes/NZ*.fa compare_genomes/*.fa > r.tropici.mers.fa
```
# Now we will move onto the second script in the process, identifying potential primers for each k-mer, ensuring that each primer is only present in one gene.  This scripts input is the k-mers unique to each genome and output is the potential primers at each end of k-mer:

```bash
#python pri-mer.py
#use this to find unique ends of kmers

import screed, sys

#Breaks down a sequence input into corresponding k-mers
def rolling_window(seq, window_size):
   for i in xrange(len(seq) - window_size + 1):
      yield seq[i:i+window_size]


#Make fname be first and last 19 bp 
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
```
# This is the output from pri-mer.py on test_NZ_mers.fa
```bash
[flaterj1@dev-intel14 Scripts]$ python pri-mer.py ../test/test_NZ_mers.fa
>0_0_refrence_genomes/NZ.fa_forward
TCAGCCAGTCCGATCCAGA
>131_0_refrence_genomes/NZ.fa_reverse
TGGGGTCCGCAACTCAAGG
>0_1_refrence_genomes/NZ.fa_forward
GATGGGCGAAACCGCGGCG
>131_1_refrence_genomes/NZ.fa_reverse
GGGCTTATTGGATCGCCGT
>0_2_refrence_genomes/NZ.fa_forward
GCGATCACCCGGTGAAGCC
>131_2_refrence_genomes/NZ.fa_reverse
GAGGAGTTCCGTCAAAGAT
>0_3_refrence_genomes/NZ.fa_forward
TCATGATGGCGGCTGCGGT
>131_3_refrence_genomes/NZ.fa_reverse
ACCGTGCTGGTGCGTTCGC
>0_4_refrence_genomes/NZ.fa_forward
ACTATGGCCTTGCACCGGG
>131_4_refrence_genomes/NZ.fa_reverse
GATCCGCGCATAGGCGGCG
```
# Now we need to get the "body" out of each k-mer (the middle is 112 bp, 150 total bp minus the two 19 bp ends)
# For this, use body.py
```python
#python body.py
#use this to feed the body (150-(19*2))bp into a body.fa. We need to make sure no primers match any portion of the k-mer. In addtion, we also include this with all other .fa files we are comparing to.

import screed, sys

#Breaks down a sequence input into corresponding k-mers
def rolling_window(seq, window_size):
   for i in xrange(len(seq) - window_size + 1):
      yield seq[i:i+window_size]

#exclude first and last 19 bp (future primer) collect only center of 150 bp sequence
for record in screed.open(sys.argv[1]):
   my_seq = record.sequence
   for n, x in enumerate(rolling_window(my_seq, 112)):
      if n == 19:
         print ">" + str(n) + "_" + record.name
         print x
```
# This is what the output of body.py looks like:
```bash
[flaterj1@dev-intel14 Scripts]$ python body.py ../test/test_NZ_mers.fa
>19_0_refrence_genomes/NZ.fa
CTATTACTTCAGCCTGCGTCCCGTGGTCGAAGCCTTCCCCGACGCCCGCGTCATCGCTGCCAGCGCCACCATCGAGGCGATCAAGGCAAATGTGCAGAAGAAGCTCGACACC
>19_1_refrence_genomes/NZ.fa
GATTATCTCGCCGTCAACATGGGGCTTGGCCTTACGCTGACGTCGCTCGCCATGGCCGTAGTGCTTGTCGGTGTGCTTGCGATCCAGTTTCTGCAGGAGCGCTACGTTCCTT
>19_2_refrence_genomes/NZ.fa
TTGAGATGCGAGATCAGCACCGTGATGACGAAACCATCGATATCACGCGAAAATCGGAAGCCGTGATTGAGGCGTGGCGGCACCGTAATGATCGCGGGTGGGGTGATGGCGT
>19_3_refrence_genomes/NZ.fa
GGTGCGTGGCTCGTCGGCCCAGACATAGCCGACCATCAGCGATGCATCCCAGACGCCATCGATGGCGTCGATCCCGGGCAGCATGTCGTAGAGGCCCTTGGCTGGTTCATCG
>19_4_refrence_genomes/NZ.fa
CAAGGACGAGAAGCCCGGTCTCTGGCTGGTTGGAGACCAAGGGATCTACGTCATGTCGAATGGAAGGCTGCGATCAGACGCCAGACCACTCGTGGTCTATGCGGAGGAATGC
```
# if we paste the k-mer, primer, and body into text edit, we can see that everything lines up how would like it to :
```bash
K-mer: >0_refrence_genomes/NZ.fa: 
TCAGCCAGTCCGATCCAGACTATTACTTCAGCCTGCGTCCCGTGGTCGAAGCCTTCCCCGACGCCCGCGTCATCGCTGCCAGCGCCACCATCGAGGCGATCAAGGCAAATGTGCAGAAGAAGCTCGACACCTGGGGTCCGCAACTCAAGG
primer:>0_0_refrence_genomes/NZ.fa_forward                                                                                          >131_0_refrence_genomes/NZ.fa_reverse
TCAGCCAGTCCGATCCAGA                                                                                                                TGGGGTCCGCAACTCAAGG
body:              CTATTACTTCAGCCTGCGTCCCGTGGTCGAAGCCTTCCCCGACGCCCGCGTCATCGCTGCCAGCGCCACCATCGAGGCGATCAAGGCAAATGTGCAGAAGAAGCTCGACACC
```
```python
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
```
# At this point we still need some test files to fully analyze this pipeline.
# Need: comparison genomes, refsoil, body
# Create: test file that is a copy of each of the above, but contains one primer. 
# This will serve as a false positive. 

##  Run subset fasta.py on all files in "full genomes" and "master genomes"
```bash
for file in compare_genomes/*.fa; do python subset_fasta.py $file >> Ederson-hpc/test/test_compare_genomes.fa; done
```
make sure file is empty or does not exist when running above ^
# Handy script to return odd lines in a file:
```bash
sed -n 1~2p file
```
This will start at line 1 and print every 2 lines, odd ones. If this were a fasta file, it would only print the carrot lines!
```bash
[flaterj1@dev-intel14 test]$ sed -n 1~2p test_compare.fa 
[flaterj1@dev-intel14 test]$ head test_compare.fa
>NC_004041.2 Rhizobium etli CFN 42 plasmid symbiotic plasmid p42d, complete sequence
>NC_007761.1 Rhizobium etli CFN 42, complete genome
>NC_007762.1 Rhizobium etli CFN 42 plasmid p42a, complete sequence
>NC_007763.1 Rhizobium etli CFN 42 plasmid p42b, complete sequence
>NC_007764.1 Rhizobium etli CFN 42 plasmid p42c, complete sequence
>NC_007765.1 Rhizobium etli CFN 42 plasmid p42e, complete sequence
>NC_007766.1 Rhizobium etli CFN 42 plasmid p42f, complete sequence
>NC_008378.1 Rhizobium leguminosarum bv. viciae plasmid pRL12 complete genome, strain 3841
>NC_008379.1 Rhizobium leguminosarum bv. viciae plasmid pRL9 complete genome, strain 3841
>NC_008380.1 Rhizobium leguminosarum bv. viciae chromosome complete genome, strain 3841
```

# Need to investigate how to test a primer for it's reverse compliment

```bash
cat input.fa | while read L; do echo $L; read L; echo "$L" | rev | tr "ATGC" "TACG" ; done
```
# I seems that for this to work on a file with multiple lines, the .fa file should first be linearized:
```bash
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < input.fa
```
