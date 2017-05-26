
#get sequences from refseq using their ID
for x in ids.txt; do curl -s 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={$x}&rettype=fasta'; done
