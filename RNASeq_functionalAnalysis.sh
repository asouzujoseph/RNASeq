cd /mnt/storage/$USER/jupyternotebooks/RNASeq

head deseq.results.tsv

head -1 deseq.results.tsv  # show the headers
grep -n CDKN1A deseq.results.tsv
grep -n BBC3 deseq.results.tsv
grep -n GDF15 deseq.results.tsv

awk '$3 != "NA" && $3 > 1 && $7 < 0.05 {print $1}' deseq.results.tsv > upRegGenes.txt

awk '$3 != "NA" && $3 < -1 && $7 < 0.05 {print $1}' deseq.results.tsv > downRegGenes.txt

wc -l upRegGenes.txt
wc -l downRegGenes.txt

cat deseq.results.tsv | sort -k 3,3gr | awk '$3 != "NA" {print $1}' | grep -v Gene > rankedGenes.txt

cat deseq.results.tsv | sort -k3,3g | awk '$3 != "NA" {print $1}' | grep -v Gene > deseq.results.sortFCasc.txt

head deseq.results.sortFCasc.txt 

# extract the gene name and logFCvalues for all genes. Write output to a file
cat deseq.results.tsv | sort -k 3,3gr | awk '$3 != "NA" {print $1, $3}' | grep -v Gene | tr ' ' '\t' > deseq.logFC.rnk

head deseq.logFC.rnk


