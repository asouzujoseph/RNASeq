# create directory for the analysis
mkdir -p /mnt/storage/$USER/jupyternotebooks/RNASeq
cd /mnt/storage/$USER/jupyternotebooks/RNASeq

# Convert SRA files to fastq format using the code:
# fastq-dump --split-files SRR number where SRR Number represents the run identifier of the sequence
# NS1 = SRR867127, NS2 = SRR867128, S1 = SRR867129 and S2 = SRR867130
# Note that --split-files is useful for paired reads. 

#For this project, this can  be done using a for loop as follows:
#    for i in SRR867127 SRR867128 SRR867129 SRR867130
#    do 
#        fastq-dump -O . --split-files $i 
#    done

# For this project, I made a symbolic link to the already pregenerated fastq files
ln -sf /mnt/nfs/data/RNA-seq/MCF7/*.fastq .

# write a for loop to perform quality check on the fastq files
mkdir -p FastQC
for i in NS1.fastq NS2.fastq S1.fastq S2.fastq 
do
    fastqc -o ./FastQC $i
done

# We will use the human genome assembly - hg19

# Alignment for NS1 data
STAR --genomeDir /mnt/nfs/mfiers/STAR/hg19_star_db \
     --genomeLoad NoSharedMemory \
     --runThreadN 2 \
     --readFilesIn NS1.fastq \
     --outFileNamePrefix NS1.


# Alignment for NS2 data
STAR --genomeDir /mnt/nfs/mfiers/STAR/hg19_star_db \
     --genomeLoad NoSharedMemory \
     --runThreadN 2 \
     --readFilesIn NS2.fastq \
     --outFileNamePrefix NS2.

# Alignment for S2 data
STAR --genomeDir /mnt/nfs/mfiers/STAR/hg19_star_db \
     --genomeLoad NoSharedMemory \
     --runThreadN 2 \
     --readFilesIn S2.fastq \
     --outFileNamePrefix S2.

# Alignment for S1 data
STAR --genomeDir /mnt/nfs/mfiers/STAR/hg19_star_db \
     --genomeLoad NoSharedMemory \
     --runThreadN 2 \
     --readFilesIn S1.fastq \
     --outFileNamePrefix S1.

# 
head -34 S1.Aligned.out.sam | grep -v '^@'

samtools sort -o NS1.bam NS1.Aligned.out.sam

samtools sort -o NS2.bam NS2.Aligned.out.sam

samtools sort -o S1.bam S1.Aligned.out.sam

samtools sort -o S2.bam S2.Aligned.out.sam

for i in NS1.bam NS2.bam S1.bam S2.bam
do 
    samtools index $i
done

ls *.bai

samtools view S1.bam | head -3

# samtools idxstats shows many reads map to each chromosome. 
# The output shows chromosome name, sequence length, # mapped read-segments and # unmapped read-segments
samtools idxstats NS1.bam

# A symbolic link is made to the hg19 annotation
ln -sf /mnt/nfs/data/RNA-seq/gencode.v19.nopseudo.plus.sort.gtf
ls -l *gtf

# Counts the number of genomic features on the mapped reads. I used gene name as reference for counting the features
# this code couns the genomic features in the four BAM files
featureCounts -Q 10 -g gene_name -a gencode.v19.nopseudo.plus.sort.gtf -o gene.counts *.bam

# Non-useful information are removed from the counts file
cut -f1,7- gene.counts | grep -v '^#' > all.gene.counts

head all.gene.counts

grep BBC3 all.gene.counts     # upregulated gene ---> positive control
grep RRM2B all.gene.counts    # upregulated gene ---> positive control
grep CDKN1A all.gene.counts   # upregulated gene ---> positive control
grep FEN1 all.gene.counts      # downregulated gene ---> negative control
grep TIMELESS all.gene.counts  # downregulated gene ---> negative control
grep DEK all.gene.counts       # downregulated gene ---> negative control
