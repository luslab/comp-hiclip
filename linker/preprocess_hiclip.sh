#/bin/bash

# Preprocess data from Sugimoto et al., 2015
# A. M. Chakrabarti
# 12th June 2020, updated 16th September 2020, updated 28th January 2021
# Updated 2nd December 2022

conda activate comp-hiclip-dev

# Paths
DATADIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/data
GITHUBDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/comp-hiclip

# Make preprocessing directory
mkdir -p /camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/preprocessed
cd /camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/preprocessed

# Download raw data from EBI
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR605/ERR605257/ERR605257.fastq.gz -o $DATADIR/ERR605257.fastq.gz

# Move UMI to header and demultiplex
umi_tools extract -p NNNXXXXNN -I $DATADIR/ERR605257.fastq.gz -S ERR605257.umi.fastq.gz
reformat.sh trimreaddescription=t in=ERR605257.umi.fastq.gz out=ERR605257.umi.trd.fastq.gz # Need to remove everything after first whitespace in read name
cutadapt -j 8 -e 0 --no-indels -g LigPlusLow="^AATA" -g LigPlusHigh="^GGTT" -g LigMinus="^GGCG" -o {name}.fastq.gz ERR605257.umi.trd.fastq.gz > demux.cutadapt.log 2>&1

rm ERR605257.umi.fastq.gz
rm ERR605257.umi.trd.fastq.gz
rm LigMinus.fastq.gz
rm unknown.fastq.gz

# Trim adapters and identify linker reads
A=AGATCGGAAGAGCGGTTCAG
B=CTGTAGGCACCATACAATG

for i in LigPlusHigh.fastq.gz LigPlusLow.fastq.gz; do

    echo ${i%%.*}

    # Trim 3' BA concatemers and A
    cutadapt -j 8 -a $B$A -m 12 -n 10 -o ${i%%.*}.BA_trimmed.fastq.gz $i > ${i%%.*}_BA.cutadapt.log
    cutadapt -j 8 -a $A -m 12 -n 10 -o ${i%%.*}.A_BA_trimmed.fastq.gz ${i%%.*}.BA_trimmed.fastq.gz > ${i%%.*}_A_BA.cutadapt.log

    # Identify reads with B in the middle and at least 12 nt arms either side
    $GITHUBDIR/linker/linker.R ${i%%.*}.A_BA_trimmed.fastq.gz ${i%%.*}.linker.fastq.gz

    # Filter all adapter removed FASTQ for linker read ids
    zcat ${i%%.*}.linker.fastq.gz | \
    grep '^@R.' | \
    sed 's/^@R.//' \
    > ${i%%.*}.hits.txt

    filterbyname.sh in=${i%%.*}.A_BA_trimmed.fastq.gz out=${i%%.*}.nolinker.fastq.gz names=${i%%.*}.hits.txt include=f

    rm ${i%%.*}.BA_trimmed.fastq.gz
    rm ${i%%.*}.A_BA_trimmed.fastq.gz
    rm ${i%%.*}.hits.txt

done