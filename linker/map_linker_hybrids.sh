#/bin/bash

# Map linker data from Sugimoto et al., 2015
# 01 February 2021

# conda activate comp-hiclip-dev

cd /camp/lab/luscomben/home/users/iosubi/projects/comp_hiclip/linker

# DATADIR=/camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip/preprocessed
WORKDIR=/camp/lab/luscomben/home/users/iosubi/projects/comp_hiclip/linker


# generate STAR index

# STAR --runMode genomeGenerate --runThreadN 8 --genomeDir STAR_GRCh38_GencodeV33_masked \
# --genomeFastaFiles /camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip/ref/human.fa


for i in LigPlusHigh.linker.fastq.gz LigPlusLow.linker.fastq.gz; do

    echo ${i%%.*}

    # Map
    STAR --runThreadN 8 --runMode alignReads --genomeDir STAR_GRCh38_GencodeV33_masked \
	--readFilesIn $i \
	--readFilesCommand gunzip -c --genomeLoad NoSharedMemory --outFileNamePrefix ${i%%.*} \
	--outSAMattributes All --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
	--outFilterType BySJout --alignIntronMin 20 --alignIntronMax 100000 --limitBAMsortRAM 60000000000 \
	--twopassMode Basic --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 --outFilterScoreMin 10
	
	# Get hybrid table
	Rscript --vanilla $WORKDIR/bam_to_dataframe.R ${i%%.*}Aligned.sortedByCoord.out.bam ${i%%.*}
done
