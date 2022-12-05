#/bin/bash

# Map linker data from Sugimoto et al., 2015
# 01 February 2021

conda activate comp-hiclip-dev

DATADIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/preprocessed
REFDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref
GITHUBDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/comp-hiclip
TOSCADIR=~/.nextflow/assets/amchakra/tosca

# ==========
# Generate STAR index
# ==========

# cd $REFDIR

# STAR \
# --runThreadN 8 \
# --runMode genomeGenerate \
# --genomeDir STAR_GRCh38_GencodeV33_masked \
# --genomeFastaFiles GRCh38.gencode_v33.fa

# ==========

mkdir -p /camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/results_linker
cd /camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/results_linker

for i in stau1_high.linker.fastq.gz stau1_low.linker.fastq.gz; do

    echo ${i%%.*}

    # Map
    STAR \
	--runThreadN 8 \
	--runMode alignReads \
	--genomeDir $REFDIR/STAR_GRCh38_GencodeV33_masked \
	--readFilesIn $DATADIR/$i \
	--readFilesCommand gunzip -c \
	--genomeLoad NoSharedMemory \
	--outFileNamePrefix ${i%%.*}. \
	--outSAMattributes All \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterType BySJout \
	--alignIntronMin 20 \
	--alignIntronMax 100000 \
	--limitBAMsortRAM 60000000000 \
	--twopassMode Basic \
	--outFilterMultimapNmax 1 \
	--outFilterMultimapScoreRange 1 \
	--outFilterScoreMin 10

	rm ${i%%.*}.Log.progress.out 
	rm ${i%%.*}._STARgenome/*
	rm ${i%%.*}._STARpass1/*
	rmdir ${i%%.*}._STARgenome
	rmdir ${i%%.*}._STARpass1

	# Get hybrid table
	Rscript --vanilla $GITHUBDIR/linker/bam_to_toscatable.R \
	${i%%.*}.Aligned.sortedByCoord.out.bam \
	${i%%.*}

	# Deduplicate using UMIs
	$TOSCADIR/bin/deduplicate_hybrids.py \
	${i%%.*}.hybrids.tsv.gz \
	${i%%.*}.hybrids.dedup.tsv.gz \
	_ \
	directional \
	> ${i%%.*}.dedup.log

	# Process hybrids.dt
	Rscript --vanilla $GITHUBDIR/linker/process_toscatable.R \
	${i%%.*}.hybrids.dedup.tsv.gz \
	$REFDIR/GRCh38.gencode_v33.tx.gtf.gz \
	$REFDIR/icount_mini_utr3/regions.gtf.gz \
	${i%%.*}

done
