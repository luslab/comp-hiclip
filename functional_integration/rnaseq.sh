
ultraplex -i ERR618765.fastq.gz -b barcodes.csv -d demux -o ERR618765 -t 8 -sb -u -inm
ultraplex -i ERR618766.fastq.gz -b barcodes.csv -d demux -o ERR618766 -t 8 -sb -u -inm
ultraplex -i ERR618767.fastq.gz -b barcodes.csv -d demux -o ERR618767 -t 8 -sb -u -inm

DATADIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid/rnaseq
RESULTSDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid
REFDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref

## LOAD REQUIRED MODULES
ml purge
ml Nextflow
ml Singularity/3.4.2
ml Graphviz
ml CAMP_proxy

export NXF_SINGULARITY_CACHEDIR=/camp/lab/luscomben/home/shared/singularity

## UPDATE PIPELINE
nextflow pull nf-core/rnaseq -r 3.1

## RUN PIPELINE
nextflow run nf-core/rnaseq -r 3.1 \
--input hek293_rnaseq_demux_samplesheet.csv \
--outdir results_demux \
--fasta $RESULTSDIR/GRCh38.primary_assembly.genome.fa.gz \
--gtf $REFDIR/gencode.v33.annotation.gtf.gz \
--gencode \
--skip_trimming \
--aligner star_salmon \
--pseudo_aligner salmon \
-profile crick