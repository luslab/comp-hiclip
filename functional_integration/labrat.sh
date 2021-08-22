conda install -c bioconda labrat
conda install -c bioconda tbb=2020.2


cd /camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid/rnaseq/labrat/db
sbatch -N 1 --mem 32G --time=24:00:00 --wrap="
LABRAT.py \
--mode makeTFfasta \
--gff gencode.v33.primary_assembly.annotation.gff3.gz \
--genomefasta GRCh38.primary_assembly.genome.fa.gz \
--lasttwoexons \
--librarytype RNAseq"

/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid/rnaseq/labrat/ref
LABRAT.py \
--mode makeTFfasta \
--gff gencodecomprehensive.v28.gff3 \
--genomefasta GRCh38.primary_assembly.genome.fa.gz \
--lasttwoexons \
--librarytype RNAseq

mkdir results_salmon
cd /camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid/rnaseq/labrat/results_salmon
DATADIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid/rnaseq/demux/

LABRAT.py \
--mode runSalmon \
--txfasta ../ref/TFseqs.fasta \
--reads1 $DATADIR/ultraplex_ERR618766_5bc_NNNCGGANN.fastq.gz,$DATADIR/ultraplex_ERR618766_5bc_NNNTGGCNN.fastq.gz,$DATADIR/ultraplex_ERR618767_5bc_NNNCGGANN.fastq.gz,$DATADIR/ultraplex_ERR618767_5bc_NNNTGGCNN.fastq.gz,$DATADIR/ultraplex_ERR618765_5bc_NNNCGGANN.fastq.gz,$DATADIR/ultraplex_ERR618765_5bc_NNNTGGCNN.fastq.gz \
--samplename untreated_1,untreated_2,knockdown_1,knockdown_2,rescue_1,rescue_2 \
--threads 8


mkdir results_psi
cd /camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid/rnaseq/labrat/results_psi

LABRAT.py \
--mode calculatepsi \
--salmondir ../results_salmon \
--sampconds ../sampconds.tsv \
--conditionA rescue \
--conditionB knockdown \
--gff ../ref/gencodecomprehensive.v28.gff3 \
--librarytype RNAseq

