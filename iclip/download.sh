#!/bin/bash

# Download TDP-43 iCLIP data for
# A. M. Chakrabarti
# 14th December 2022

ml Aspera-Connect/3.6.1

#!/usr/bin/env bash
ascp -QT -l 300m -P33001 -i /camp/apps/eb/software/Aspera-Connect/3.6.1/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR153/000/ERR1530360/ERR1530360.fastq.gz . && mv ERR1530360.fastq.gz ERR1530360_iCLIP_TARDBP_293Flp_TARDBP_LC_flag_GFP_IP_Hs_NNNTTTCNN_20151023_JU1_10.fastq.gz
ascp -QT -l 300m -P33001 -i /camp/apps/eb/software/Aspera-Connect/3.6.1/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR153/001/ERR1530361/ERR1530361.fastq.gz . && mv ERR1530361.fastq.gz ERR1530361_iCLIP_TARDBP_293Flp_TARDBP_LC_flag_GFP_IP_Hs_NNNACTCNN_20151023_JU1_11.fastq.gz

conda activate comp-hiclip-dev
umi_tools extract -p NNNXXXXNN -I ERR1530360_iCLIP_TARDBP_293Flp_TARDBP_LC_flag_GFP_IP_Hs_NNNTTTCNN_20151023_JU1_10.fastq.gz -S ERR1530360_iCLIP_TARDBP_293Flp_TARDBP_LC_flag_GFP_IP_Hs_NNNTTTCNN_20151023_JU1_10.umi.fastq.gz 
umi_tools extract -p NNNXXXXNN -I ERR1530361_iCLIP_TARDBP_293Flp_TARDBP_LC_flag_GFP_IP_Hs_NNNACTCNN_20151023_JU1_11.fastq.gz -S ERR1530361_iCLIP_TARDBP_293Flp_TARDBP_LC_flag_GFP_IP_Hs_NNNACTCNN_20151023_JU1_11.umi.fastq.gz 

cutadapt -j 8 -u 4 -o ERR1530360_iCLIP_TARDBP_293Flp_TARDBP_LC_flag_GFP_IP_Hs_NNNTTTCNN_20151023_JU1_10.rbc.fastq.gz ERR1530360_iCLIP_TARDBP_293Flp_TARDBP_LC_flag_GFP_IP_Hs_NNNTTTCNN_20151023_JU1_10.umi.fastq.gz 
cutadapt -j 8 -u 4 -o ERR1530361_iCLIP_TARDBP_293Flp_TARDBP_LC_flag_GFP_IP_Hs_NNNACTCNN_20151023_JU1_11.rbc.fastq.gz ERR1530361_iCLIP_TARDBP_293Flp_TARDBP_LC_flag_GFP_IP_Hs_NNNACTCNN_20151023_JU1_11.umi.fastq.gz 