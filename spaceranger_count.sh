#!/bin/bash
# spaceranger count for FFPE sample
# $1: output name
# $2: prefix of fastq files
# $3: image file
# $4: capture area
# $5: loupe alignment
spaceranger count --id=$1 \
                   --transcriptome=./tenx_spaceranger/refdata-gex-mm10-2020-A \
                   --probe-set=./software/spaceranger-1.3.1/probe_sets/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --fastqs=./fastq \
                   --sample=$2 \
                   --image=./image/$3 \
                   --slide=V11F09-051 \
                   --area=$4 \
                   --localcores=64 \
                   --localmem=512 \
                   --loupe-alignment=./image/$5
