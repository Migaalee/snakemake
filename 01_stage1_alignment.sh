#!/bin/bash

echo "Aligning fastq files..."

mkdir -p 01_stage1_alignment
cd 01_stage1_alignment

# Script to perform alignment with bwa and generate stats, based on ENCODE pipeline
# See https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#
# Adapted by hywel 20190410

# ========================================
# Inputs
# ========================================

BWA_INDEX_NAME="../$1"
FASTQ_FILE_1="../$2"

FASTQ_FILE_1_NAME=`basename "$FASTQ_FILE_1"`
echo "FASTQ_FILE_1_NAME: $FASTQ_FILE_1_NAME"

OFPREFIX=${FASTQ_FILE_1_NAME%.fastq}
echo "OFPREFIX: $OFPREFIX"

# ========================================
# General constants
# ========================================

NTHREADS=1

# ========================================
# Map reads to create raw SAM file
# ========================================

echo "Aligning ${FASTQ_FILE_1} file to genome..."

SAI_FILE_1="${OFPREFIX}.sai"
echo "SAI_FILE_1: $SAI_FILE_1"

RAW_BAM_PREFIX="${OFPREFIX}.raw.srt"
echo "RAW_BAM_PREFIX: $RAW_BAM_PREFIX"

RAW_BAM_FILE="${RAW_BAM_PREFIX}.bam" # To be stored
echo "RAW_BAM_FILE: $RAW_BAM_FILE"

bwa aln -q 5 -l 32 -k 2 -t ${NTHREADS} ${BWA_INDEX_NAME} ${FASTQ_FILE_1} > ${SAI_FILE_1}

bwa samse ${BWA_INDEX_NAME} ${SAI_FILE_1} ${FASTQ_FILE_1} | samtools view -Su - | samtools sort -o ${RAW_BAM_FILE} -

rm ${SAI_FILE_1}
