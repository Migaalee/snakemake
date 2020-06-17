#!/bin/bash

echo "Filtering aligned files..."

mkdir -p 02_stage1_post_alignment_filtering
cd 02_stage1_post_alignment_filtering

# Script to post-process aligned bam files, based on ENCODE pipeline
# See https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#
# Adapted by hywel 20190412

# ========================================
# Inputs
# ========================================

RAW_BAM_FILE="../$1"

RAW_BAM_NAME=`basename "$RAW_BAM_FILE"`
echo "RAW_BAM_NAME: $RAW_BAM_NAME"

OFPREFIX=${RAW_BAM_NAME%.raw.srt.bam}
echo "OFPREFIX: $OFPREFIX"

# ========================================
# General constants
# ========================================

MAPQ_THRESH=30

# =============================
# Remove  unmapped, mate unmapped
# not primary alignment, reads failing platform
# Remove low MAPQ reads
# ============================

FILT_BAM_PREFIX="${OFPREFIX}.filt.srt"
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"

samtools view -F 1804 -q ${MAPQ_THRESH} -b ${RAW_BAM_FILE} -o ${FILT_BAM_FILE}

# ======================
# Mark duplicates
# ======================

DUP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.dupmark.bam"
DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc" # QC file

picard MarkDuplicates INPUT=${FILT_BAM_FILE} OUTPUT=${DUP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false

# ============================
# Remove duplicates
# Index final position sorted BAM
# ============================

FINAL_BAM_PREFIX="${OFPREFIX}.filt.nodup.srt"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai" # To be stored

samtools view -F 1804 -b ${DUP_FILT_BAM_FILE} -o ${FINAL_BAM_FILE}

# Index Final BAM file
samtools index ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX_FILE}

# Clean up
rm ${FILT_BAM_FILE} ${DUP_FILT_BAM_FILE}
