#!/bin/bash

echo "Creating tagAlign file..."

mkdir -p 04_stage2_create_tagalign_file
cd 04_stage2_create_tagalign_file

# Script to create a tagAlign file, based on ENCODE pipeline
# See https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#
# Adapted by hywel 20190412

# ========================================
# Inputs
# ========================================

FINAL_BAM_FILE="../$1"

FINAL_BAM_FILE_NAME=`basename "$FINAL_BAM_FILE"`
echo "FINAL_BAM_FILE_NAME: $FINAL_BAM_FILE_NAME"

FINAL_BAM_PREFIX=${FINAL_BAM_FILE_NAME%.bam}
echo "FINAL_BAM_PREFIX: $FINAL_BAM_PREFIX"

# ===================
# Create tagAlign file
# ===================

# Create SE tagAlign file
FINAL_TA_FILE="${FINAL_BAM_PREFIX}.SE.tagAlign.gz"

bedtools bamtobed -i ${FINAL_BAM_FILE} | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -nc > ${FINAL_TA_FILE}
