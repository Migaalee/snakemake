#!/bin/bash

echo "Creating self pseudoreplicates..."

mkdir -p 05_stage2_create_self_pseudoreplicates
cd 05_stage2_create_self_pseudoreplicates

# Script to create self-pseudoreplicates, based on ENCODE pipeline
# See https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#
# Adapted by hywel 20190412

# ========================================
# Inputs
# ========================================

FINAL_TA_FILE="../$1"

FINAL_TA_FILE_NAME=`basename "$FINAL_TA_FILE"`
echo "FINAL_TA_FILE_NAME: $FINAL_TA_FILE_NAME"

OFPREFIX=${FINAL_TA_FILE_NAME%.filt.nodup.srt.SE.tagAlign.gz}
echo "OFPREFIX: $OFPREFIX"

# ========================
# Create pseudoReplicates
# =======================

PR_PREFIX="${OFPREFIX}.filt.nodup"
PR1_TA_FILE="${PR_PREFIX}.pr1.SE.tagAlign.gz"
PR2_TA_FILE="${PR_PREFIX}.pr2.SE.tagAlign.gz"

# Get total number of read pairs
nlines=$( zcat ${FINAL_TA_FILE} | wc -l )
nlines=$(( (nlines + 1) / 2 ))

# Shuffle and split BED file into 2 equal parts
zcat ${FINAL_TA_FILE} | shuf --random-source=${FINAL_TA_FILE}  | split -d -l ${nlines} - ${PR_PREFIX} # Will produce ${PR_PREFIX}00 and ${PR_PREFIX}01

# Convert reads into standard tagAlign file
gzip -nc "${PR_PREFIX}00" > ${PR1_TA_FILE}
rm "${PR_PREFIX}00"
gzip -nc "${PR_PREFIX}01" > ${PR2_TA_FILE}
rm "${PR_PREFIX}01"
