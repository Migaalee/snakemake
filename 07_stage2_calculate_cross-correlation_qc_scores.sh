#!/bin/bash

echo "Calculating cross-correlation qc scores..."

mkdir -p 07_stage2_calculate_cross-correlation_qc_scores
cd 07_stage2_calculate_cross-correlation_qc_scores

# phantompeakqualtools needs a special conda environment as it depends on readline 6.2

source ~/miniconda3/bin/activate phantompeakqualtools

# Script to calculate cross-correlation QC scores, based on ENCODE pipeline
# See https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#
# Adapted by hywel 20190418

# ========================================
# Inputs
# ========================================

FINAL_TA_FILE="../$1"

FINAL_TA_FILE_NAME=`basename "$FINAL_TA_FILE"`
echo "FINAL_TA_FILE_NAME: $FINAL_TA_FILE_NAME"

FINAL_TA_PREFIX=${FINAL_TA_FILE_NAME%.tagAlign.gz}
echo "FINAL_TA_PREFIX: $FINAL_TA_PREFIX"

# ========================================
# General constants
# ========================================

NTHREADS=1

# =================================
# Subsample tagAlign file
# ================================

NREADS=15000
SUBSAMPLED_TA_FILE="${FINAL_TA_PREFIX}.nodup.filt.sample.$((NREADS / 1000)).SE.tagAlign.gz"

zcat ${FINAL_TA_FILE} | grep -v 'chrM' | shuf -n ${NREADS} --random-source=${FINAL_TA_FILE} | gzip -nc > ${SUBSAMPLED_TA_FILE}

# ========================================
# Calculate cross-correlation QC scores
# ========================================

CC_SCORES_FILE="${FINAL_TA_PREFIX}.cc.qc"
CC_PLOT_FILE="${FINAL_TA_PREFIX}.cc.plot.pdf"

# CC_SCORE FILE format
# Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag

run_spp.R -c=${SUBSAMPLED_TA_FILE} -p=${NTHREADS} -filtchr=chrM -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE}

sed -r 's/,[^\t]+//g' ${CC_SCORES_FILE} > temp
mv temp ${CC_SCORES_FILE}

