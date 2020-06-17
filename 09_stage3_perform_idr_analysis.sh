#!/bin/bash

echo "Performing IDR analysis..."

mkdir -p 09_stage3_perform_idr_analysis
cd 09_stage3_perform_idr_analysis

# Script to perform IDR analysis, based on ENCODE pipeline
# See https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#
# Adapted by hywel 20190419

# ========================================
# Inputs
# ========================================

REP1_VS_REP2="$1"
REP1_PEAK_FILE="../$2"
REP2_PEAK_FILE="../$3"
POOLED_PEAK_FILE="../$4"

# ========================================
# Outputs
# ========================================

IDR_OUTPUT=${REP1_VS_REP2}.idr_output.txt

# =============================
# Perform IDR analysis.
# Generate a plot and IDR output with additional columns including IDR scores.
# =============================

IDR_THRESH=0.05

idr --samples ${REP1_PEAK_FILE} ${REP2_PEAK_FILE} --peak-list ${POOLED_PEAK_FILE} --input-file-type narrowPeak --output-file ${IDR_OUTPUT} --rank signal.value --soft-idr-threshold ${IDR_THRESH} --plot --use-best-multisummit-IDR

# =============================
# Get peaks passing IDR threshold of 5%
# =============================

IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}')

awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ${IDR_OUTPUT} | sort | uniq | sort -k7n,7n | gzip -nc > ${REP1_VS_REP2}.IDR0.05.narrowPeak.gz

NPEAKS_IDR=$(zcat ${REP1_VS_REP2}.IDR0.05.narrowPeak.gz | wc -l)

echo ${NPEAKS_IDR} > ${REP1_VS_REP2}.IDR0.05.npeaks.txt
