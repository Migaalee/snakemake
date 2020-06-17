#!/bin/bash

echo "Generating pooled dataset..."

mkdir -p 06_stage2_generate_pooled_dataset
cd 06_stage2_generate_pooled_dataset

# Script to generate pooled datasets, based on ENCODE pipeline
# See https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#
# Adapted by hywel 20190416

# ========================================
# Inputs
# ========================================

REP1_TA_FILE="../$1"
REP2_TA_FILE="../$2"

REP1_PR1_TA_FILE="../$3"
REP1_PR2_TA_FILE="../$4"

REP2_PR1_TA_FILE="../$5"
REP2_PR2_TA_FILE="../$6"

# ========================
# Create pooled datasets
# =======================

POOLED_TA_FILE=rep0.filt.nodup.srt.SE.tagAlign.gz

zcat ${REP1_TA_FILE} ${REP2_TA_FILE} | gzip -nc > ${POOLED_TA_FILE}

# ========================
# Create pooled pseudoreplicates
# =======================

PPR1_TA_FILE=rep0.filt.nodup.srt.pr1.SE.tagAlign.gz
PPR2_TA_FILE=rep0.filt.nodup.srt.pr2.SE.tagAlign.gz

zcat ${REP1_PR1_TA_FILE} ${REP2_PR1_TA_FILE} | gzip -nc > ${PPR1_TA_FILE}
zcat ${REP1_PR2_TA_FILE} ${REP2_PR2_TA_FILE} | gzip -nc > ${PPR2_TA_FILE}
