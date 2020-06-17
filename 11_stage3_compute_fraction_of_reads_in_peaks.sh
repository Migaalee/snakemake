#!/bin/bash

echo "Computing fraction of reads in peaks (FRIP)..."

mkdir -p 11_stage3_compute_fraction_of_reads_in_peaks
cd 11_stage3_compute_fraction_of_reads_in_peaks

# Script to compute the fraction of reads in peaks (FRIP), see  https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#
# Written by hywel 20190419

# ========================================
# Inputs
# ========================================

FRIP_PREFIX="$1"
TA_FILE="../$2"
IDR_PEAK_FILE="../$3"
CC_QC_LOG="../$4"
CHRSIZEFILE="../$5"

#echo ${FRIP_PREFIX}
#echo ${TA_FILE}
#echo ${IDR_PEAK_FILE}
#echo ${CC_QC_LOG}
#echo ${CHRSIZEFILE}

# ========================================
# Outputs
# ========================================

FRIP=${FRIP_PREFIX}.frip.txt

# ========================================
# Compute fraction of reads in peaks
# ========================================

# get estimated fragment length from cross-corr. analysis log
FRAGLEN=$(cat ${CC_QC_LOG} | awk '{print $3}')
HALF_FRAGLEN=$(( (FRAGLEN+1)/2 )) # rounding to integer

val1=$(bedtools slop -i ${TA_FILE} -g $CHRSIZEFILE -s -l -$HALF_FRAGLEN -r $HALF_FRAGLEN | \
    awk '{if ($2>=0 && $3>=0 && $2<=$3) print $0}' | \
    bedtools intersect -a stdin -b ${IDR_PEAK_FILE} -wa -u | wc -l)
val2=$(zcat $TA_FILE | wc -l)
awk 'BEGIN {print '${val1}'/'${val2}'}' > ${FRIP}
