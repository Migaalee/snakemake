#!/bin/bash

echo "Calling peaks..."

mkdir -p 08_stage3_call_peaks
cd 08_stage3_call_peaks

# MACS2 needs a special conda environment as it depends on python 2

source ~/miniconda3/bin/activate macs2

# Script to call peaks, based on ENCODE pipeline
# See https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#


# ========================================
# Inputs
# ========================================

CHIP_TA_PREFIX="$1"
CHIP_TA_PATH="../$2"
CONTROL_TA_PATH="../$3"
FRAGLEN_PATH="../$4"

# ========================================
# Outputs
# ========================================

PEAK_OUTPUT_DIR=`pwd`
FRAGLEN=`cut -f3 "$FRAGLEN_PATH" | sed 's/[^0-9]//g'`

# ========================================
# Testing
# ========================================

echo "$CHIP_TA_PREFIX"
ls "$CHIP_TA_PATH" "$CONTROL_TA_PATH" "$FRAGLEN_PATH"
echo "$PEAK_OUTPUT_DIR"
echo "$FRAGLEN"
echo `which macs2`

# ========================================
# MACS2
# ========================================

GENOMESIZE="dm" 
NPEAKS=500000 # capping number of peaks called from MACS2

# ===========================================
# Generate narrow peaks and preliminary signal tracks
# ============================================

macs2 callpeak -t ${CHIP_TA_PATH} -c ${CONTROL_TA_PATH} -f BED -n ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX} -g ${GENOMESIZE} -p 1e-2 --nomodel --shift 0 --extsize ${FRAGLEN} --keep-dup all -B --SPMR

# Sort by Col8 in descending order and replace long peak names in Column 4 with Peak_<peakRank>
sort -k 8gr,8gr ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | head -n ${NPEAKS}  | gzip -nc > ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.narrowPeak.gz

# remove additional files
rm -f ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_peaks.xls ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_peaks.narrowPeak ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_summits.bed
