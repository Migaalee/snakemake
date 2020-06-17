#!/bin/bash

# Check that prerequisite scripts have been run
if [[ ! -L srcdata || ! -L refdata ]]; then
    echo "Please run setup.sh in the current directory before running this file"
    exit 1
fi

if [[ ! -d 01_stage1_alignment || ! -d 02_stage1_post_alignment_filtering || ! -d 03_stage1_generate_bedgraph ]]; then
    echo "Please run stage1.sh in the current directory before running this file"
    exit 1
fi

if [[ ! -d 04_stage2_create_tagalign_file || ! -d 05_stage2_create_self_pseudoreplicates || ! -d 06_stage2_generate_pooled_dataset || ! -d 07_stage2_calculate_cross-correlation_qc_scores ]]; then
    echo "Please run stage2.sh in the current directory before running this file"
    exit 1
fi

# Activate conda environment for this script
source ~/miniconda3/bin/activate chipseq

# This works after conda activation
set -ueo pipefail

#------------ Start of Stage 3 ------------#

echo -e "---- Stage 3 ----\n"
echo "started in `realpath .` @ `date`"

#- Initialising variables from previous step

CONTROL_TA_FILE=04_stage2_create_tagalign_file/control.filt.nodup.srt.SE.tagAlign.gz
REP1_TA_FILE=04_stage2_create_tagalign_file/rep1.filt.nodup.srt.SE.tagAlign.gz
REP2_TA_FILE=04_stage2_create_tagalign_file/rep2.filt.nodup.srt.SE.tagAlign.gz
REP0_TA_FILE=06_stage2_generate_pooled_dataset/rep0.filt.nodup.srt.SE.tagAlign.gz

REP1_PR1_TA_FILE=05_stage2_create_self_pseudoreplicates/rep1.filt.nodup.pr1.SE.tagAlign.gz
REP1_PR2_TA_FILE=05_stage2_create_self_pseudoreplicates/rep1.filt.nodup.pr2.SE.tagAlign.gz
REP2_PR1_TA_FILE=05_stage2_create_self_pseudoreplicates/rep2.filt.nodup.pr1.SE.tagAlign.gz
REP2_PR2_TA_FILE=05_stage2_create_self_pseudoreplicates/rep2.filt.nodup.pr2.SE.tagAlign.gz
PPR1_TA_FILE=06_stage2_generate_pooled_dataset/rep0.filt.nodup.srt.pr1.SE.tagAlign.gz
PPR2_TA_FILE=06_stage2_generate_pooled_dataset/rep0.filt.nodup.srt.pr2.SE.tagAlign.gz

REP1_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep1.filt.nodup.srt.SE.cc.qc
REP2_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep2.filt.nodup.srt.SE.cc.qc
REP0_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep0.filt.nodup.srt.SE.cc.qc

REP1_PR1_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep1.filt.nodup.pr1.SE.cc.qc
REP2_PR1_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep2.filt.nodup.pr1.SE.cc.qc
REP1_PR2_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep1.filt.nodup.pr2.SE.cc.qc
REP2_PR2_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep2.filt.nodup.pr2.SE.cc.qc
PPR1_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep0.filt.nodup.srt.pr1.SE.cc.qc
PPR2_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep0.filt.nodup.srt.pr2.SE.cc.qc

#- Calling peaks using MACS2

echo -e "\n-- step 8: peak calling with MACS2 --\n"

# ChIP replicates vs control

08_stage3_call_peaks_macs2.sh "rep1" "$REP1_TA_FILE" "$CONTROL_TA_FILE" "$REP1_CC_QC_FILE"
REP1_PEAK_FILE=08_stage3_call_peaks/rep1.narrowPeak.gz

08_stage3_call_peaks_macs2.sh "rep2" "$REP2_TA_FILE" "$CONTROL_TA_FILE" "$REP2_CC_QC_FILE"
REP2_PEAK_FILE=08_stage3_call_peaks/rep2.narrowPeak.gz

# Pooled ChIP vs control

08_stage3_call_peaks_macs2.sh "rep0" "$REP0_TA_FILE" "$CONTROL_TA_FILE" "$REP0_CC_QC_FILE"
POOLED_PEAK_FILE=08_stage3_call_peaks/rep0.narrowPeak.gz

# Pseudoreplicate 1 vs control

08_stage3_call_peaks_macs2.sh "rep1_pr1" "$REP1_PR1_TA_FILE" "$CONTROL_TA_FILE" "$REP1_PR1_CC_QC_FILE"
REP1_PR1_PEAK_FILE=08_stage3_call_peaks/rep1_pr1.narrowPeak.gz

08_stage3_call_peaks_macs2.sh "rep2_pr1" "$REP2_PR1_TA_FILE" "$CONTROL_TA_FILE" "$REP2_PR1_CC_QC_FILE"
REP2_PR1_PEAK_FILE=08_stage3_call_peaks/rep2_pr1.narrowPeak.gz

# Pseudoreplicate 2 vs control

08_stage3_call_peaks_macs2.sh "rep1_pr2" "$REP1_PR2_TA_FILE" "$CONTROL_TA_FILE" "$REP1_PR2_CC_QC_FILE"
REP1_PR2_PEAK_FILE=08_stage3_call_peaks/rep1_pr2.narrowPeak.gz

08_stage3_call_peaks_macs2.sh "rep2_pr2" "$REP2_PR2_TA_FILE" "$CONTROL_TA_FILE" "$REP2_PR2_CC_QC_FILE"
REP2_PR2_PEAK_FILE=08_stage3_call_peaks/rep2_pr2.narrowPeak.gz

# Pooled pseudoreplicate 1 vs control

08_stage3_call_peaks_macs2.sh "ppr1" "$PPR1_TA_FILE" "$CONTROL_TA_FILE" "$PPR1_CC_QC_FILE"
PPR1_PEAK_FILE=08_stage3_call_peaks/ppr1.narrowPeak.gz

# Pooled pseudoreplicate 2 vs control

08_stage3_call_peaks_macs2.sh "ppr2" "$PPR2_TA_FILE" "$CONTROL_TA_FILE" "$PPR2_CC_QC_FILE"
PPR2_PEAK_FILE=08_stage3_call_peaks/ppr2.narrowPeak.gz

#- IDR analysis

echo -e "\n-- step 9: IDR analysis --\n"

09_stage3_perform_idr_analysis.sh "rep1_vs_rep2" "$REP1_PEAK_FILE" "$REP2_PEAK_FILE" "$POOLED_PEAK_FILE"
NT_PATH=09_stage3_perform_idr_analysis/rep1_vs_rep2.IDR0.05.npeaks.txt

09_stage3_perform_idr_analysis.sh "rep1_pr1_vs_rep1_pr2" "$REP1_PR1_PEAK_FILE" "$REP1_PR2_PEAK_FILE" "$REP1_PEAK_FILE"
REP1_PR_IDR_PEAK_FILE=09_stage3_perform_idr_analysis/rep1_pr1_vs_rep1_pr2.IDR0.05.narrowPeak.gz
N1_PATH=09_stage3_perform_idr_analysis/rep1_pr1_vs_rep1_pr2.IDR0.05.npeaks.txt

09_stage3_perform_idr_analysis.sh "rep2_pr1_vs_rep2_pr2" "$REP2_PR1_PEAK_FILE" "$REP2_PR2_PEAK_FILE" "$REP2_PEAK_FILE"
REP2_PR_IDR_PEAK_FILE=09_stage3_perform_idr_analysis/rep2_pr1_vs_rep2_pr2.IDR0.05.narrowPeak.gz
N2_PATH=09_stage3_perform_idr_analysis/rep2_pr1_vs_rep2_pr2.IDR0.05.npeaks.txt

09_stage3_perform_idr_analysis.sh "ppr1_vs_ppr2" "$PPR1_PEAK_FILE" "$PPR2_PEAK_FILE" "$POOLED_PEAK_FILE"
NP_PATH=09_stage3_perform_idr_analysis/ppr1_vs_ppr2.IDR0.05.npeaks.txt

#- Calculating IDR ratios

echo -e "\n-- step 10: calculation of IDR ratios --\n"

10_stage3_calculate_idr_ratios.sh "$NT_PATH" "$NP_PATH" "$N1_PATH" "$N2_PATH"
IDR_RATIOS=10_stage3_calculate_idr_ratios/idr_ratios.tsv

#- Computing fraction of reads in peaks (FRIP)

echo -e "\n-- step 11: computing fraction of reads in peaks --\n"

CHRSIZEFILE=refdata/index/dm6.genome

11_stage3_compute_fraction_of_reads_in_peaks.sh "rep1" "$REP1_TA_FILE" "$REP1_PR_IDR_PEAK_FILE" "$REP1_CC_QC_FILE" "$CHRSIZEFILE"
REP1_FRIP=11_stage3_compute_fraction_of_reads_in_peaks/rep1.frip.txt

11_stage3_compute_fraction_of_reads_in_peaks.sh "rep2" "$REP2_TA_FILE" "$REP2_PR_IDR_PEAK_FILE" "$REP2_CC_QC_FILE" "$CHRSIZEFILE"
REP2_FRIP=11_stage3_compute_fraction_of_reads_in_peaks/rep2.frip.txt

echo -e "\nStage 3 finished @ `date`"
