#!/bin/bash

# Check that prerequisite scripts have been run
if [[ ! -L srcdata || ! -L refdata ]]; then
    echo "Please run setup.sh in the current directory before running this file"
    exit 1
fi

if [[ ! -d 01_stage1_alignment || ! -d 02_stage1_post_alignment_filtering ]]; then
    echo "Please run stage1.sh in the current directory before running this file"
    exit 1
fi

# Activate conda environment for this script
source ~/miniconda3/bin/activate chipseq

# This works after conda activation
set -ueo pipefail

#------------ Start of Stage 2 ------------#

echo -e "---- Stage 2 ----\n"
echo "started in `realpath .` @ `date`"

#- Initialising variables from previous step

CONTROL_FILT_BAM=02_stage1_post_alignment_filtering/control.filt.nodup.srt.bam
REP1_FILT_BAM=02_stage1_post_alignment_filtering/rep1.filt.nodup.srt.bam
REP2_FILT_BAM=02_stage1_post_alignment_filtering/rep2.filt.nodup.srt.bam

#- Creating gzipped tagAlign files for each dataset

echo -e "\n-- step 4: creation of tagAlign.gz files --\n"

04_stage2_create_tagalign_file.sh "$CONTROL_FILT_BAM"
CONTROL_TA_FILE=04_stage2_create_tagalign_file/control.filt.nodup.srt.SE.tagAlign.gz

04_stage2_create_tagalign_file.sh "$REP1_FILT_BAM"
REP1_TA_FILE=04_stage2_create_tagalign_file/rep1.filt.nodup.srt.SE.tagAlign.gz

04_stage2_create_tagalign_file.sh "$REP2_FILT_BAM"
REP2_TA_FILE=04_stage2_create_tagalign_file/rep2.filt.nodup.srt.SE.tagAlign.gz

#- Creating self pseudoreplicates

echo -e "\n-- step 5: creation of self pseudoreplicates --\n"

05_stage2_create_self_pseudoreplicates.sh "$CONTROL_TA_FILE"
CONTROL_PR1_TA_FILE=05_stage2_create_self_pseudoreplicates/control.filt.nodup.pr1.SE.tagAlign.gz
CONTROL_PR2_TA_FILE=05_stage2_create_self_pseudoreplicates/control.filt.nodup.pr2.SE.tagAlign.gz

05_stage2_create_self_pseudoreplicates.sh "$REP1_TA_FILE"
REP1_PR1_TA_FILE=05_stage2_create_self_pseudoreplicates/rep1.filt.nodup.pr1.SE.tagAlign.gz
REP1_PR2_TA_FILE=05_stage2_create_self_pseudoreplicates/rep1.filt.nodup.pr2.SE.tagAlign.gz

05_stage2_create_self_pseudoreplicates.sh "$REP2_TA_FILE"
REP2_PR1_TA_FILE=05_stage2_create_self_pseudoreplicates/rep2.filt.nodup.pr1.SE.tagAlign.gz
REP2_PR2_TA_FILE=05_stage2_create_self_pseudoreplicates/rep2.filt.nodup.pr2.SE.tagAlign.gz

#- Generating pooled dataset

echo -e "\n-- step 6: generation of pooled dataset --\n"

06_stage2_generate_pooled_dataset.sh "$REP1_TA_FILE"     "$REP2_TA_FILE" \
                                     "$REP1_PR1_TA_FILE" "$REP1_PR2_TA_FILE" \
                                     "$REP2_PR1_TA_FILE" "$REP2_PR2_TA_FILE"

REP0_TA_FILE=06_stage2_generate_pooled_dataset/rep0.filt.nodup.srt.SE.tagAlign.gz
PPR1_TA_FILE=06_stage2_generate_pooled_dataset/rep0.filt.nodup.srt.pr1.SE.tagAlign.gz
PPR2_TA_FILE=06_stage2_generate_pooled_dataset/rep0.filt.nodup.srt.pr2.SE.tagAlign.gz

#- Calculating cross-correlation qc scores

echo -e "\n-- step 7: calculation of cross-correlation QC scores --\n"

07_stage2_calculate_cross-correlation_qc_scores.sh "$CONTROL_TA_FILE"
07_stage2_calculate_cross-correlation_qc_scores.sh "$REP1_TA_FILE"
07_stage2_calculate_cross-correlation_qc_scores.sh "$REP2_TA_FILE"

07_stage2_calculate_cross-correlation_qc_scores.sh "$REP1_PR1_TA_FILE"
07_stage2_calculate_cross-correlation_qc_scores.sh "$REP2_PR1_TA_FILE"
07_stage2_calculate_cross-correlation_qc_scores.sh "$REP1_PR2_TA_FILE"
07_stage2_calculate_cross-correlation_qc_scores.sh "$REP2_PR2_TA_FILE"

07_stage2_calculate_cross-correlation_qc_scores.sh "$REP0_TA_FILE"
07_stage2_calculate_cross-correlation_qc_scores.sh "$PPR1_TA_FILE"
07_stage2_calculate_cross-correlation_qc_scores.sh "$PPR2_TA_FILE"

REP1_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep1.filt.nodup.srt.cc.qc
REP2_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep2.filt.nodup.srt.cc.qc
REP0_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep0.filt.nodup.srt.cc.qc

REP1_PR1_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep1.filt.nodup.pr1.SE.cc.qc
REP2_PR1_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep2.filt.nodup.pr1.SE.cc.qc
REP1_PR2_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep1.filt.nodup.pr2.SE.cc.qc
REP2_PR2_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep2.filt.nodup.pr2.SE.cc.qc

PPR1_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep0.filt.nodup.srt.pr1.SE.cc.qc
PPR2_CC_QC_FILE=07_stage2_calculate_cross-correlation_qc_scores/rep0.filt.nodup.srt.pr2.SE.cc.qc

echo -e "\nStage 2 finished @ `date`"
