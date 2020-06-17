#!/bin/bash

# Check that prerequisite scripts have been run
if [[ ! -L srcdata || ! -L refdata ]]; then
    echo "Please run setup.sh in the current directory before running this file"
    exit 1
fi

# Activate conda environment for this script
source ~/miniconda3/bin/activate chipseq

# This works after conda activation
set -ueo pipefail

#------------ Start of Stage 1 ------------#

echo -e "---- Stage 1 ----\n"
echo "started in `realpath .` @ `date`"

#- Align the control and two reps to the INDEX. Each one outputs a BAM file

INDEX="refdata/index/dm6"

echo -e "\n-- step 1: alignment --\n"

01_stage1_alignment.sh "$INDEX" srcdata/control.fastq
CONTROL_SRT_BAM=01_stage1_alignment/control.raw.srt.bam

01_stage1_alignment.sh "$INDEX" srcdata/rep1.fastq
REP1_SRT_BAM=01_stage1_alignment/rep1.raw.srt.bam

01_stage1_alignment.sh "$INDEX" srcdata/rep2.fastq
REP2_SRT_BAM=01_stage1_alignment/rep2.raw.srt.bam

#- Post alignment filtering. Outputs a filtered, deduplicated version of each BAM file

echo -e "\n-- step 2: post alignment filtering --\n"

02_stage1_post_alignment_filtering.sh "$CONTROL_SRT_BAM"
CONTROL_FILT_BAM=02_stage1_post_alignment_filtering/control.filt.nodup.srt.bam

02_stage1_post_alignment_filtering.sh "$REP1_SRT_BAM"
REP1_FILT_BAM=02_stage1_post_alignment_filtering/rep1.filt.nodup.srt.bam

02_stage1_post_alignment_filtering.sh "$REP2_SRT_BAM"
REP2_FILT_BAM=02_stage1_post_alignment_filtering/rep2.filt.nodup.srt.bam

#- Generating bedgraph files for visualisation

echo -e "\n-- step 3: generation of bedgraph files --\n"

CHROMOSOME_LENGTHS_PATH=refdata/index/dm6.genome

03_stage1_generate_bedgraph.sh "$CHROMOSOME_LENGTHS_PATH" "$CONTROL_FILT_BAM"
CONTROL_BEDGRAPH_PLUS_STRAND=03_stage1_generate_bedgraph/control.filt.nodup.srt.plus_strand.bedgraph
CONTROL_BEDGRAPH_MINUS_STRAND=03_stage1_generate_bedgraph/control.filt.nodup.srt.minus_strand.bedgraph

03_stage1_generate_bedgraph.sh "$CHROMOSOME_LENGTHS_PATH" "$REP1_FILT_BAM"
REP1_BEDGRAPH_PLUS_STRAND=03_stage1_generate_bedgraph/rep1.filt.nodup.srt.plus_strand.bedgraph
REP1_BEDGRAPH_MINUS_STRAND=03_stage1_generate_bedgraph/rep1.filt.nodup.srt.minus_strand.bedgraph

03_stage1_generate_bedgraph.sh "$CHROMOSOME_LENGTHS_PATH" "$REP2_FILT_BAM"
REP2_BEDGRAPH_PLUS_STRAND=03_stage1_generate_bedgraph/rep2.filt.nodup.srt.plus_strand.bedgraph
REP2_BEDGRAPH_MINUS_STRAND=03_stage1_generate_bedgraph/rep2.filt.nodup.srt.minus_strand.bedgraph

echo -e "\nStage 1 finished @ `date`"
