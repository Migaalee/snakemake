#!/bin/bash

echo "Generating bedgraph files..."

mkdir -p 03_stage1_generate_bedgraph
cd 03_stage1_generate_bedgraph

# Script to create genome coverage plots in bedgraph format, based on ENCODE pipeline
# Written by hywel 20190417

# ========================================
# Inputs
# ========================================

# TODO - unused so remove this
CHROMOSOME_LENGTHS_PATH="../$1"
INPUT_BAM_PATH="../$2"

INPUT_BAM_NAME=`basename "$INPUT_BAM_PATH"`
INPUT_BAM_PREFIX=${INPUT_BAM_NAME%.bam}

# ========================================
# bedtools genomecov
# ========================================

bedtools genomecov -ibam ${INPUT_BAM_PATH} -bg -strand + > ${INPUT_BAM_PREFIX}".plus_strand.bedgraph"

bedtools genomecov -ibam ${INPUT_BAM_PATH} -bg -strand - > ${INPUT_BAM_PREFIX}".minus_strand.bedgraph"
