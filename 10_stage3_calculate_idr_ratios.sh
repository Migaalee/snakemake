#!/bin/bash

echo "Calculating IDR ratios..."

mkdir -p 10_stage3_calculate_idr_ratios
cd 10_stage3_calculate_idr_ratios

# Script to calculate IDR ratios, see  https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#
# Written by hywel 20190419

# ========================================
# Inputs
# ========================================

NT_PATH="../$1"
NP_PATH="../$2"

N1_PATH="../$3"
N2_PATH="../$4"

# ========================================
# Outputs
# ========================================

IDR_RATIOS_PATH=idr_ratios.tsv

# ========================================
# Calculate IDR ratios
# ========================================

echo -e "NP\tNF\tN1\tN2\tRescue ratio\tSelf-consistency ratio" > ${IDR_RATIOS_PATH}

cat ${NP_PATH} ${NT_PATH} ${N1_PATH} ${N2_PATH} | tr '\n' '\t' | awk -F'\t' '{rrmax=$1>$2?$1:$2; rrmin=$2>$1?$1:$2; rrres=rrmin>0?(rrmax/rrmin):"NaN"; scrmax=$3>$4?$3:$4; scrmin=$4>$3?$3:$4; scrres=scrmin>0?(scrmax/scrmin):"NaN"; print $1"\t"$2"\t"$3"\t"$4"\t"rrres"\t"scrres}' >> ${IDR_RATIOS_PATH}


#echo "NT:"
#echo `cat ${NT_PATH}`
#
#echo "N1:"
#echo `cat ${N1_PATH}`
#
#echo "N2:"
#echo `cat ${N2_PATH}`
#
#echo "NP:"
#echo `cat ${NP_PATH}`

#calculate_ratio_cmd () {
    
#    cat "$1" "$2" | tr '\n' '\t' | awk -F'\t' '{max=$1>$2?$1:$2; min=$2>$1?$1:$2; res=min>0?(max/min):"NaN"; print res}'

#}

#echo "Rescue ratio: $(calculate_ratio_cmd ${NT_PATH} ${NP_PATH})"

#echo "Self-consistency ratio: $(calculate_ratio_cmd ${N1_PATH} ${N2_PATH})"
