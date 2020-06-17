#!/bin/bash

echo -e "Full pipeline started in `realpath .` @ `date`\n"

setup.sh
stage1.sh
stage2.sh
stage3.sh

echo -e "\nFull pipeline finished @ `date`"
