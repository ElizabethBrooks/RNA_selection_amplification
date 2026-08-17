#!/bin/bash

# created by: Elizabeth Brooks

# script to run the sequence of scripts that identify conserved regions
# usage: bash 13_conservation.sh

# status message 
echo "Beginning analysis of conserved regions..."

# run the analysis steps
#bash 13a_conservation.sh
bash 13b_conservation.sh
bash 13c_conservation_top10.sh
bash 13d_conservation_families.sh
qsub 13e_conservation_r0.sh 1
qsub 13e_conservation_r0.sh 2
qsub 13e_conservation_r0.sh 3
qsub 13e_conservation_r0.sh 4

# status message
echo "Conservation analysis running!"
