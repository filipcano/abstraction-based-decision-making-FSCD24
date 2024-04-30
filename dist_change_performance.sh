#!/bin/bash

#compile
g++ -O2 -o dist_change.o dist_change.cc

# Setting the memory limit to 20GB
ulimit -v $((20 * 1024 * 1024))

for j in 50 150 250 350
do
   ./dist_change.o --compute_optimized --simulate --n_simulations=100 <inputs/established_clientele_input_N100_Dcenteredbino_l-5_u4_n6_k1_T500_B${j}_S10.txt >collected_data/established_clientele_output_N100_Dcenteredbino_l-5_u4_n6_k1_T500_B${j}_S10.txt
   echo "finished with $j"

done
