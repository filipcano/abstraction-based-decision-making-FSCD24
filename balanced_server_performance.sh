#!/bin/bash

#compile
g++ -O2 -o balanced_server.o balanced_server.cc

# Setting the memory limit to 20GB
ulimit -v $((20 * 1024 * 1024))

for j in 5
do
   ./balanced_server.o --compute_optimized --simulate --n_simulations=100 <inputs/balanced_server_input_T1000_X20_pa051_pb049_B$j.txt >collected_data/output_balanced_server_T1000_X20_pa51_pb049_B$j.txt
   echo "finished with $j"

done