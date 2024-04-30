#!/bin/bash

#compile
g++ -O2 -o max_response.o max_response.cc

# Setting the memory limit to 20GB
ulimit -v $((20 * 1024 * 1024))

for j in 5 10 15 20 25
do
   ./max_response.o --compute_optimized --simulate --n_simulations=100 <inputs/max_response_input_T100_X6_pa03_pb075_B$j.txt >collected_data/max_response_output_T100_X6_pa03_pb075_B5$j.txt
   echo "finished with $j"

done
