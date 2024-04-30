#!/bin/bash

#compile
g++ -O2 -o dist_change.o dist_change.cc

# Setting the memory limit to 50GB
ulimit -v $((200 * 1024 * 1024))

for j in 1 2 3 4 5
do
   for i in 1 2 3 4 5 6
   do
      # Running each command with a 30-second time limit
      timeout 3600s ./dist_change.o --compute_naive <inputs/established_clientele_input_$i.txt
      ret_code=$?
      if [ $ret_code -eq 124 ]; then
         echo "Command timed out for input $i"
         break
      elif [ $ret_code -ne 0 ]; then
         echo "Command potentially failed due to memory limit for input $i"
         break
      fi
   done

   for i in 1 2 3 4 5 6 7 8 9 10 100 250 500 750 1000 1250 1500 1750 2000 2250 2500 2750 3000 3250 3500 3750 4000 4250 4500
   do
      # Running each command with a 30-second time limit
      timeout 3600s ./dist_change.o --compute_optimized <inputs/established_clientele_input_$i.txt
      ret_code=$?
      if [ $ret_code -eq 124 ]; then
         echo "Command timed out for input $i"
         break
      elif [ $ret_code -ne 0 ]; then
         echo "Command potentially failed due to memory limit for input $i"
         break
      fi
   done

   for i in 1 2 3 4 5 6 7 8 9 10 100 250 500 750 1000 1250 1500 1750 2000 2250 2500 2750 3000 3250 3500 3750 4000 4250 4500
   do
      # Running each command with a 30-second time limit
      timeout 3600s ./dist_change.o --compute_mid <inputs/established_clientele_input_$i.txt
      ret_code=$?
      if [ $ret_code -eq 124 ]; then
         echo "Command timed out for input $i"
         break
      elif [ $ret_code -ne 0 ]; then
         echo "Command potentially failed due to memory limit for input $i"
         break
      fi
   done
done