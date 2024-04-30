#!/bin/bash

#compile
g++ -O2 -o max_response.o max_response.cc

# Setting the memory limit to 200GB
ulimit -v $((200 * 1024 * 1024))

for j in 1 2 3
do
   for i in 1 2 3 4 5 6 7 8 9 10 11 12
   do
      # Running each command with a 30-second time limit
      timeout 3600s ./max_response.o --compute_naive <inputs/max_response_input_$i.txt
      ret_code=$?
      if [ $ret_code -eq 124 ]; then
         echo "Command timed out for input $i"
         break
      elif [ $ret_code -ne 0 ]; then
         echo "Command potentially failed due to memory limit for input $i"
         break
      fi
   done

   for i in 1 2 3 4 5 6 7 8 9 10 25 50 75 100 125 150 175 200 225 250 275 300 325 350 375 400 425
   do
      # Running each command with a 30-second time limit
      timeout 3600s ./max_response.o --compute_optimized <inputs/max_response_input_$i.txt
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