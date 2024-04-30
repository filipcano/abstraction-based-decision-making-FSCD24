#!/bin/bash

#compile
g++ -O2 -o balanced_server.o balanced_server.cc

# Setting the memory limit to 200GB
ulimit -v $((200 * 1024 * 1024))

for j in 1 2 3
do
   for i in 1 2 3 4 5 6 7 8 9 10
   do
      # Running each command with a 30-second time limit
      timeout 3600s ./balanced_server.o --compute_naive <inputs/max_response_input_$i.txt
      echo inputs/max_response_input_$i.txt
      ret_code=$?
      if [ $ret_code -eq 124 ]; then
         echo "Command timed out for input $i"
         break
      elif [ $ret_code -ne 0 ]; then
         echo "Command potentially failed due to memory limit for input $i"
         break
      fi
   done

   for i in 1 2 3 4 5 6 7 8 9 10 5000 10000 15000 20000 25000 30000 35000 40000 45000 50000 55000 60000
   do
      # Running each command with a 30-second time limit
      timeout 3600s ./balanced_server.o --compute_optimized <inputs/max_response_input_$i.txt
      ret_code=$?
      if [ $ret_code -eq 124 ]; then
         echo "Command timed out for input $i"
         break
      elif [ $ret_code -ne 0 ]; then
         echo "Command potentially failed due to memory limit for input $i"
         break
      fi
   done

   for i in 1 2 3 4 5 6 7 8 9 10 100 250 500 750 1000 1250 1500 1750 2000 2250 2500 2750 3000
   do
      # Running each command with a 30-second time limit
      timeout 3600s ./balanced_server.o --compute_mid <inputs/max_response_input_$i.txt
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