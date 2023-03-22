#!/bin/sh
# Script for running Project 4B


cd build/
make
J=1.0
B=0.0
part="4B"

for exp in $(seq 3 5);
do
  filename=../Data/Part${part}/ising_$((2**exp))_${part}.dat
  echo $filename
  ./ISING $((2**exp)) $J $B $filename &
done
wait

cd ..
