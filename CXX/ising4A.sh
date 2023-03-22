#!/bin/sh
# Script for running Project 4A


cd build/
make
J=-1.0
N=16
part="4A"

for B in {1..10..2};
do
  filename=../Data/Part${part}/ising_${N}_${B}_${part}.dat
  echo $filename
  ./ISING $N $J $B $filename &
done
wait

cd ..
