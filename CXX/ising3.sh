#!/bin/sh
# Script for running Part3


cd build/
make
J=1.0
B=0.0
part=3

for exp in $(seq 3 5);
do
  filename=../Data/Part${part}/ising_$((2**exp))_${part}.dat
  echo -e $filename
  ./ISING $((2**exp)) $J $B $filename &
done
wait

cd ..
