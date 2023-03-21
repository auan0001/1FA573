#!/bin/sh

cd build/
make

for i in $(seq 3 5);
do
  filename=ising_$((2**i))_$(date "+%Y.%m.%d-%H.%M.%S").dat
  echo $filename
  ./ISING $((2**i)) $filename
done

cd ..
