#!/bin/bash

g++ -o ising p3_auan.cxx -march=native -O3 && ./ising && python Plots/plot.py
# awk -f tp.awk run.txt > runt.txt
# watch -n 1 cat run.txt
