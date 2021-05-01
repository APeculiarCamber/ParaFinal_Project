#!/bin/bash

rm slurm-*
module load spectrum-mpi xl_r
mpicc Quickhull.c -o quickhull

for i in `seq 24 30`; do
  for j in `seq 0 6`; do
    rank=$(echo "2^${j}" | bc)
    nodes=1
    if [[ $rank -gt 32 ]]
    then
      nodes=2
      rank=32
    fi
    exp=$(echo "2^${i}" | bc)
    echo $nodes $rank $exp
    ./run-quick-bulk.sh $nodes $rank $exp
  done
done