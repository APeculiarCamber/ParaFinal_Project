rm slurm-*
module load spectrum-mpi xl_r
mpicc Quickhull.c -o quickhull

./run-quick-bulk.sh 1 1 16777216 24
./run-quick-bulk.sh 1 2 16777216 24
./run-quick-bulk.sh 1 4 16777216 24
./run-quick-bulk.sh 1 8 16777216 24
./run-quick-bulk.sh 1 16 16777216 24
./run-quick-bulk.sh 1 32 16777216 24
./run-quick-bulk.sh 2 32 16777216 24

./run-quick-bulk.sh 1 1 33554432 25
./run-quick-bulk.sh 1 2 33554432 25
./run-quick-bulk.sh 1 4 33554432 25
./run-quick-bulk.sh 1 8 33554432 25
./run-quick-bulk.sh 1 16 33554432 25
./run-quick-bulk.sh 1 32 33554432 25
./run-quick-bulk.sh 2 32 33554432 25

./run-quick-bulk.sh 1 1 67108864 26
./run-quick-bulk.sh 1 2 67108864 26
./run-quick-bulk.sh 1 4 67108864 26
./run-quick-bulk.sh 1 8 67108864 26
./run-quick-bulk.sh 1 16 67108864 26
./run-quick-bulk.sh 1 32 67108864 26
./run-quick-bulk.sh 2 32 67108864 26

./run-quick-bulk.sh 1 1 134217728 27
./run-quick-bulk.sh 1 2 134217728 27
./run-quick-bulk.sh 1 4 134217728 27
./run-quick-bulk.sh 1 8 134217728 27
./run-quick-bulk.sh 1 16 134217728 27
./run-quick-bulk.sh 1 32 134217728 27
./run-quick-bulk.sh 2 32 134217728 27

./run-quick-bulk.sh 1 1 268435456 28
./run-quick-bulk.sh 1 2 268435456 28
./run-quick-bulk.sh 1 4 268435456 28
./run-quick-bulk.sh 1 8 268435456 28
./run-quick-bulk.sh 1 16 268435456 28
./run-quick-bulk.sh 1 32 268435456 28
./run-quick-bulk.sh 2 32 268435456 28

./run-quick-bulk.sh 1 1 536870912 29
./run-quick-bulk.sh 1 2 536870912 29
./run-quick-bulk.sh 1 4 536870912 29
./run-quick-bulk.sh 1 8 536870912 29
./run-quick-bulk.sh 1 16 536870912 29
./run-quick-bulk.sh 1 32 536870912 29
./run-quick-bulk.sh 2 32 536870912 29

./run-quick-bulk.sh 1 1 1073741818 30
./run-quick-bulk.sh 1 2 1073741824 30
./run-quick-bulk.sh 1 4 1073741824 30
./run-quick-bulk.sh 1 8 1073741824 30
./run-quick-bulk.sh 1 16 1073741824 30
./run-quick-bulk.sh 1 32 1073741824 30
./run-quick-bulk.sh 2 32 1073741824 30
