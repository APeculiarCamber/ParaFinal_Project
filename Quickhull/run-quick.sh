module load spectrum-mpi xl_r
mpicc Quickhull.c -o quickhull
sbatch -N $1 --ntasks-per-node=$2 --gres=gpu:1 -t 10 ./quickhull.sh
