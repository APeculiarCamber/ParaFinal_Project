module load spectrum-mpi xl_r
mpicc Monotone_chain.c -o mono_chain
sbatch -N $1 --ntasks-per-node=$2 --gres=gpu:1 -t 10 ./slurm_chain.sh
