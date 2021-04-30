module load spectrum-mpi xl_r
mpicc Quickhull.c -o quickhull -D DEBUG
sbatch -N $1 --ntasks-per-node=$2 --partition=dcs --gres=gpu:1 -t 10 ./slurm_quick.sh $3 $4
