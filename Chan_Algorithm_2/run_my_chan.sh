module load spectrum-mpi xl_r
mpicc My_Chan_Algo.c -o chan_algo
sbatch -N $1 --ntasks-per-node=$2 --gres=gpu:4 -t 30 ./slurm_chan.sh $3 $4