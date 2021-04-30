module load spectrum-mpi xl_r
mpicc Chan_Algo.c -o chan_algo -D DEBUG
sbatch -N $1 --ntasks-per-node=$2 --partition=dcs --gres=gpu:4 -t 10 ./slurm_chan.sh $3 $4
