module load spectrum-mpi xl_r
mpicc MC_Chan_Algo.c -o chan_algo -D DEBUG
sbatch -N $1 --ntasks-per-node=$2 --partition=dcs --gres=gpu:4 -t 30 ./slurm_chan.sh $3 $4
