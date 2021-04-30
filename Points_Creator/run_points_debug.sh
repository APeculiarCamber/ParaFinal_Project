module load xl_r spectrum-mpi cuda/10.2
mpixlc -g PointsCreator.c -c -o pc-mpi.o -D DEBUG
nvcc -g -G -arch=sm_70 PointsCreator.cu -c -o pc-cuda.o -D DEBUG
mpicc -g pc-mpi.o pc-cuda.o -o pc-exe -L/usr/local/cuda-10.2/lib64/ -lcudadevrt -lcudart -lstdc++ -D DEBUG
sbatch -N $1 --ntasks-per-node=$2 --partition=dcs --gres=gpu:4 -t 10 ./slurm_pc.sh $3 $4 $5 $6 $7 $8 $9
# run_points.sh <number nodes> <number ranks> <number of points to make PER RANK> <thread count PER RANK> <left> <lower> <right> <upper>
# left < right, lower < upper. The thread count shouldn't exceeed 32. The number of ranks shouldn't exceed 32.
# Keep in mind that number of points is PER RANK, so you will create (number of ranks) * (number of points) TOTAL POINTS.
# For data collection use 1024 for K and 1024 * 1024 for M.  
