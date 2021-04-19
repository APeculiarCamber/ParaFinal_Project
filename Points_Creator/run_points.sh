./pc-make.sh
sbatch -N $1 --ntasks-per-node=$2 --gres=gpu:4 -t 10 ./slurm_pc.sh $3 $4 $5 $6 $7 $8 $9
# run_points.sh <number nodes> <number ranks> <number of points to make PER RANK> <thread count PER RANK> <left> <lower> <right> <upper>
# left < right, lower < upper. The thread count shouldn't exceeed 32. The number of ranks shouldn't exceed 32.
# Keep in mind that number of points is PER RANK, so you will create (number of ranks) * (number of points) TOTAL POINTS.
# For data collection use 1024 for K and 1024 * 1024 for M.  