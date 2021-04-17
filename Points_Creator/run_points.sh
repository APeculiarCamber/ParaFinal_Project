./pc-make.sh
sbatch -N $1 --ntasks-per-node=$2 --gres=gpu:4 -t 10 ./slurm_pc.sh $3 $4 $5 $6 $7 $8 $9
