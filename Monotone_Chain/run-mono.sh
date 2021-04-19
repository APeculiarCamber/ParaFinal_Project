module load spectrum-mpi xl_r
mpicc Monotone_chain.c -o mono_chain
sbatch -N $1 --ntasks-per-node=$2 --gres=gpu:4 -t 10 ./slurm_chain.sh $3 $4
# ./run-mono.sh <number of nodes> <number of mpi ranks PER NODE> <TOTAL NUMBER OF POINTS> <Oversampling>
# Mpi ranks should not exceed 32 per NODE, so when running 64 ranks, use ./run-mono.sh 2 32 ... and NOT ./run-mono.sh 1 64 ...
# In contrast to Points Creator, the number of points here is TOTAL.
# For monotone chain's sample sorting, it splits the data into bins by randomly selecting pivots; <oversampling> specifies the amount 
# of pivots per rank that should be used. Oversampling allows for more even distribution of bins. 
# A value of like ~256 might be good, or you could linearly vary it with number of points per rank.
 