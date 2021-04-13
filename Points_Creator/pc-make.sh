module load xl_r spectrum-mpi cuda/10.2
mpixlc -g PointsCreator.c -c -o pc-mpi.o -D DEBUG
nvcc -g -G -arch=sm_70 PointsCreator.cu -c -o pc-cuda.o -D DEBUG
mpicc -g pc-mpi.o pc-cuda.o -o pc-exe -L/usr/local/cuda-10.2/lib64/ -lcudadevrt -lcudart -lstdc++
