#include<unistd.h>
#include<stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

struct Point {
    float x;
    float y;
};
typedef struct Point Point;

MPI_Datatype MPI_POINT;
void createPointType() {
    int blockLengths[2] = {1,1};
    MPI_Datatype types[2] = {MPI_FLOAT, MPI_FLOAT};
    MPI_Aint offsets[2];

    offsets[0] = 0;
    offsets[1] = sizeof(float);

    MPI_Type_create_struct(2, blockLengths, offsets, types, &MPI_POINT);
    MPI_Type_commit(&MPI_POINT);
}


extern void c_cudaAlloc(Point ** points, size_t numPoints);
extern void c_callKernel(int numBlocks, int threadsCount, Point * p_data, size_t numPoints, int seed,
        float leftX, float lowerY, float rightX, float upperY);
extern void c_callCircleKernel(int numBlocks, int threadsCount, Point * p_data, 
        size_t numPoints, int seed, float radius);
extern void setCudaDeviceByRank(int myrank);
extern void freeCudaMemory(Point * p_data);

void writeToMPIFile(Point * p_data, int myrank, int numranks, int numPoints) {
    MPI_File mpiFile;
    MPI_Status stat;

    int rc = MPI_File_open(MPI_COMM_WORLD, "points.bin", 
        MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpiFile);

    rc = MPI_File_write_at(mpiFile, sizeof(Point) * myrank * numPoints, p_data,
        numPoints, MPI_POINT, &stat);

    MPI_Barrier(MPI_COMM_WORLD);
    rc = MPI_File_close(&mpiFile);
}

int comparePoints(const void* a, const void* b) {
  int yComp = (((Point*)a)->y > ((Point*)b)->y) - (((Point*)a)->y < ((Point*)b)->y);
  int xComp = (((Point*)a)->x > ((Point*)b)->x) - (((Point*)a)->x < ((Point*)b)->x);
  return xComp == 0 ? yComp : xComp;
}


void localSort(Point * points, size_t numPoints) {
    qsort (points, numPoints, sizeof(Point), comparePoints);
}

int main(int argc, char* argv[])
{
    int myrank, numranks;
    MPI_Init(&argc, &argv); // init MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);
    createPointType();
    setCudaDeviceByRank(myrank);

    if (argc != 4 && argc != 7) {
        if (myrank == 0)
            printf("FORMAT FOR SQUARE:\n\t%s <num points per proc> <num threads> <left x> <lower y> <right x> <upper y>\n", argv[0]);
            printf("FORMAT FOR CIRCLE:\n\t%s <num points per proc> <num threads> <radius>\n", argv[0]);
        return false;
    }

    // parse args
    size_t numPoints = 0;
    if (1 != sscanf(argv[1], "%zu", &numPoints))
        return false;
    int threadsCount = atoi(argv[2]);
    float leftX, lowerY, rightX, upperY, radius;
    char sampleSquare = argc == 7;
    if (sampleSquare) { 
    	leftX = atof(argv[3]);
    	lowerY = atof(argv[4]);
    	rightX = atof(argv[5]);
    	upperY = atof(argv[6]);
	} else {
		radius = atof(argv[3]);
	}

    printf("Number elements for rank %d is %zu with %d threads.\n", myrank, numPoints, threadsCount);
    
    Point * p_data;
    c_cudaAlloc(&p_data, numPoints);
    if (!p_data)
        printf("Rank %d failed to allocation.\n", myrank);

    size_t numBlocks = ((numPoints) + (threadsCount - 1)) / threadsCount;
    numBlocks = (numBlocks > 65535) ? 65535 : numBlocks;
    int seed = (14129 * myrank) + numranks;
    if (sampleSquare) {
    	c_callKernel(numBlocks, threadsCount, p_data, numPoints, seed,
        	leftX, lowerY, rightX, upperY);
	} else {
		c_callCircleKernel(numBlocks, threadsCount, p_data, numPoints, seed,
        	radius);
	}

    writeToMPIFile(p_data, myrank, numranks, numPoints);

#ifdef DEBUG
    for (int r = 0; r < numranks; ++r) {
        if (r == myrank) {
            for (int p = 0; p < numPoints; ++p)
                printf("{%f, %f}\n", p_data[p].x, p_data[p].y);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    freeCudaMemory(p_data);
    MPI_Finalize();
    return true;
}
