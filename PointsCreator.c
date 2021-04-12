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

void readFromMPIFile(int myrank, int numranks, int numPoints) {
    MPI_File mpiFile;
    MPI_Status stat;
    Point * points = malloc(numPoints * sizeof(Point));
    int rc = MPI_File_open(MPI_COMM_WORLD, "points.bin", 
        MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiFile);

    rc = MPI_File_read_at(mpiFile, sizeof(Point) * myrank * numPoints, points,
        numPoints, MPI_POINT, &stat);

    MPI_Barrier(MPI_COMM_WORLD);
    rc = MPI_File_close(&mpiFile);

    for (int r = 0; r < numranks; ++r) {
        if (r == myrank) {
            printf("*For Rank %d:\n", myrank);
            for (int p = 0; p < numPoints; ++p)
                printf("{%f, %f}\n", points[p].x, points[p].y);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

// TODO : add MPI and parallel file processing for even larger files
int main(int argc, char* argv[])
{
    if (argc < 7) {
        printf("FORMAT: %s <num points per proc> <num threads> <left x> <lower y> <right x> <upper y>\n", argv[0]);
        return false;
    }
    int myrank, numranks;
    MPI_Init(&argc, &argv); // init MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);
    createPointType();

    size_t numPoints = 0;
    if (1 != sscanf(argv[1], "%zu", &numPoints))
        return false;
    int threadsCount = atoi(argv[2]);
    float leftX = atof(argv[3]);
    float lowerY = atof(argv[4]);
    float rightX = atof(argv[5]);
    float upperY = atof(argv[6]);
    printf("Number elements for rank %d is %zu with %d threads.\n", myrank, numPoints, threadsCount);
    
    Point * p_data;
    c_cudaAlloc(&p_data, numPoints);
    if (!p_data)
        printf("Rank %d failed to allocation.\n", myrank);

    size_t numBlocks = ((numPoints) + (threadsCount - 1)) / threadsCount;
    numBlocks = (numBlocks > 65535) ? 65535 : numBlocks;
    int seed = 10000 * myrank + numranks;
    c_callKernel(numBlocks, threadsCount, p_data, numPoints, seed,
        leftX, lowerY, rightX, upperY);

    writeToMPIFile(p_data, myrank, numranks, numPoints);
#ifdef DEBUG
    for (int r = 0; r < numranks; ++r) {
        if (r == myrank) {
            printf("*For Rank %d:\n", myrank);
            for (int p = 0; p < numPoints; ++p)
                printf("{%f, %f}\n", p_data[p].x, p_data[p].y);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    printf("\n");
    readFromMPIFile(myrank, numranks, numPoints);
#endif
    return true;
}
