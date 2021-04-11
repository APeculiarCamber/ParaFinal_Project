#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdbool.h>
#include <string.h>
#include "mpi.h"
// *************************************
// TODO : do samplesort instead and send data on request, this allows proper division of memory


struct Point {
    float x;
    float y;
};

typedef struct Point Point;

// USING MERGE SORT

void readInData(char * fileName, Point * points, 
    int myrank, size_t stride, size_t numPoints) {
    MPI_File mpiFile;
    MPI_Status stat;

    int rc = MPI_File_open(MPI_COMM_WORLD, fileName, 
        MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiFile);
    rc = MPI_File_read_at(mpiFile, (myrank * stride) * 2, points,
                  numPoints * 2, MPI_FLOAT, &stat);
    rc = MPI_File_close(&mpiFile);
}

void localSort(Point * points) {
    // TODO : inplace local sort
}

void merge(Point * intermed, Point * results, int step) {
    // TODO : merge
}


Point * mergeAll() {
    int parent, child;
    Point * intermed, * merged;

    int currentStep = 2;
    int ranksLeft = numranks;
    
    while (ranksLeft > 0) { 
        parent = (id / currentStep) * currentStep;

        if (parent == id) {
            child = id + (currentStep / 2);

            intermed = realloc(merged, size * 2 * sizeof(Point));

            // TODO : very notable issue : requires size of array which fits on 1 process...
            MPI_Recv(intermed + size, size * 2, MPI_FLOAT, child, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            merged = (Point*) malloc (size * 2 * sizeof(Point));
            merge(intermed, merged, size);
            size = size * 2;
            
            free(intermed); 

            currentStep = currentStep + currentStep;

        } else {
              MPI_Send(merged, size * 2, MPI_FLOAT, parent, 0, MPI_COMM_WORLD);
              free(merged);  
              return NULL;
        }
    }
    return merged;
}

int main(int argc, char* argv[])
{
    if (argc < 4) {
        printf("FORMAT: %s <num points> <input file> <oversampling>\n", argv[0]);
        return false;
    }

    int myrank, numranks;
    MPI_Init(&argc, &argv); // init MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);

    size_t numPoints = 0;
    if (1 != sscanf(argv[1], "%zu", &numPoints))
        return false;
    int oversampling = atoi(argv[3]);
    printf("Number elements is %zu\n", numPoints);
    size_t myPointsNum = numPoints / numranks;
    // add stragglers to the last rank
    myPointsNum += (numranks == myrank + 1) * (numPoints % myPointsNum);


    Point * points = malloc(myPointsNum * sizeof(Point));
    readInData(argv[2], points, myrank, numPoints / numranks, myPointsNum);
    localSort(points);

    merge(myrank, numrank, points);

    return true;
}
