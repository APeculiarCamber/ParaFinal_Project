#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdbool.h>
#include <string.h>
#include "mpi.h"

#define LEAD_RANK 0
#define STARTING_VECTOR_SIZE 64

typedef struct Point {
    float x;
    float y;
} Point;

typedef struct Vector {
    Point * pts;
    int size;
} Vector;

void ensureReturnCode(int rc, char* point) {
    if (rc != MPI_SUCCESS) {
        printf("An error (%d) occured when %s\n", rc, point);
        exit(-1);
    }
}

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


void readInData(char * fileName, Point * points, 
    int myrank, size_t stride, size_t numPoints) {
    MPI_File mpiFile;
    MPI_Status stat;

    int rc = MPI_File_open(MPI_COMM_WORLD, fileName, 
        MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiFile);
    ensureReturnCode(rc, "opening file");
    rc = MPI_File_read_at(mpiFile, (myrank * stride), points,
                  numPoints, MPI_POINT, &stat);
    ensureReturnCode(rc, "reading file");
    rc = MPI_File_close(&mpiFile);
    ensureReturnCode(rc, "closing file");
}

int comparePoints(const void* a, const void* b) {
  int yComp = (((Point*)a)->y > ((Point*)b)->y) - (((Point*)a)->y < ((Point*)b)->y);
  int xComp = (((Point*)a)->x > ((Point*)b)->x) - (((Point*)a)->x < ((Point*)b)->x);
  return xComp == 0 ? yComp : xComp;
}

Point * vectorInit(int * cap, int * size) {
    *cap = STARTING_VECTOR_SIZE;
    *size = 0;
    return malloc(STARTING_VECTOR_SIZE * sizeof(Point));
}

Point * vectorPop(Point * points, int * cap, int * size) {
    *size = (*size) - 1;
    return points;
}

Point * vectorAdd(Point * points, int * cap, int * size, Point added) {
    if ((*cap) <= (*size)) {
        *cap = (*cap) * 2;
        points = realloc(points, (*cap) * sizeof(Point));
    }
    points[(*size)++] = added;
    return points;
}

Vector ConcatHulls(Vector p1, Vector p2) {
    Vector concat;
    concat.size = p1.size + p2.size - 2;
    concat.pts = malloc((concat.size) * sizeof(Point));
    for (int p = 0; p < p1.size - 1; ++p)
        concat.pts[p] = p1.pts[p];
    for (int p = 0; p < p2.size - 1; ++p)
        concat.pts[p + p1.size - 1] = p2.pts[p];

    return concat;
}

float cross(Point a, Point b, Point c) {
    return (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x); 
}

void prepareToSendToLeadRank(Point * myPoints, int numPoints, int myrank, int numranks) {
    // send and block on small size
    MPI_Send(&numPoints, 1, MPI_UNSIGNED, LEAD_RANK, 0, MPI_COMM_WORLD);
    MPI_Send(myPoints, numPoints, MPI_POINT, LEAD_RANK, 0, MPI_COMM_WORLD);
    MPI_Send(&numPoints, 1, MPI_UNSIGNED, LEAD_RANK, 0, MPI_COMM_WORLD);
    MPI_Send(myPoints, numPoints, MPI_POINT, LEAD_RANK, 0, MPI_COMM_WORLD);
}

void outputHull(Point * hull, unsigned hullSize) {
    printf("\n\nTHE CONVEX HULL IS:\n");
    for (int p = 0; p < hullSize; ++p) {
        printf("{%f, %f}\n", hull[p].x, hull[p].y);
    }
}

int main(int argc, char* argv[])
{
    if (argc < 4) {
        printf("FORMAT: %s <num points> <input file>\n", argv[0]);
        return false;
    }

    int myrank, numranks;
    MPI_Init(&argc, &argv); // init MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);
    createPointType();

    size_t totalPoints = 0;
    if (1 != sscanf(argv[1], "%zu", &totalPoints))
        return false;

    printf("Number elements is %zu\n", totalPoints);
    size_t numPoints = totalPoints / numranks;
    // add stragglers to the last rank
    if (numranks == myrank + 1) {
        numPoints += totalPoints % numPoints;
    }
    Point * points = malloc(numPoints * sizeof(Point));
    readInData(argv[2], points, myrank, (totalPoints / numranks) * sizeof(Point), numPoints);



    // TODO : PERFORM A TON OF FREES
    MPI_Barrier(MPI_COMM_WORLD);
    return true;
}
