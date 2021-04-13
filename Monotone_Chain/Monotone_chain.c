#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdbool.h>
#include <string.h>
#include "mpi.h"

#define LEAD_RANK 0
#define STARTING_VECTOR_SIZE 64

// TODO : file offsets appear to be in bytes?

struct Point {
    float x;
    float y;
};
typedef struct Point Point;

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


void localSort(Point * points, size_t numPoints) {
    qsort (points, numPoints, sizeof(Point), comparePoints);
}

void pickPivots(Point * points, size_t numPoints, Point * pivots, int oversampling) {
    if (oversampling >= numPoints) {
        printf("ERROR: oversampling too many for sub-data of size %ld\n", numPoints);
        exit(-1);
    }
    for (int o = 0; o < oversampling; ++o) {
        int r = rand() % numPoints;
        pivots[o] = points[r];
        // swap to preserve contents and allow rand decr
        Point p = points[r];
        points[r] = points[numPoints - 1];
        points[numPoints - 1] = p;
        numPoints -= 1;
    }
}

void extractPivotBins(Point * pivots, int numranks, int oversampling, Point * finalPivots) {
    for (int n = 0; n < numranks - 1; ++n) {
        finalPivots[n] = pivots[(oversampling * (n+1)) - 1];
    }
}

int binSearch(Point * points, int numPoints, Point p) {
    // TODO : I've got a bad track record with bin search, test this
    int i;
    for (i = 0; i < numPoints; ++i)
        if (comparePoints(&p, &points[i]) <= 0)
            return i;
    return i;

    /*
    int l, r, m;
    l = 0, r = numPoints - 1;

    while (l <= r) {
        int m = l + (r - l) / 2;
  
        if (comparePoints(&points[m], &p) == 0) {
            return m;
        }
  
        if (comparePoints(&points[m], &p) < 0)
            l = m + 1;
        else
            r = m - 1;
    }
    return l;*/
}

void sortIntoBins(Point * points, Point * bins, int * steps, int * sizes,
    Point * finalPivots, size_t numPoints, int numranks, int myrank) {

    int * stepsTemp = calloc(numranks, sizeof(int));
    // count size of bins
    for (int p = 0; p < numPoints; ++p) {
        sizes[binSearch(finalPivots, numranks - 1, points[p])]++;
    }

    for (int s = 1; s < numranks; ++s) {
        stepsTemp[s] = (steps[s] = steps[s-1] + sizes[s-1]);
    }
    
    int bin;
    for (int p = 0; p < numPoints; ++p) {
        bin = binSearch(finalPivots, numranks - 1, points[p]);
        bins[stepsTemp[bin]] = points[p];
        stepsTemp[bin]++;
    }
    free(stepsTemp);
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

Point * ConcatHulls(Point * p1, int size1, Point * p2, int size2, unsigned * hullSize) {
    *hullSize = size1 + size2 - 2;
    Point * concat = malloc((*hullSize) * sizeof(Point));
    for (int p = 0; p < size1 - 1; ++p)
        concat[p] = p1[p];
    for (int p = 0; p < size2 - 1; ++p)
        concat[p + size1 - 1] = p2[p];

    return concat;
}

float cross(Point a, Point b, Point c) {
    return (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x); 
}

// TODO : test on case where x values of 2 floats are identical...
// TODO : improve message passing to block on first size get, them use non-blocking to perform ops while recving??????
Point * performHullContstruction(Point * myPoints, unsigned myNumPoints, int numranks, unsigned * hullSize) {
    int bottomCap, bottomSize;
    Point * bottomHull = vectorInit(&bottomCap, &bottomSize);

    // TRAVERSE THE TOP
    Point * borrowedPoints = myPoints;
    unsigned int numPoints = myNumPoints;
    for (int r = 0; r < numranks; ++r) {
        for (int p = 0; p < numPoints; ++p) {
            while (bottomSize >= 2 && cross(bottomHull[bottomSize - 2], bottomHull[bottomSize - 1], borrowedPoints[p]) <= 0) {
                bottomHull = vectorPop(bottomHull, &bottomCap, &bottomSize);
            }
            bottomHull = vectorAdd(bottomHull, &bottomCap, &bottomSize, borrowedPoints[p]);
        }

        if (r + 1 < numranks) {
            printf("Rank 0: Retrieving data from rank %d\n", r + 1);
            // get the next array

            if (borrowedPoints != myPoints)
                free(borrowedPoints);
            MPI_Recv(&numPoints, 1, MPI_UNSIGNED, r + 1, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            borrowedPoints = malloc(numPoints * sizeof(Point));
            MPI_Recv(borrowedPoints, numPoints, MPI_POINT, r + 1, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    // TRAVERSE THE BOTTOM, TODO : going the other direction ( > 0) might work, but the math is poorly understood rn
    int topCap, topSize;
    Point * topHull = vectorInit(&topCap, &topSize);
    
    for (int r = numranks - 1; r >= 0; --r) {
        for (int p = numPoints - 1; p >= 0; --p) {
            while (topSize >= 2 && cross(topHull[topSize - 2], topHull[topSize - 1], borrowedPoints[p]) <= 0) {
               topHull = vectorPop(topHull, &topCap, &topSize);
            }
            topHull = vectorAdd(topHull, &topCap, &topSize, borrowedPoints[p]);
        }


        // get the next array
        if (r - 1 == 0) {
            free(borrowedPoints);
            numPoints = myNumPoints;
            borrowedPoints = myPoints;
        } else {
        MPI_Recv(&numPoints, 1, MPI_UNSIGNED, r - 1, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        free(borrowedPoints);
        borrowedPoints = malloc(numPoints * sizeof(Point));
        MPI_Recv(borrowedPoints, numPoints, MPI_POINT, r - 1, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    Point * concat = ConcatHulls(topHull, topSize, bottomHull, bottomSize, hullSize);
    // free(borrowedPoints);
    free(topHull);
    free(bottomHull);
    return concat;
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
        printf("FORMAT: %s <num points> <input file> <oversampling>\n", argv[0]);
        return false;
    }

    int myrank, numranks;
    MPI_Init(&argc, &argv); // init MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);
    createPointType();

    size_t totlalPoints = 0;
    if (1 != sscanf(argv[1], "%zu", &totlalPoints))
        return false;
    int oversampling = atoi(argv[3]);
    printf("Number elements is %zu\n", totlalPoints);
    size_t numPoints = totlalPoints / numranks;
    // add stragglers to the last rank
    numPoints += ((numranks == myrank + 1) * (totlalPoints % numPoints));


    Point * points = malloc(numPoints * sizeof(Point));
    readInData(argv[2], points, myrank, (totlalPoints / numranks) * sizeof(Point), numPoints);

    // get pivots
    Point * pivots = malloc(numranks * oversampling * sizeof(Point));
    Point * localPivots = malloc(oversampling * sizeof(Point));
    pickPivots(points, numPoints, localPivots, oversampling);
    // All Gather pivots
    int rc = MPI_Allgather(
        localPivots, oversampling, MPI_POINT, 
        pivots, oversampling, MPI_POINT, 
        MPI_COMM_WORLD);

    localSort(pivots, numranks * oversampling);

    Point * finalPivots = malloc((numranks - 1) * sizeof(Point));
    extractPivotBins(pivots, numranks, oversampling, finalPivots);
    free(pivots);

    // place into bins based on pivots, TODO : see about in-place
    int * steps = calloc(numranks, sizeof(int));
    int * sizes = calloc(numranks, sizeof(int));
    Point * bins = calloc(numPoints, sizeof(Point));
    sortIntoBins(points, bins, steps, sizes, finalPivots, numPoints, numranks, myrank);

    // send respective bin sizes to others
    int * allMyBinSizes = malloc(numranks * sizeof(int));
    int * allMyBinDispl = malloc(numranks * sizeof(int));
    MPI_Alltoall(sizes, 1, MPI_INT, allMyBinSizes, 1, MPI_INT, MPI_COMM_WORLD);

    // exhange bins with other processes
    size_t total = 0;
    for (int i = 0; i < numranks; i++) {
        allMyBinDispl[i] = total;
        total += allMyBinSizes[i];
    }

    Point * myREALPoints = malloc(total * sizeof(Point));
    printf("Rank %d's bin size is %lu\n", myrank, total);
    MPI_Alltoallv(bins, sizes, steps, MPI_POINT, myREALPoints,
                  allMyBinSizes, allMyBinDispl, MPI_POINT, MPI_COMM_WORLD);
    
    // sort my bin
    localSort(myREALPoints, total);

    for (int r = 0; r < numranks; ++r) {
        if (r == myrank) {
            for (int i = 0; i < total; ++i)
                printf("Rank %d: {%f, %f}\n", myrank, myREALPoints[i].x, myREALPoints[i].y);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    printf("\n");

    // perform the monotone chain
    if (myrank == LEAD_RANK) {
        unsigned int hullSize;
        Point * hull = performHullContstruction(myREALPoints, total, numranks, &hullSize);
        outputHull(hull, hullSize);
    } else {
        prepareToSendToLeadRank(myREALPoints, total, myrank, numranks);
    }

    // TODO : PERFORM A TON OF FREES
    MPI_Barrier(MPI_COMM_WORLD);
    return true;
}