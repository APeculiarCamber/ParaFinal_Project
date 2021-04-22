#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdbool.h>
#include <string.h>
#include "mpi.h"

#define LEAD_RANK 0
#define STARTING_VECTOR_SIZE 64
#define CLOCK_RATE 512000000.0 // for clock to seconds conversion

typedef unsigned long long ticks;

// Get clock time
static __inline__ ticks getticks(void) {
    unsigned int tbl, tbu0, tbu1;
    do {
        __asm__ __volatile__ ("mftbu %0" : "=r"(tbu0));
        __asm__ __volatile__ ("mftb %0" : "=r"(tbl));
        __asm__ __volatile__ ("mftbu %0" : "=r"(tbu1));
    } while (tbu0 != tbu1);
    return ((((unsigned long long)tbu0) << 32) | tbl);
}

// Custom Struct Point
struct Point {
    float x;
    float y;
};
typedef struct Point Point;

// Create a custom MPI_Datatype MPI_POINT
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

// print an error and exit if a return code is not MPI_SUCCESS
void ensureReturnCode(int rc, char* point) {
    if (rc != MPI_SUCCESS) {
        printf("An error (%d) occured when %s\n", rc, point);
        exit(-1);
    }
}

// read in this MPI rank's chunk of the data
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

//
int comparePoints(const void* a, const void* b) {
  Point * ap = (Point*)a;
  Point * bp = (Point*)b;
  int yComp = (ap->y > bp->y) - (ap->y < bp->y);
  int xComp = (ap->x > bp->x) - (ap->x < bp->x);
  return xComp == 0 ? yComp : xComp;
}

// perform a local serial sort of points
void localSort(Point * points, size_t numPoints) {
    qsort (points, numPoints, sizeof(Point), comparePoints);
}

// Pick the proposed pivots for samplesort, oversample by (oversampling - 1) 
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

// find the pivots for samplesort
void extractPivotBins(Point * pivots, int numranks, int oversampling, Point * finalPivots) {
    for (int n = 0; n < numranks - 1; ++n) {
        finalPivots[n] = pivots[(oversampling * (n+1)) - 1];
    }
}

// Perform a binary search for the hypothetical location of p within points
int binSearch(Point * points, int numPoints, Point p) {

    int l = 0;
    int r = numPoints - 1;
    int m = 0;
    while (l <= r) {
        m = l + ((r - l) / 2);
        if (comparePoints(&points[m], &p) == 0)
            return m;
        if (comparePoints(&points[m], &p) < 0)
            l = m + 1;
        else
            r = m - 1;
    }
    return l;
}

// Sort a rank's local points in bins that can be redistributed to other ranks as part of sample sort
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

// VECTOR OPERATION : CONSTRUCT
Point * vectorInit(int * cap, int * size) {
    *cap = STARTING_VECTOR_SIZE;
    *size = 0;
    return malloc(STARTING_VECTOR_SIZE * sizeof(Point));
}

// VECTOR OPERATION : POP
Point * vectorPop(Point * points, int * cap, int * size) {
    *size = (*size) - 1;
    return points;
}

// VECTOR OPERATION : ADD
Point * vectorAdd(Point * points, int * cap, int * size, Point added) {
    if ((*cap) <= (*size)) {
        *cap = (*cap) * 2;
        points = realloc(points, (*cap) * sizeof(Point));
    }
    points[(*size)++] = added;
    return points;
}

// Concat top and bottom hulls together, the last elements of both are duplicates and should not be added
Point * ConcatHulls(Point * hull1, int size1, Point * hull2, int size2, unsigned * hullSize) {
    *hullSize = size1 + size2 - 2;
    Point * concat = malloc((*hullSize) * sizeof(Point));
    for (int p = 0; p < size1 - 1; ++p)
        concat[p] = hull1[p];
    for (int p = 0; p < size2 - 1; ++p)
        concat[p + size1 - 1] = hull2[p];

    return concat;
}

// triple cross
float cross(Point a, Point b, Point c) {
    return (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x); 
}

// reverse the contents of a points array
void reversePointsArray(Point * pts, int size) {
    int l = 0, r = size - 1;
    while (l < r) {
        Point t = pts[l];
        pts[l++] = pts[r];
        pts[r--] = t;
    }
}

// The LEAD RANK, performs hull construction, calling each rank in sorted order
Point * performHullContstruction(Point * myPoints, unsigned myNumPoints, int numranks, unsigned * hullSize) {
    int bottomCap, bottomSize;
    Point * bottomHull = vectorInit(&bottomCap, &bottomSize);
    int topCap, topSize;
    Point * topHull = vectorInit(&topCap, &topSize);
    // TRAVERSE THE TOP
    Point * oldPoints;
    Point * currentPoints = NULL;
    Point * nextPoints = myPoints;
    unsigned int numPoints = myNumPoints;
    unsigned int nextNumPoints = myNumPoints;
    MPI_Request dataRequest;
    for (int r = 0; r < numranks; ++r) {
        // set the next set of points to work on, start a nonblocking recv for the next block
        oldPoints = currentPoints;
        currentPoints = nextPoints;
        numPoints = nextNumPoints;
        if (r + 1 < numranks) {
            // get the next array (will become the working array after this iteration)
            if (oldPoints) free(oldPoints);
            // get rank points size, THIS RECV UNBLOCKS AN MPI RANK, CAUSING IT TO SEND ALL ITS POINTS TO LEAD RANK
            MPI_Recv(&nextNumPoints, 1, MPI_UNSIGNED, r + 1, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            nextPoints = malloc(nextNumPoints * sizeof(Point));
            // non-blocking, let the rank send while we perform hull construction using the current array
            MPI_Irecv(nextPoints, nextNumPoints, MPI_POINT, r + 1, 0,
                 MPI_COMM_WORLD, &dataRequest);
        }

        // construct portion of the hull for the current points
        for (int p = 0; p < numPoints; ++p) {
            while (bottomSize >= 2 && cross(bottomHull[bottomSize - 2], bottomHull[bottomSize - 1], currentPoints[p]) <= 0) {
                bottomHull = vectorPop(bottomHull, &bottomCap, &bottomSize);
            }
            bottomHull = vectorAdd(bottomHull, &bottomCap, &bottomSize, currentPoints[p]);

            while (topSize >= 2 && cross(topHull[topSize - 2], topHull[topSize -1], currentPoints[p]) >= 0) {
                topHull = vectorPop(topHull, &topCap, &topSize);
            }
            topHull = vectorAdd(topHull, &topCap, &topSize, currentPoints[p]);
        }

        // wait til we fully recieve the next ordered points set
        if (r + 1 < numranks)
            MPI_Wait(&dataRequest, MPI_STATUS_IGNORE);
    }
    // top hull must be reversed
    reversePointsArray(topHull, topSize);
    Point * concat = ConcatHulls(topHull, topSize, bottomHull, bottomSize, hullSize);

    free(currentPoints);
    free(topHull);
    free(bottomHull);
    return concat;
}

void prepareToSendToLeadRank(Point * myPoints, int numPoints, int myrank, int numranks) {
    // send and block on message on size of this rank's points
    MPI_Send(&numPoints, 1, MPI_UNSIGNED, LEAD_RANK, 0, MPI_COMM_WORLD);
    // send all points after the LEAD RANK makes a request for this rank's size
    MPI_Send(myPoints, numPoints, MPI_POINT, LEAD_RANK, 0, MPI_COMM_WORLD);
}

// Output times and results, store results in a file
void outputResults(Point * hull, unsigned hullSize, 
    ticks startFromReadin, ticks startFromAlgo, ticks finishTime) {
    printf("The time from file IO to end was %f\n", 
        (double)(finishTime - startFromReadin) / CLOCK_RATE);
    printf("The time from after file IO to end was %f\n", 
        (double)(finishTime - startFromAlgo) / CLOCK_RATE);

    // write hull to file
    FILE *f = fopen("mono_hull_results.bin", "wb");
    fwrite(hull, sizeof(Point), hullSize, f);
    fclose(f);

#ifdef DEBUG
    printf("\nTHE CONVEX HULL IS:\n");
    for (int p = 0; p < hullSize; ++p) {
        printf("{%f, %f}\n", hull[p].x, hull[p].y);
    }
#endif
}

// makes pivots and exchanges points with other ranks based on them
// WARNING : FREES THE POINTS ARRAY PASSED IN
// WARNING : SHOULD NOT BE USED IF THERE IS ONLY 1 RANK
Point * makeAndExchangeByPivots(Point * points, unsigned numPoints, unsigned oversampling,
    int myrank, int numranks, size_t * outTotal) {
    // get pivots
    Point * pivots = malloc(numranks * oversampling * sizeof(Point));
    Point * localPivots = malloc(oversampling * sizeof(Point));
    pickPivots(points, numPoints, localPivots, oversampling);
    // All Gather pivots
    int rc = MPI_Allgather(
        localPivots, oversampling, MPI_POINT, 
        pivots, oversampling, MPI_POINT, 
        MPI_COMM_WORLD);
    free(localPivots);
    // local sort all the gathered pivots
    localSort(pivots, numranks * oversampling);
    // extract the real pivots from the oversampled, sorted array
    Point * finalPivots = malloc((numranks - 1) * sizeof(Point));
    extractPivotBins(pivots, numranks, oversampling, finalPivots);
    free(pivots);
    // place into bins based on pivots
    int * steps = calloc(numranks, sizeof(int));
    int * sizes = calloc(numranks, sizeof(int));
    Point * bins = calloc(numPoints, sizeof(Point));
    sortIntoBins(points, bins, steps, sizes, finalPivots, numPoints, numranks, myrank);
    free(finalPivots); free(points);
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
    free(steps); free(sizes); free(bins);
    free(allMyBinDispl); free(allMyBinSizes);
    // sort my bin
    localSort(myREALPoints, total);
    *outTotal = total;
    return myREALPoints;
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

    size_t totalPoints = 0;
    if (1 != sscanf(argv[1], "%zu", &totalPoints))
        return false;
    unsigned oversampling = atoi(argv[3]);
    if (myrank == LEAD_RANK)
        printf("*Number of points is %zu run over %d ranks with %d oversampling.\n",
         totalPoints, numranks, oversampling);
    size_t numPoints = totalPoints / numranks;
    // add stragglers to the last rank
    numPoints += ((numranks == myrank + 1) * (totalPoints % numPoints));

    // read in the data and start the timers
    Point * points = malloc(numPoints * sizeof(Point));
    ticks startFromReadin = getticks();
    readInData(argv[2], points, myrank, (totalPoints / numranks) * sizeof(Point), numPoints);
    ticks startFromAlgo = getticks();

    // make points and exchange them with other ranks, if there are no other ranks, just local sort
    Point * myREALPoints = points;
    size_t total = numPoints;
    if (numranks > 1)
        myREALPoints = makeAndExchangeByPivots(points, numPoints, oversampling, myrank, numranks, &total);
    else
        localSort(myREALPoints, numPoints);

#ifdef DEBUG
    for (int r = 0; r < numranks; ++r) {
        if (r == myrank) {
            for (int i = 0; i < total; ++i)
                printf("Rank %d: {%f, %f}\n", myrank, myREALPoints[i].x, myREALPoints[i].y);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    printf("\n");
#endif

    MPI_Barrier(MPI_COMM_WORLD);
    // perform the monotone chain
    if (myrank == LEAD_RANK) {
        unsigned int hullSize;
        Point * hull = performHullContstruction(myREALPoints, total, numranks, &hullSize);
        ticks finishTime = getticks();
        outputResults(hull, hullSize, startFromReadin, startFromAlgo, finishTime);
    } else {
        prepareToSendToLeadRank(myREALPoints, total, myrank, numranks);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    return true;
}
