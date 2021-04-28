#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdbool.h>
#include<string.h>
#include<math.h>
#include"mpi.h"

#define STARTING_VECTOR_SIZE 64

struct Point{
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


// print an error and exit if a return code is not MPI_SUCCESS
void ensureReturnCode(int rc, char* point) {
    if (rc != MPI_SUCCESS) {
        printf("An error (%d) occured when %s\n", rc, point);
        exit(-1);
    }
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

void readInData(char* fileName, Point * points, int myrank, size_t stride, size_t numPoints){
    MPI_File mpiFile;
    MPI_Status stat;

    int rc = MPI_File_open(MPI_COMM_WORLD, fileName, 
        MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiFile);
    ensureReturnCode(rc, "opening file");
    // TODO : problem with reading 268435456
    if (numPoints == 268435456) {
    rc = MPI_File_read_at(mpiFile, (myrank * stride), points,
                  (numPoints / 2), MPI_POINT, &stat);
    ensureReturnCode(rc, "reading file");
    rc = MPI_File_read_at(mpiFile, 
        (myrank * stride) + ((numPoints / 2) * sizeof(Point)), 
        points + (numPoints / 2),
        (numPoints / 2), MPI_POINT, &stat);
    } else {
        rc = MPI_File_read_at(mpiFile, (myrank * stride), points,
                  numPoints, MPI_POINT, &stat);
    }
    ensureReturnCode(rc, "reading file");
    rc = MPI_File_close(&mpiFile);
    ensureReturnCode(rc, "closing file");
}

// triple cross
float cross(Point a, Point b, Point c) {
    return (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x); 
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


// reverse the contents of a points array
void reversePointsArray(Point * pts, int size) {
    int l = 0, r = size - 1;
    while (l < r) {
        Point t = pts[l];
        pts[l++] = pts[r];
        pts[r--] = t;
    }
}

// Concat top and bottom hulls together, the last elements of both are duplicates and should not be added
Point * ConcatHulls(Point * hull1, int size1, Point * hull2, int size2, int * hullSize) {
    *hullSize = size1 + size2 - 2;
    Point * concat = malloc((*hullSize) * sizeof(Point));
    for (int p = 0; p < size1 - 1; ++p)
        concat[p] = hull1[p];
    for (int p = 0; p < size2 - 1; ++p)
        concat[p + size1 - 1] = hull2[p];

    return concat;
}


Point * makeSubHull_Monotone(Point * pts, size_t numPoints, int * subHullSize) {
	localSort(pts, numPoints);

	int bottomCap, bottomSize;
    Point * bottomHull = vectorInit(&bottomCap, &bottomSize);
    int topCap, topSize;
    Point * topHull = vectorInit(&topCap, &topSize);
    // TRAVERSE
    for (int p = 0; p < numPoints; ++p) {
        while (bottomSize >= 2 && cross(bottomHull[bottomSize - 2], bottomHull[bottomSize - 1], pts[p]) <= 0) {
            bottomHull = vectorPop(bottomHull, &bottomCap, &bottomSize);
        }
        bottomHull = vectorAdd(bottomHull, &bottomCap, &bottomSize, pts[p]);

        while (topSize >= 2 && cross(topHull[topSize - 2], topHull[topSize -1], pts[p]) >= 0) {
            topHull = vectorPop(topHull, &topCap, &topSize);
        }
        topHull = vectorAdd(topHull, &topCap, &topSize, pts[p]);
    }
    // top hull must be reversed
    reversePointsArray(topHull, topSize);
    Point * concat = ConcatHulls(topHull, topSize, bottomHull, bottomSize, subHullSize);

    free(topHull);
    free(bottomHull);
    return concat;
}

#define COLINEAR 0
#define CLOCKWISE 1
#define COUNTER_CLOCKWISE 2
int orientation(Point a, Point b, Point c)
{
    float val = (b.y - a.y) * (c.x - b.x) -
              (b.x - a.x) * (c.y - b.y);
      return (val >= 0) ? 
      	(val == 0 ? COLINEAR : CLOCKWISE)
      	: COUNTER_CLOCKWISE;
}
  
Point getFurthestLeft(Point * pts, int ptsSize) {
	Point l = pts[0];
	for (int p = 0; p < ptsSize; ++p)
		if (l.x > pts[p].x)
			l = pts[p];
	return l;
}

char pointsEqual(Point a, Point b) {
	return a.x == b.x && a.y == b.y;
}

Point localJarvisMarch(Point * pts, int ptsSize, Point pivot) {
    Point winning = !pointsEqual(pts[0], pivot) ? pts[0] : pts[1 % ptsSize];
    for (int i = 0; i < ptsSize; i++)
       if (orientation(pivot, pts[i], winning) == COUNTER_CLOCKWISE)
           winning = pts[i];
    return winning;
}

Point * performJarvisMarch(Point * pts, int ptsSize, int * hullSize, int myrank, int numranks) {
	
	Point * rankWinners = (myrank == 0 ? malloc(numranks * sizeof(Point)) : NULL);

	// get the furthest left point of ALL ranks
	Point leftPoint = getFurthestLeft(pts, ptsSize);
	MPI_Gather(&leftPoint, 1, MPI_POINT, rankWinners, 1, MPI_POINT, 0, MPI_COMM_WORLD);
	if (myrank == 0) {
		leftPoint = getFurthestLeft(rankWinners, numranks);
#ifdef DEBUG 
		printf("The left point is {%f, %f}\n", leftPoint.x, leftPoint.y); 
#endif
	}
	MPI_Bcast(&leftPoint, 1, MPI_POINT, 0, MPI_COMM_WORLD);

	int hullCap;
	Point * hull =  vectorInit(&hullCap, hullSize);
	Point pivot = leftPoint;
	Point localWinner;
	int i = 0;
	do {
		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Bcast(&pivot, 1, MPI_POINT, 0, MPI_COMM_WORLD);
		if (pointsEqual(pivot, leftPoint) && i != 0)  {
			break; // this is functional as we only apply assignments, but cannot handle repeats
		}

		// only the leading rank needs this, but it's not problematic to allow them all to see the hull
		hull = vectorAdd(hull, &hullCap, hullSize, pivot);
		localWinner = localJarvisMarch(pts, ptsSize, pivot);

		// the leading rank 0 should gather local winners and perform local jarvis march again
		MPI_Gather(&localWinner, 1, MPI_POINT, rankWinners, 1, MPI_POINT, 0, MPI_COMM_WORLD);
		if (myrank == 0)
			pivot = localJarvisMarch(rankWinners, numranks, pivot); // reduce rank winners to pivot
		i++;
	} while (i < ptsSize * numranks || myrank != 0);
	if (i == ptsSize * numranks) {
		printf("AN ERROR OCCURED: did not perform jarvis as expected.\n");
		exit(-1);
	}

	return hull;
}

int main(int argc, char*argv[]){

	int myrank, numranks;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &numranks);
	createPointType();

	size_t totalPoints = 0;
	if (1 != sscanf(argv[1], "%zu", &totalPoints))
		return false;

	size_t numPoints = totalPoints / numranks;
	numPoints += ((numranks == myrank + 1)*(totalPoints%numPoints));

	Point * points = malloc(numPoints * sizeof(Point));
	ticks start = getticks();
	readInData(argv[2], points, myrank, (totalPoints/numranks)*sizeof(Point), numPoints);

#ifdef DEBUG
	for (int i = 0; i < numPoints; ++i) {
		printf("P: %f, %f\n", points[i].x, points[i].y);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// construct my subhull
	int subHullSize;
	Point * subhull = makeSubHull_Monotone(points, numPoints, &subHullSize);

	// if there's only 1 rank, the subhull is it, we're done
	if (numranks == 1) { 
		ticks End = getticks();
		printf("The final time was %f\n", ((double)(End - start)) / CLOCK_RATE);
		printf("%llu to %llu\n", start, End);
#ifdef DEBUG
		if (myrank == 0) {
			printf("The final convex hull was of size %d\n", subHullSize);
			for (int i = 0; i < subHullSize; ++i) {
				printf("{%f, %f}\n", subhull[i].x, subhull[i].y);
			}
		}
#endif
		MPI_Finalize();
		return 0;
	}
	int* hull_sizes = NULL;
	if (myrank == 0)
		hull_sizes = malloc(numranks * sizeof(int));


	// Gather subhull sizes so we can gather and then evenly scatter them
	MPI_Gather(&subHullSize, 1, MPI_INT, hull_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int * disp_array = NULL;
	Point * hull_gather = NULL;
	size_t gather_size = 0;

	if (myrank == 0) {
#ifdef DEBUG
		for (int r = 0; r < numranks; ++r) {
			printf("The size for rank %d was %d\n", r, hull_sizes[r]);
		}
#endif
		for (int i = 0; i < numranks; i++)
			gather_size += hull_sizes[i];
		disp_array = calloc(numranks, sizeof(int));
		disp_array[0] = 0;
		for (int i = 1; i < numranks; i++)
			disp_array[i] = disp_array[i - 1] + hull_sizes[i - 1];

		hull_gather = calloc(gather_size, sizeof(Point));
		printf("My gather size for %zu and %d ranks is %zu\n", totalPoints, numranks, gather_size);
	}

	// Gather all subhulls to main process
	MPI_Gatherv(subhull, subHullSize, MPI_POINT, 
		hull_gather, hull_sizes, disp_array, MPI_POINT, 
		0, MPI_COMM_WORLD);

	// this gather_size / numranks will only function correctly for rank 1 so broadcast the result
	int scatterSize = gather_size / numranks;
	if (myrank == 0 && scatterSize == 0) {
		printf("ERROR: Scatter size was zero, cannot continue.\n");
		MPI_Finalize();
		return 1;
	}
	// broadcast the scattered size, any stragglers will go to rank 0
	MPI_Bcast(&scatterSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// rank 0 getting extras in the division works best, rank0Bonus == 0 for all other ranks
	int rank0Bonus = (myrank == 0) * (gather_size % numranks);
	Point * scatteredPoints = malloc((scatterSize + rank0Bonus) * sizeof(Point));
	MPI_Scatter(hull_gather + rank0Bonus, scatterSize, MPI_POINT,
    scatteredPoints + rank0Bonus, scatterSize, MPI_POINT, 0, MPI_COMM_WORLD);
    // apply the remaining bonus remainder points to rank 0
    for (int i = 0; i < rank0Bonus; ++i) {
    	scatteredPoints[i] = hull_gather[i];
    }
    scatterSize += rank0Bonus;
    // make free calls
    if (myrank == 0) { free(hull_gather); free(disp_array); free(hull_sizes); }

    // once each process has its set, perform parallel jarvis march
	int finalHullSize;
	Point * finalHull = performJarvisMarch(
		scatteredPoints, scatterSize, &finalHullSize, myrank, numranks);

	if (myrank == 0) {
		ticks end = getticks();
		printf("%llu to %llu\n", start, end);
		printf("The final time was %f for %zu and %d ranks\n", 
			((double)(end - start)) / CLOCK_RATE, totalPoints, numranks);
		printf("The final convex hull was of size %d\n", finalHullSize);
#ifdef DEBUG
		for (int i = 0; i < finalHullSize; ++i) {
			printf("{%f, %f}\n", finalHull[i].x, finalHull[i].y);
		}
#endif
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
