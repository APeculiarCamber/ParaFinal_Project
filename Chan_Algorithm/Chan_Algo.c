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

struct Point_stack{
	struct Point* stack;
	int size;
	int count;
};

struct retType{
	struct Point* final;
	int count;
};

typedef struct Point_stack Polar_s;
typedef struct Point Point;
typedef struct retType retType;

#define CLOCK_RATE 512000000.0
typedef unsigned long long ticks;

static __inline__ ticks getticks(void){
	unsigned int tb1, tbu0, tbu1;
	do{
		__asm__ __volatile__ ("mftbu %0" : "=r"(tbu0));
		__asm__ __volatile__ ("mftb %0" : "=r"(tb1));
		__asm__ __volatile__ ("mftbu %0" : "=r"(tbu1));
	}while(tbu0 != tbu1);
	return ((((unsigned long long)tbu0) << 32) | tb1);
}

void readInData(char* fName, Point * points, int myrank, size_t stride, size_t numPoints){
	MPI_File mpiFile;
	MPI_Status stat;
	//after open, have offset be based upon rank
	//printf("starting point is %zu at %d\n", myrank*numPoints*sizeof(Point), myrank);
	int rc = MPI_File_open(MPI_COMM_WORLD, fName, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiFile);
	//printf("rc is %d\n", rc);
	rc = MPI_File_read_at(mpiFile, myrank*numPoints*sizeof(Point), points, 2*numPoints, MPI_FLOAT, &stat);
	//printf("rc is %d\n", rc);
	rc = MPI_File_close(&mpiFile);
	//printf("rc is %d\n", rc);
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


MPI_Datatype MPI_POINT;
void createPointType(){
	int blockLenghts[2] = {1,1};
	MPI_Datatype types[2] = {MPI_FLOAT, MPI_FLOAT};
	MPI_Aint offsets[2];

	offsets[0] = 0;
	offsets[1] = sizeof(float);
	MPI_Type_create_struct(2, blockLenghts, offsets, types, &MPI_POINT);
	MPI_Type_commit(&MPI_POINT);
}

float orientation(Point p1, Point p2, Point p3){
	return (p2.y - p1.y)*(p3.x-p2.x) - (p3.y-p2.y)*(p2.x-p1.x);
}

int compare (const void *a, const void * b){
	Point *point_a = (Point *)a;
	Point *point_b = (Point *)b;

	if (point_a->x == point_b->x){
		return point_a->y < point_b->y;
	}
	if (point_a->x < point_b->x){
		return -1;
	}
	return 1;
}

retType perform_chan(Point* hull_set, int m){
	int current_size = 8;
	retType final;

	final.final = calloc(current_size, sizeof(Point));

	qsort(hull_set, m, sizeof(Point), compare);

	/*for (int i = 0; i < m ;i++){
		printf("%d hull_set is x:%f, y:%f\n", i, hull_set[i].x, hull_set[i].y);
	}*/

	final.final[0] = hull_set[0];

	int i;
	//count = count+1;

	final.count = 1;

	Point next = hull_set[1];

	while(next.x != final.final[0].x && next.y != final.final[0].y){
		next = hull_set[1];
		/*printf("new iteration\n");
		for (i = 0; i < final.count; i++){
			printf("%d current hull is x:%f, y:%f\n", i, final.final[i].x, final.final[i].y);
		}
		printf("\n");*/
		int curr = final.count;
		//printf("next is x:%f, y:%f\n", next.x, next.y);
		//printf("prior is x:%f, y:%f\n", final.final[curr-1].x, final.final[curr-1].y);
		for (i = 0 ; i < m; i++){
			Point current = hull_set[i];
			if (current.x != final.final[curr-1].x && current.y != final.final[curr-1].y){
				if (current.x != next.x && current.y != next.y){
					float curr_orientation = orientation(final.final[curr-1], current, next);
					if (curr_orientation < 0){
						next = current;
					}
				}
			}
		}
		final.count = final.count+1;
		if (final.count == current_size){
			//printf("before realloc\n");
			int j;
			/*for (j = 0; j < final.count; j++){
				printf("%d current hull is x:%f, y:%f\n", j, final.final[j].x, final.final[j].y);
			}*/
			current_size = current_size*2;
			final.final = realloc(final.final, current_size*sizeof(Point));
			//printf("after realloc\n");
			/*for (j = 0; j < final.count; j++){
				printf("%d current hull is x:%f, y:%f\n", j, final.final[j].x, final.final[j].y);
			}*/
		}
		final.final[final.count-1] = next;
	}

	return final;
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
	int oversample = atoi(argv[3]);
	size_t numPoints = totalPoints/numranks;
	numPoints += ((numranks == myrank + 1)*(totalPoints%numPoints));

	Point * points = malloc(numPoints * sizeof(Point));
	ticks start = getticks();
	readInData(argv[2], points, myrank, (totalPoints/numranks)*sizeof(Point), numPoints);

	int subHullSize;

	Point * subhull = makeSubHull_Monotone(points, numPoints, &subHullSize);

	//printf("%d after monotone\n", myrank);

	int* hull_sizes = NULL;
	if (myrank == 0){
		hull_sizes = malloc(numranks *sizeof(int));
	}

	MPI_Gather(&subHullSize, 1, MPI_INT, hull_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int * disp_array = NULL;
	Point * hull_gather;
	size_t gather_size = 0;

	//printf("%d after gather\n", myrank);

	if (myrank == 0){
		for (int i = 0; i < numranks;i++){
			gather_size += hull_sizes[i];
		}

		//printf("after gather size");

		disp_array = calloc(numranks, sizeof(int));
		disp_array[0] = 0;
		for (int i = 1; i < numranks; i++){
			disp_array[i] = disp_array[i-1]+hull_sizes[i-1];
		}

		//printf("after disp");
		hull_gather = calloc(gather_size, sizeof(Point));
	}
	

	MPI_Gatherv(subhull, subHullSize, MPI_POINT, hull_gather, hull_sizes, disp_array, MPI_POINT, 0, MPI_COMM_WORLD);

	//printf("after gatherv\n");

	if (myrank == 0){
		retType final_hull = perform_chan(hull_gather, (int) gather_size);
		ticks end = getticks();
		printf("%llu to %llu\n", start,end);
		printf("The final time was %f for %zu and %d ranks\n",((double)(end-start))/CLOCK_RATE, totalPoints,numranks);
		printf("The final convex hull was of size %d\n", final_hull.count);

		free(hull_gather);
		free(disp_array);
		free(hull_sizes);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}