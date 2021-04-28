#define _GNU_SOURCE
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdbool.h>
#include<string.h>
#include<math.h>
#include"mpi.h"
#include""

struct Point{
	float x;
	float y;
};

struct Tuple{
	struct Point p;
	float angle;
};

struct Polar_vector{
	struct Tuple * polar_list;
	int size;
	int count;
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

typedef struct Tuple Tuple;
typedef struct Polar_vector Polar_l;
typedef struct Point_stack Polar_s;
typedef struct Point Point;
typedef struct retType retType;
typedef unsigned long long ticks;

static __inline__ ticks getticks(void){
	unsigned int tb1, tbu0, tbu1;
	do{
		__asm__ __volatile__ ("mftbu %0" : "=r"(tbu0));
		__asm__ __volatile__ ("mftbu %0" : "=r"(tb1));
		__asm__ __volatile__ ("mftbu %0" : "=r"(tbu1));
	}while(tbu0 != tbu1);
	return ((((unsigned long long)tbu0) << 32) | tb1);
}

void readInData(char* fName, Point * points, int myrank, size_t stride, size_t numPoints){
	MPI_File mpiFile;
	MPI_Status stat;
	//after open, have offset be based upon rank
	printf("starting point is %zu at %d\n", myrank*numPoints*sizeof(Point), myrank);
	int rc = MPI_File_open(MPI_COMM_WORLD, fName, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiFile);
	printf("rc is %d\n", rc);
	rc = MPI_File_read_at(mpiFile, myrank*numPoints*sizeof(Point), points, 2*numPoints, MPI_FLOAT, &stat);
	printf("rc is %d\n", rc);
	rc = MPI_File_close(&mpiFile);
	printf("rc is %d\n", rc);
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

void createSamples(Point* points, int num_sets,Point *** point_set, int num_points){
	//TODO: create a sample of points through this
	/*
	Pseudo code
	Randomly pick group of points

	*/
	*point_set = calloc(num_sets, sizeof(Point*));

	int size_of_points = num_points/num_sets;

	int i = 0;
	for (int i = 0; i < num_sets; i++){
		*point_set[i] = points+size_of_points*i;
	}
}

float return_polar_angle(Point p1, Point p2){
	float d_x = (p1.y)-(p2.x);
	float d_y = (p1.y)-(p2.y);
	return atan2(d_y,d_x);
}

Point bottom_left(Point * set_points, int num_points, int* index){
	float min_y = 1000000;
	Point min_point;
	int i;
	for (i = 0; i < num_points; i++){
		if (set_points[i].y < min_y){
			//printf("min y: %f, min x: %f\n", min_point.y, min_point.x);
			min_y = set_points[i].y;
			min_point = set_points[i];
			*index = i;
		}
		else if(set_points[i].y == min_y){
			if (set_points[i].x < min_point.x){
				min_point = set_points[i];
				*index = i;
			}
		}
	}
	return min_point;
}

float orientation(Point p1, Point p2, Point p3){
	return (p2.y - p1.y)*(p3.x-p2.x) - (p3.y-p2.y)*(p2.x-p1.x);
}

bool clockwise_turn(Point p1, Point p2, Point p3){
	float or = orientation(p1,p2,p3);

	if (or < 0){
		return false;
	}
	return true;
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

/*
Point* remove_repeat(Point* current_list, Point min_point, int numPoints, int myrank, int* amount){
	Point* without_repeat;
	int current_size = 8;

	without_repeat = calloc(current_size, sizeof(Point));
	
	without_repeat[0] = min_point;

	int count = 1;
	int i;
	for (i = 1; i < numPoints; i++){
		int j;
		Point current = current_list[i];
		int index = item_exists(without_repeat, current, min_point, numPoints);
		if (index == -1){
			if (count == current_size){
				current_size = current_size*2;
				without_repeat = realloc(without_repeat, current_size*sizeof(Point));
			}
			without_repeat[count] = current;
			count += 1;
		}
		else{
			if (distance_compare(current, without_repeat[index], min_point) < 0){
				without_repeat[index] = current;
			}
		}
	}
	*amount = count;
	return without_repeat;
}
*/

void pop_from_stack(Polar_s* ps){
	//Point empty;
	//ps->stack[ps->size] = empty;
	ps->size = ps->size-1;
}

void add_to_stack(Polar_s* ps, Point new_point, int current_size){
	ps->size = ps->size+1;
	ps->stack[ps->size] = new_point;
}

retType perform_chan(Point* hull_set, int m){
	int current_size = 8;
	retType final;

	final.final = calloc(current_size, sizeof(Point));

	qsort(hull_set, m, sizeof(Point), compare);

	for (int i = 0; i < m ;i++){
		printf("%d hull_set is x:%f, y:%f\n", i, hull_set[i].x, hull_set[i].y);
	}

	final.final[0] = hull_set[0];

	int i;
	//count = count+1;

	final.count = 1;

	Point next = hull_set[1];

	while(next.x != final.final[0].x && next.y != final.final[0].y){
		next = hull_set[1];
		printf("new iteration\n");
		for (i = 0; i < final.count; i++){
			printf("%d current hull is x:%f, y:%f\n", i, final.final[i].x, final.final[i].y);
		}
		printf("\n");
		int curr = final.count;
		printf("next is x:%f, y:%f\n", next.x, next.y);
		printf("prior is x:%f, y:%f\n", final.final[curr-1].x, final.final[curr-1].y);
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
			current_size = current_size*2;
			realloc(final.final, current_size*sizeof(Polar_s));
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
	readInData(argv[2], points, myrank, (totalPoints/numranks)*sizeof(Point), numPoints);

	int i;

	int location = 0;

	//Point min_point = bottom_left(points, numPoints, &location);

	//qsort_r(points, numPoints, sizeof(struct Point), compare, &min_point);W

	//Point* final_point_set = remove_repeat(points, min_point, numPoints, myrank, &count);

	if (myrank == 0){
		retType hull = perform_chan(points, numPoints);

		int final_size = 0;

		//Point* set = calloc(hull.size, sizeof(Point));

		//printf("%d size is %d\n", myrank, hull.size);

		for (i = 0; i < hull.count; i++){
			printf("%d hull is x:%f, y%f\n", i,hull.final[i].x, hull.final[i].y);
		}
		
	}

	//MPI_Reduce(&hull.size, &final_size, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	int* hull_sizes = calloc(numranks, sizeof(int));

	MPI_Gather(&hull.size, 1, MPI_INT, hull_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

	for (i = 0; i < numranks; i++){
		final_size = final_size+hull_sizes[i];
	}

	int* disp_array = calloc(numranks, sizeof(int));

	for (int i = 0; i < numranks; i++){
		int sum = 0;
		for (int j = 0; j < i; j++){
			sum = sum+hull_sizes[j];
		}
		disp_array[i] = sum;
	}

	Point* hull_gather = calloc(final_size, sizeof(Point));

	MPI_Gatherv(set, hull.size, MPI_POINT, hull_gather, hull_sizes, disp_array, MPI_POINT, 0, MPI_COMM_WORLD);

	return 0;
}