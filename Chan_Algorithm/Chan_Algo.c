#define _GNU_SOURCE
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdbool.h>
#include<string.h>
#include<math.h>
#include"mpi.h"

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

Tuple* polar_points(Point * set_points, int num_points, Point min_point, int myrank){
	int i;
	Tuple* point_list = calloc(num_points-1, sizeof(Tuple));
	int count = 0;
	for (i = 0; i < num_points; i++){
		if (set_points[i].y != min_point.y && set_points[i].x != min_point.x){
			Point temp = set_points[i];
			float temp_polar = return_polar_angle(temp, min_point);
			Tuple temp_tuple;
			temp_tuple.p = temp;
			temp_tuple.angle = temp_polar;
			point_list[count] = temp_tuple;
			count = count+1;
		}
	}
	return point_list;
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

int distance_compare(Point p1, Point p2, Point min_point){
	float d1 = sqrt(pow((p1.x-min_point.x),2)+pow((p1.y-min_point.y),2));
	float d2 = sqrt(pow((p2.x-min_point.x),2)+pow((p2.y-min_point.y),2));

	return d1-d2;
}

int compare (const void *a, const void * b, void * c){
	Point *point_a = (Point *)a;
	Point *point_b = (Point *)b;
	Point *point_c = (Point *)c;

	float or = (point_b->y - point_a->y)*(point_c->x-point_b->x) - (point_c->y-point_b->y)*(point_b->x-point_a->x);

	if (or == 0){
		float d1 = sqrt(pow((point_b->x-point_a->x),2)+pow((point_b->y-point_a->y),2));
		float d2 = sqrt(pow((point_c->x-point_a->x),2)+pow((point_c->y-point_a->y),2));
		if (d1 > d2){
			return -1;
		}
		else{
			return 1;
		}
	}

	if (or > 0){
		return 1;
	}
	return -1;
}

int item_exists(Point* exist_list, Point current, Point min_point, int numPoints){
	int i;
	for (i = 0; i < numPoints; i++){
		float or = orientation(min_point, current, exist_list[i]);
		if (or == 0){
			return i;
		}
	}
	return -1;
}

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

void pop_from_stack(Polar_s* ps){
	Point empty;
	ps->stack[ps->size] = empty;
	ps->size = ps->size-1;
}

void add_to_stack(Polar_s* ps, Point new_point, int current_size){
	ps->size = ps->size+1;
	ps->stack[ps->size] = new_point;
}

Polar_s Perform_Graham(Point* point_set, int set_size, int myrank){
	Polar_s hull;
	int current_size = 8;
	hull.stack = calloc(current_size, sizeof(Point));
	hull.size = 3;

	hull.stack[0] = point_set[0];
	hull.stack[1] = point_set[1];
	hull.stack[2] = point_set[2];

	int j;

	/*printf("before iteration\n");
	for (j = 0; j < 3; j++){
		printf("%d point is x:%f, y%f \n", myrank, hull.stack[j].x, hull.stack[j].y);
	}*/
	

	int i;
	for (i = 3; i < set_size; i++){
		Point prev = point_set[hull.size-1];
		Point curr = hull.stack[hull.size];
		Point next = point_set[i];
		if (myrank == 0){
			printf("next is x:%f, y:%f\n", next.x, next.y);
			printf("new iteration\n");
			for (j = 0; j < hull.size; j++){
				printf("%d current hull x:%f, y:%f. myrank %d\n", j, hull.stack[j].x, hull.stack[j].y, myrank);
			}
			printf("\n");
		}
		bool right_turn = clockwise_turn(prev, curr, next);
		if (right_turn == false){
			if (hull.size == current_size){
				current_size = current_size*2;
				realloc(hull.stack, current_size*sizeof(Polar_s));
			}
			add_to_stack(&hull, next, hull.size);
		}
		else{
			while (right_turn == true){
				pop_from_stack(&hull);
				int current = hull.size;
				curr = hull.stack[current];
				prev = hull.stack[current-1];
				if (myrank == 0){
					printf("after removal next is %f %f\n", next.x, next.y);
					printf("after removal curr is x:%f, y:%f\n", curr.x, curr.y);
					printf("after removal prev is x:%f, y:%f\n", prev.x, prev.y);
					printf("\n");
					for (j = 0; j < hull.size; j++){
						printf("%d current hull x:%f, y:%f. myrank %d\n", j, hull.stack[j].x, hull.stack[j].y, myrank);
					}
					printf("\n");
				}
				right_turn = clockwise_turn(prev, curr, next);
			}
			add_to_stack(&hull, next, hull.size);
		}
	}
	return hull;
}

retType perform_chan(Point* hull_set, Point min_point, int m){
	int current_size = 8;
	retType final;

	final.final = calloc(current_size, sizeof(Point));

	Point initial;
	initial.x = -INFINITY;
	initial.y = 0;

	final.final[0] = min_point;

	int i;
	int count = 1;
	float max_orientation = 0;
	Point max_point;
	for (i = 0 ; i < m; i++){
		Point current = hull_set[i];
		if (count == 1){
			float curr_orientation = orientation(initial, min_point, current);
			if (curr_orientation > max_orientation){
				max_point = current;
			}
		}
	}

	final.final[count] = max_point;
	//count = count+1;

	final.count = 1;

	while(max_point.x != final.final[0].x && max_point.y != final.final[0].y){
		Point max_point;
		for (i = 0 ; i < m; i++){
			Point current = hull_set[i];
			float curr_orientation = orientation(final.final[i-1], final.final[i], current);
			if (curr_orientation > max_orientation){
				max_point = current;
			}
		}
		count = count+1;
		final.count = final.count+1;
		if (count == current_size){
			current_size = current_size*2;
			realloc(final.final, current_size*sizeof(Polar_s));
		}
		final.final[count] = max_point;	
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

	Point min_point = bottom_left(points, numPoints, &location);

	Point temp = points[0];
	points[0] = min_point;
	points[location] = temp;

	qsort_r(points, numPoints, sizeof(struct Point), compare, &min_point);

	int j;

	int count = 0;

	Point* final_point_set = remove_repeat(points, min_point, numPoints, myrank, &count);
	
	if (myrank == 0){
		for (int i = 0; i < count; i++){
			printf("%f, %f\n", final_point_set[i].x, final_point_set[i].y);
		}
	}

	if (count < 3){
		return 0;
	}

	Polar_s hull = Perform_Graham(final_point_set, count, myrank);

	int final_size = 0;

	Point* set = calloc(hull.size, sizeof(Point));

	//printf("%d size is %d\n", myrank, hull.size);

	for (i = 0; i < hull.size; i++){
		//printf("outer hull is x:%f, y:%f, rank is %d\n", hull.stack[i].x, hull.stack[i].y, myrank);
		set[i] = hull.stack[i];
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

	/*
	if (myrank == 0){

		for (i = 0; i < final_size; i++){
			printf("point is x: %f, y:%f\n", hull_gather[i].x, hull_gather[i].y);
		}
	
		printf("final_size is %d\n", final_size);
		

		Point bottom_point = bottom_left(hull_gather, final_size);

		retType final_hull = perform_chan(hull_gather, bottom_point, final_size);

		/*
		for (int i = 0; i < final_hull.count; i++){
			printf("final hull is x:%f, y:%f", final_hull.final[i].x, final_hull.final[i].y);
		}
		
	}
	*/
	return 0;
}