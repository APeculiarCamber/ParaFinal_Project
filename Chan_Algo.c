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
};

struct Point_stack{
	struct Point* stack
	int size;
};

typedef struct Tuple Tuple;
typedef struct Polar_vector Polar_l;
typedef struct Point_stack Polar_s;
typedef struct Point Point;

void readInData(char* fName, Point * points, int myrank, size_t stride, size_t numPoints){
	MPI_File mpiFile;
	MPI_Status stat;
	int rc = MPI_File_open(MPI_COMM_WORLD, fName, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiFile);
	rc = MPI_File_read_at(mpiFile, (myrank*stride)*2, points, numPoints*2, MPI_FLOAT, &stat);
	rc = MPI_File_close(&mpiFile);
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

void createSamples(Point* points, int num_sets,Point ** point_set){
	//TODO: create a sample of points through this
	/*
	Pseudo code
	Randomly pick group of points

	*/
}

bool clockwise_turn(Point p1, Point p2, Point p3){
	float slope_l1 = ((p2.y)-(p1.y))/((p2.x)-(p1.x));
	float slope_l2 = ((p3.y)-(p2.y))/((p3.x)-(p2.x));

	if (slope_l1 >= slope_l2){
		return false;
	}
	else if(slope_l1 < slope_l2){
		return true;
	}
}

float orientation(Point p1, Point p2, Point p3){
	return (p2.y - p1.y)*(p3.x-p2.x) - (p3.y-p2.y)*(p2.x-p1.x);
}

float return_polar_angle(Point p1, Point p2){
	float d_x = (p1.y)-(p2.x);
	float d_y = (p1.y)-(p2.y);
	return atan2(d_y,d_x);
}

Point bottom_left(Point * set_points, int num_points){
	float min_y = 1000000;
	Point min_point;
	int i;
	for (i = 0; i < num_points; i++){
		if (set_points[i].y < min_y){
			printf("min y: %f, min x: %f\n", min_point.x, min_point.y);
			min_y = set_points[i].y;
			min_point = set_points[i];
		}
		else if(set_points[i].x < min_point.x){
			min_point = set_points[i];
		}
	}
	return min_point;
}

Polar_l polar_points(Point * set_points, int num_points, Point min_point){
	int i;
	Polar_l point_list;
	point_list.polar_list = calloc(num_points-1, sizeof(Tuple));
	point_list.size = num_points-1;
	int count = 0;
	for (i = 0; i < num_points; i++){
		if (set_points[i].y != min_point.y && set_points[i].x != min_point.x){
			Point temp = set_points[i];
			float temp_polar = return_polar_angle(temp, min_point);
			Tuple temp_tuple;
			temp_tuple.p = temp;
			temp_tuple.angle = temp_polar;
			point_list.polar_list[count] = temp_tuple;
			count = count+1;
		}
	}
	return point_list;
}

void* add_to_polar_l(Polar_l* c_list){
	c_list->size = c_list->size*2;
	c_list->polar_list = (realloc(c_list->polar_list, c_list->size*sizeof(Tuple)));
}

int item_exists(Polar_l exist_list, float angle){
	int i;
	for (i = 0; i < exist_list.size; i++){
		if (angle == exist_list[i].angle){
			return i;
		}
	}
	return -1;
}

int distance_compare(Point p1, Point p2, min_point){
	float d1 = sqrt((p1.x-min_point.x)**2+(p1.y-min_point.y)**2)
	float d2 = sqrt((p2.x-min_point.x)**2+(p2.y-min_point.y)**2)

	return d1-d2
}

Polar_l remove_repeat(Polar_l current_list, Point min_point){
	Polar_l without_repeat;
	int current_size = 8;

	Polar_l without_repeat.polar_list = calloc(current_size, sizeof(Tuple));

	int count = 0;
	int i;
	for (i = 0; i < current_list.size; i++){
		Tuple current = current_list.polar_list[i];
		int index = item_exists(current_size,angle);
		if (index == -1){
			count += 1;
			if (count == current_size){
				add_to_polar_l(&without_repeat);
			}
			without_repeat.polar_list[count] = current;
		}
		else{
			if (distance_compare(current.p, without_repeat[index].p, min_point) < 0){
				without_repeat[index] = current;
			}
		}
	}
	return without_repeat;
}

int compare (const void *a, const void * b){
	Tuple *tupleA = (Tuple *a);
	Tuple *tupleB = (Tuple *b);

	return (tupleA->angle - tupleB->angle);
}

void pop_from_stack(Point* point_stack){
	hull->size = size-1;
	Point empty;
	hull[size] = empty;
}

void add_to_stack(Point* point_stack, Point new_point, int current_size){
	hull->size = size+1;
	hull[size-1] = new_point;
}

Polar_s Perform_Graham(Point* point_set, int set_size){
	Polar_s hull;
	int current_size = 8;
	hull.stack = calloc(current_size, sizeof(Point));
	hull.size = 3;

	hull.stack[0] = point_set[0];
	hull.stack[1] = point_set[1];
	hull.stack[2] = point_set[2];

	int i;
	for (i = 3; i < set_size; i++){
		prev = point_set[i-1];
		curr = point_set[i];
		next = point_set[i+1];
		bool right_turn = clockwise_turn(prev, curr, next);
		if (right_turn == false){
			if (hull.size == current_size){
				current_size = current_size*2;
				realloc(hull.stack, current_size*sizeof(Polar_s));
			}
			add_to_stack(&hull, next);
		}
		else{
			while (right_turn(prev,curr,next) == true){
				pop_from_stack(&hull);
				int current = hull.size-1;
				curr = hull.stack[current];
				prev = hull.stack[current-1];
			}
			add_to_stack(&hull, next);
		}
	}
	return hull;
}

Point* perform_chan(Point_set* hull_set, Point min_point, int m){
	int current_size = 8;
	Point* final = calloc(current_size, sizeof(Point));

	Point initial;
	initial.x = -INFINITY;
	initial.y = 0;

	final[0] = min_point;

	int i;
	int count = 1;
	float max_orientation = 0;
	Point max_point;
	for (i = 0 ; i < m; i++){
		Point current = hull_set[i]
		if (count == 1){
			float curr_orientation = orientation(initial, p1, current);
			if (curr_orientation > max_orientation){
				max_point = current;
			}
		}
	}

	final[count] = current;
	count = count+1;

	while(max_point != final[0]){
		Point max_point
		for (i = 0 ; i < m; i++){
			Point current = hull_set[i]
			if (count == 1){
				float curr_orientation = orientation(initial, p1, current);
				if (curr_orientation > max_orientation){
					max_point = current;
				}
			}
		}
		count = count+1
		if (count == current_size){
			current_size = current_size*2
			realloc(final, current_size*sizeof(Polar_s));
		}
		final[count] = max_point;
	}

	return final;
}

int main(int argc, char*argv[]){
	/*
	In main, first grab the data

	split into squares

	perform Graham on different hulls

	then find outer 
	*/
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

	Point min_point = bottom_left(points, numPoints);

	Polar_l p_points = polar_points(points, numPoints, min_point);

	qsort(p_points, numPoints-1, Tuple, compare);

	Polar_l final_polar_set = remove_repeat(p_points, min_point);

	if (final_polar_set.size < 3){
		return 0
	}

	Point* final_point_set = calloc(final_polar_set.size, sizeof(Point));
	int i;
	for (i = 0; i < final_set.size; i++){
		final_point_set[i] = final_polar_set.polar_list[i];
	}

	Polar_s hull = Perform_Graham(final_point_set, final_polar_set.size);

	return 0;
}