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

float return_polar_angle(Point p1, Point p2){
	float d_x = (p1.y)-(p2.x);
	float d_y = (p1.y)-(p2.y);
	return atan2(d_y,d_x);
}

void Perform_Graham(Point* point_set){
	/*
	have list of poitns
	create empty stack

	find lowest left most point
	sort coordinates by polar angle

	loop through points
		keep removing pionts while orientation is not counterclockwise
		
	once have set, repeat process and find actual outer hull

	*/


}

Point bottom_left(Point * set_points, int num_points){
	float min_y = 1000000;
	Point min_point
	int i;
	for (i = 0; i < num_points; i++){
		printf("%f\n", (set_points[i].x));
		if (set_points[i].y < min_y){
			min_point = set_points[i]
		}
		else if(set_points[i].x < min_point.x){
			min_point = set_points[i]
		}
	}
	return min_point
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

	bottom_left(points, numPoints);

	return 0;
}