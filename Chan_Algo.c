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
}

typedef struct Point Point;

void readInData(char* fName, Point * points, int myrank, size_t stride, size_t numPoints){
	MPI_File mpiFile;
	MPI_Status, stat;
	int rc = MPI_File_open(MPI_COMM_WORLD, fName, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpi_file);
	rc = MPI_File_read_at(mpiFile, (myrank*stride)*2, points, numPoints*2, MPI_FLOAT, &stat);
	rc = MPI_file_colse(&mpiFile)
}

void createSamples(Point* points, int num_sets,Point ** point_set){
	//TODO: create a sample of points through this
	/*
	Pseudo code
	Randomly pick group of points

	*/
}

bool clockwise_turn(Point p1, Point p2, Point p3){
	slope_l1 = (p2.y-p1.y)/(p2.x-p1.y);
	slope_l2 = (p3.y-p2.y)/(p3.x-p2.x);

	if (slope_l1 >= slope_l2){
		return false;
	}
	else if(slope_l1 < slope_l2){
		return true;
	}
}

float return_polar_angle(Point p1, Point p2){
	float d_x = p1.x-p2.x;
	float d_y = p1.y-p2.y;
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

int main(int argc, char*argv[]){
	/*
	In main, first grab the data

	split into squares

	perform Graham on different hulls

	then find outer 
	*/

	return;
}