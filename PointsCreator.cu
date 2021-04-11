#include<unistd.h>
#include<stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<cuda.h>
#include<cuda_runtime.h>
#include <curand_kernel.h>

struct Point {
    float x;
    float y;
};

__global__ void data_kernel(Point* p_data, size_t numPoints,
    float leftX, float upperY, float rightX, float lowerY) {

    curandState_t state;
    curand_init (0x2340238, 0, 0, &state);
    float sizeX = rightX - leftX;
    float sizeY = upperY - lowerY;

    size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    float x, y;
    while (index < numPoints) {
        x = curand_uniform (&state);
        y = curand_uniform (&state);
        p_data[index].x = (x * sizeX) + leftX;
        p_data[index].y = (y * sizeY) + lowerY;
        index += blockDim.x * gridDim.x;
    }
}


// TODO : add MPI and parallel file processing for even larger files
int main(int argc, char* argv[])
{
    if (argc < 7) {
        printf("FORMAT: %s <num points> <num threads> <left x> <upper y> <right x> <lower y>", argv[0]);
        return false;
    }
    size_t numPoints = 0;
    if (1 != sscanf(argv[1], "%zu", &numPoints))
        return false;
    int threadsCount = atoi(argv[2]);
    float leftX = atof(argv[3]);
    float upperY = atof(argv[4]);
    float rightX = atof(argv[5]);
    float lowerY = atof(argv[6]);
    printf("Number elements is %zu\n", numPoints);
    
    Point * p_data;
    cudaMallocManaged((void**)&p_data, numPoints * sizeof(Point));


    size_t numBlocks = ((numPoints) + (threadsCount - 1)) / threadsCount;
    numBlocks = (numBlocks > 65535) ? 65535 : numBlocks;
    data_kernel<<<numBlocks, threadsCount>>>(p_data, numPoints,
        leftX, upperY, rightX, lowerY);
#ifdef DEBUG
    for (int p = 0; p < numPoints; ++p)
        printf("{%f, %f}\n", p_data[p].x, p_data[p].y);
#endif
    FILE* f_ptr = fopen("points.bin", "wb");
    fwrite(p_data, sizeof(Point), numPoints, f_ptr);

    return true;
}
