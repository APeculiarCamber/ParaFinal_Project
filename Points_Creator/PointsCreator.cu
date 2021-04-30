#include<unistd.h>
#include<stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<cuda.h>
#include<cuda_runtime.h>
#include <curand_kernel.h>

#define PI 3.14159265359

struct Point {
    float x;
    float y;
};
typedef struct Point Point;

// LOCAL
__global__ void data_kernel(Point* p_data, size_t numPoints, int seed,
    float leftX, float lowerY, float rightX, float upperY) {

    curandState_t state;
    curand_init (seed + (128 * blockIdx.x) + (60000 * threadIdx.x), 0, 0, &state);
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

// LOCAL
__global__ void data_circle_kernel(Point* p_data, size_t numPoints, int seed,
    float radius) {

    curandState_t state;
    curand_init (seed + (128 * blockIdx.x) + (60000 * threadIdx.x), 0, 0, &state);

    size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    while (index < numPoints) {
        float r = radius * sqrt(curand_uniform (&state));
        float t = curand_uniform (&state) * 2 * PI;

        p_data[index].x = radius + (r * cos(t));
        p_data[index].y = radius + (r * sin(t));
        index += blockDim.x * gridDim.x;
    }
}


// Prevent name mangling
extern "C" {
    void c_cudaAlloc(Point ** points, size_t numPoints);

    void c_callKernel(int numBlocks, int threadsCount, Point * p_data, size_t numPoints, int seed,
        float leftX, float lowerY, float rightX, float upperY);

    void c_callCircleKernel(int numBlocks, int threadsCount, Point * p_data, 
        size_t numPoints, int seed, float radius);

    void setCudaDeviceByRank(int myrank);

    void freeCudaMemory(Point * pts);
}

void c_cudaAlloc(Point ** points, size_t numPoints) {
    cudaMallocManaged(points, numPoints * sizeof(Point));
}

void setCudaDeviceByRank(int myrank) {
    cudaError_t cE = cudaSuccess;
    int cudaDeviceCount = -1;
    int assignedCudaDevice = -1;

    if( (cE = cudaGetDeviceCount( &cudaDeviceCount)) != cudaSuccess )
    {
        printf(" Unable to determine cuda device count, error is %d, count is %d\n",
            cE, cudaDeviceCount );
        exit(-1);
    }
    if( (cE = cudaSetDevice( myrank % cudaDeviceCount )) != cudaSuccess )
    {
        printf(" Unable to have rank %d set to cuda device %d, error is %d \n",
            myrank, (myrank % cudaDeviceCount), cE);
        exit(-1);
    }

    if( (cE = cudaGetDevice( &assignedCudaDevice )) != cudaSuccess )
    {
        printf(" Unable to have rank %d set to cuda device %d, error is %d \n", 
            myrank, (myrank % cudaDeviceCount), cE);
        exit(-1);
    }

    if( assignedCudaDevice != (myrank % cudaDeviceCount) )
    {
        printf("MPI Rank %d: assignedCudaDevice %d NOT EQ to (myrank(%d) mod cudaDeviceCount(%d)) \n",
            myrank, assignedCudaDevice, myrank, cudaDeviceCount );
        exit(-1);
    }
    printf("Rank %d: My Cuda Device is %d\n", myrank, assignedCudaDevice);
}

void freeCudaMemory(Point * pts) {
    cudaFree(pts);
}

void c_callKernel(int numBlocks, int threadsCount, Point * p_data, size_t numPoints, int seed,
        float leftX, float lowerY, float rightX, float upperY) {
    data_kernel<<<numBlocks, threadsCount>>>(p_data, numPoints, seed,
        leftX, lowerY, rightX, upperY);
    cudaDeviceSynchronize();
}

void c_callCircleKernel(int numBlocks, int threadsCount, Point * p_data, size_t numPoints, int seed,
        float radius) {
    data_circle_kernel<<<numBlocks, threadsCount>>>(p_data, numPoints, seed,
        radius);
    cudaDeviceSynchronize();
}
