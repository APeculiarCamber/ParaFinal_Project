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
typedef struct Point Point;

// LOCAL
__global__ void data_kernel(Point* p_data, size_t numPoints, int seed,
    float leftX, float lowerY, float rightX, float upperY) {

    curandState_t state; // TODO : better setting of seed
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


// Prevent name mangling
extern "C" {
    void c_cudaAlloc(Point ** points, size_t numPoints);

    void c_callKernel(int numBlocks, int threadsCount, Point * p_data, size_t numPoints, int seed,
        float leftX, float lowerY, float rightX, float upperY);
}

void c_cudaAlloc(Point ** points, size_t numPoints) {
    cudaMallocManaged(points, numPoints * sizeof(Point));
}

void c_callKernel(int numBlocks, int threadsCount, Point * p_data, size_t numPoints, int seed,
        float leftX, float lowerY, float rightX, float upperY) {
    data_kernel<<<numBlocks, threadsCount>>>(p_data, numPoints, seed,
        leftX, lowerY, rightX, upperY);
    cudaDeviceSynchronize();
}