#ifndef _CURYSQ_UTIL_H
#define _CURYSQ_UTIL_H

/**
   @file
*/

#include "vec.h"

#include <cuda.h>
#include <cuda_runtime.h>

/**
   @brief Asserts successful completion of CUDA function
*/
#define assertCudaSuccess {						\
	cudaError_t err = cudaGetLastError();				\
	if(err != cudaSuccess) {error(cudaGetErrorString(err));}	\
    }

#ifdef __CUDACC__

/** 
 * Return CTA rank, i.e. thread index
 * 
 * 
 * @return thread index
 */
__device__ short ctaRank() {
    return (threadIdx.x + blockDim.x*(threadIdx.y + blockDim.y*threadIdx.z));
}

/** 
 * Return CTA size, i.e. total number of threads
 * 
 * @return CTA size
 */
__device__ short ctaSize() {
    return (blockDim.x*blockDim.y*blockDim.z);
}

template<class T> __device__ void swap(T &a, T &b) {
    T q = a;
    a = b;
    b = q;
}

/** 
 * Reduce local variables in a warp
 * the reduced value is stored in shmem[warp*16]
 * 
 * @param r local variable
 * @param[out] shmem shared memory
 */
template<class T> __device__ void warpReduce(T &r, T *shmem) {
    short rank = ctaRank();
    short size = ctaSize();

    // sum within the warp using tree 
    if(rank % 2 == 0) shmem[rank/2] = r;
    if(rank % 2 != 0) shmem[rank/2] += r;

    if(rank % 4 == 0 && rank+2 < size) {
	shmem[rank/2] += shmem[rank/2+1];
    }
    if(rank % 8 == 0 && rank+4 < size) {
	shmem[rank/2] += shmem[rank/2+2];
    }
    if(rank % 16 == 0 && rank+8 < size) {
	shmem[rank/2] += shmem[rank/2+4];
    }
    if(rank % 32 == 0 && rank+16 < size) {
	shmem[rank/2] += shmem[rank/2+8];
    }

}


/** 
 * Reduce local variable in CTA
 * 
 * @param r local value
 * @param shmem shared memory
 * 
 * @return reduced value for rank 0, zero value otherwise
 */
template<class T> __device__ T ctaReduce(T &r, T *shmem){
    // first reduce without sync
    warpReduce(r, shmem);

    __syncthreads();

    // rank 0 adds across warps
    T q = T(0);

    if (ctaRank() == 0) {
	int n = ctaSize()/32 + (ctaSize()% 32 != 0);
	q = Vec_add(n, shmem, 16);
	__syncthreads();
    }

    return q;
}

#endif // __CUDACC__

#endif // _CURYSQ_UTIL_H

