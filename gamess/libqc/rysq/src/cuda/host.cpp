#include <cuda.h>
#include <cuda_runtime.h>
#include <assert.h>

#include <rysq/cuda.hpp>

#include "externals/cuda/host/core.hpp"
#include "externals/cuda/host/core.cpp"
#include "externals/cuda/assert.h"

using namespace rysq::cuda;

#define cuda_throw_error( ) {						\
	cudaError_t status = cudaGetLastError();			\
	if (status != cudaSuccess)					\
	    throw std::runtime_error(cudaGetErrorString(status));	\
    }

std::vector<int> rysq::cuda::get_devices() {
    int count;
    cudaGetDeviceCount(&count);
    cuda_throw_error( );    

    std::vector<int> devices;    
    for (int i = 0; i < count; ++i) {
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, i);
	if ((prop.major == 1 && prop.minor >= 3) || 
	    (prop.major > 1)) {
	    // std::cout << "device " << i << std::endl;
	    devices.push_back(i);
	}
    }
    return devices;
}

void rysq::cuda::set_device(size_t device) {
    cudaSetDevice(device);
    cuda_throw_error( );
}

void rysq::cuda::thread_exit() {
    cudaThreadExit();
    cuda_throw_error( );
}

device_ptr<void> rysq::cuda::malloc(size_t size) {
    void * ptr;
    cudaMalloc(&ptr, size);
    cuda_throw_error( );
    return device_ptr<void>(ptr);
}

void rysq::cuda::free(device_ptr<void> ptr) {
    cudaFree(ptr.data());
    cuda_throw_error( );
}

void rysq::cuda::memcpy(device_ptr<const void> from, void *to, size_t size) {
    cudaMemcpy(to, from.data(), size, cudaMemcpyDeviceToHost);
    cuda_throw_error( );
}

void rysq::cuda::memcpy(const void *from, device_ptr<void> to, size_t size) {
    // std::cout  << "memcpy" << size<< std::endl;
    cudaMemcpy(to.data(), from, size, cudaMemcpyHostToDevice);
    cuda_throw_error( );

}

void rysq::cuda::memset(device_ptr<void> ptr, char byte, size_t size) {
    cudaMemset(ptr.data(), ptr, size);
    cuda_throw_error( );
}

