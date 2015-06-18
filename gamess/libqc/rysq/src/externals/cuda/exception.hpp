#ifndef CUDA_EXCEPTION_HPP_
#define CUDA_EXCEPTION_HPP_

#include <string>
#include <stdexcept>
#include <cuda.h>

namespace cuda {

    struct runtime_error : std::runtime_error {
	runtime_error(const std:: string &message)
	    : std::runtime_error(message) {}
    };

    struct configuration_error : cuda::runtime_error {
	configuration_error(const std::string &message)
	    : cuda::runtime_error(message) {}
    };

    void check_status() {
	cudaError_t status = cudaGetLastError();
	if (status == cudaSuccess) return;

	std::string message = cudaGetErrorString(status);
	if (status == cudaErrorInvalidConfiguration)
	    throw configuration_error(message);
	else
	    throw runtime_error(message);
    }

}

#endif /* _CUDA_EXCEPTION_HPP_ */
