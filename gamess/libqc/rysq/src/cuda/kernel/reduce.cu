#include <algorithm>

#include "externals/cuda/types.hpp"
#include "externals/cuda/assert.h"
#include "externals/cxx/utility.hpp"

#include "cuda/matrix.hpp"

template<typename T>
__global__
static void kernel(const rysq::matrix_data_array<T> A, T scale, rysq::matrix<T> B) {
    T r = 0;
    int i  = threadIdx.x + blockIdx.x* blockDim.x;
    int j = blockIdx.y;
    int ij = A. layout().element_at(i,j);
    if (i <  A.layout().size1) {
	for (int k = 0; k < A.size(); ++k) {
	    r += A[k][ij];
	}  
	B(i,j) += scale*r;
    }
}


template<typename T>
void reduce(const rysq::matrix_data_array<T> A, T scale, rysq::matrix<T> B){
    size_t size1 = A.layout().size1, size2 = A.layout().size2;
    dim3 block(std::min<size_t>(size1, 128), 1, 1);
    dim3 grid(cxx::utility::qceiling<size_t>(size1, block.x), size2, 1);
    // std::cout <<  block << grid << " " << size1 << std::endl;
    // std:: cout << block << grid << A.size() << std::endl;
    ::kernel<<< grid, block >>>(A, scale, B); 

    cudaThreadSynchronize();
    cuda_assert();
}

void rysq::cuda::reduce(const rysq::matrix_data_array<double> A,
			double scale, rysq::matrix<double> B) {
    ::reduce<double>(A, scale, B);
}

