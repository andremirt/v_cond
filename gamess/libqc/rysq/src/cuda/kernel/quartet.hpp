#ifndef _RYSQ_CUDA_KERNEL_QUARTET_HPP_
#define _RYSQ_CUDA_KERNEL_QUARTET_HPP_

#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include "cuda/detail.hpp"
#include "externals/cuda/device.hpp"
#include "externals/cuda/types.hpp"


namespace rysq {
    namespace cuda {


	namespace kernel {

	    typedef ::cuda::device::Block thread_block;
	    struct State;
	    struct Bra;
	    struct Ket;
	    struct Quartet;

	}
    }
}

namespace rysq {

    struct cuda::kernel::State {
	double r[2][3], dr[3], dr2;
	ushort K;
	// State(const ushort &K) : K(K) {}
	// static size_t size(ushort K) { return K*9; }
	// size_t size()const { return K*9; }
    protected:
	__host__ __device__
	size_t data(double *data,
		    double* &e, double* &A, double* &A1,
		    double* &rA, double* &rAi) {
	    e = data;
	    A = data + K;
	    A1 = data + 2*K;
	    rA = data + 3*K;
	    rAi = data + 6*K;
	    return 9*K;
	}
	__device__
	void primitives(const double2 *p,
			double *e, double *A, double *A1,
			double *rA, double *rAi,
			const thread_block &block) {
	    // compute bra values
	    for (ushort Kij = block.rank; Kij < this->K; Kij += block.size) {
		double2 aij = p[Kij];
		double a = aij.x + aij.y;
		double a1 = 1.0/a;
		// printf("%f\n", aij.x + aij.y);
		A[Kij] = a;
		A1[Kij] = a1;
		e[Kij] = exp(-aij.x*aij.y*a1*this->dr2);
		for (ushort i = 0; i < 3; ++i) {
		    rA[i+Kij*3] = a1*(aij.x*this->r[0][i] + aij.y*this->r[1][i]);
		    rAi[i+Kij*3] = rA[i+Kij*3] - this->r[0][i];
		}
	    }
	}
    };


    struct cuda::kernel::Bra : State {
	ushort mi, mj, m;
	double *e, *rA;
	double *A, *rAi;
	double *A1;
	__host__ __device__ ushort mij() const { return mi*mj; }
	__host__ __device__ size_t data(double *data) {
	    return State::data(data, e, A, A1, rA, rAi);
	}
	__device__ void primitives(const double2 *p,
				   const thread_block &block) {
	    State::primitives(p, e, A, A1, rA, rAi, block);
	}
    };


    struct cuda::kernel::Ket : State {
	ushort nk, nl, n;
	double *e, *rB;
	double *B, *rBi;
	double *B1;
	__host__ __device__ ushort nkl() const { return nk*nl; }
	__host__ __device__ size_t data(double *data) {
	    return State::data(data, e, B, B1, rB, rBi);
	}
	__device__ void primitives(const double2 *p,
				   const thread_block &block) {
	    State::primitives(p, e, B, B1, rB, rBi, block);
	}
    };


    struct cuda::kernel::Quartet {
	Quartet(const rysq::Quartet<rysq::Shell> &quartet) : centers_() {
	    initialize(quartet);
	} 
	Quartet(const detail::Quartet &quartet,
		const detail::Centers centers,
		unsigned short permutation)
	    : centers_(centers)
	{
	    permutation_ = permutation;
	    initialize(quartet);
	}
	__device__ size_t size() const { return size_; }
	__device__ void centers(const int *index, State &bra, State &ket,
				const thread_block &block) {
	    // compute bra and ket distances
	    for (ushort i = block.rank; i < 3; i +=  block.size) {
		bra.r[0][i] = centers_[index[0]].elems[i];
		bra.r[1][i] = centers_[index[1]].elems[i];
		ket.r[0][i] = centers_[index[2]].elems[i];
		ket.r[1][i] = centers_[index[3]].elems[i];
		bra.dr[i] = bra.r[0][i] - bra.r[1][i];
		ket.dr[i] = ket.r[0][i] - ket.r[1][i];
	    }
	    if (block.rank == 0) {
		bra.dr2 = cxx::math::dot(bra.dr);
		ket.dr2 = cxx::math::dot(ket.dr);
	    }
	    __syncthreads();
	}
	__host__ __device__ void initialize(Bra &bra, Ket &ket) {
	    bra.K = K(0)*K(1);
	    bra.mi = L(0) + 1;
	    bra.mj = L(1) + 1;
	    bra.m = bra.mi + bra.mj - 1;
	    ket.K = K(2)*K(3);
	    ket.nk = L(2) + 1;
	    ket.nl = L(3) + 1;
	    ket.n = ket.nk + ket.nl - 1;
	}
	__host__ __device__ ushort K(int i) const { return ((K_ >> i*4) & 0xf); }
	__device__ size_t K() const { return K(0)*K(1)*K(2)*K(3); }
	__host__ __device__ ushort L(int i) const { return abs(type(i)); }
	// __host__ __device__ ushort L(int i) const { return ((L_ >> i*4) & 0xf); }
	__host__ __device__ rysq::type type(int i) const {
	    return rysq::type(((type_ >> i*4) & 0xf) - 1);
	}
	__host__ __device__ bool hybrid(int i) const { return type(i) == rysq::SP; }
	__host__ __device__ ushort size(int i) const {
	    return ((L(i)+2)*(L(i)+1))/2 + hybrid(i);
	}
	__host__ __device__ ushort permutation()const { return permutation_; }
    private:
	unsigned short permutation_;
	unsigned short size_;
	unsigned short K_, type_; //L_, type_;
	const detail::Centers centers_;
	template<class Q>
	void initialize(const Q &quartet) {
	    size_ = quartet.size();
	    // permutation_ = quartet.permutation();
	    // L_ = 0;
	    type_ = 0;
	    K_ = 0;
	    for (int i = 0; i < 4; ++i) {		        
		// L_ = L_ | (quartet[i].L << i*4);
		type_ = type_|((quartet[i].type+1) << i*4);
		K_ = K_ | (quartet[i].K << i*4);
	    }
	}
    };

}

#endif /* _RYSQ_CUDA_KERNEL_QUARTET_HPP_ */
