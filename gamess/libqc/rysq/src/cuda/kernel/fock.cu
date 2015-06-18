#include <algorithm>

#include <boost/preprocessor/tuple/rem.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/utility/enable_if.hpp>
#include "boost/mpl/bool.hpp"

#include "cuda/kernel/quadrature.hpp"

#include "externals/cxx/sugar/quartet.hpp"

#include "externals/cuda/device/device.hpp"
#include "externals/cuda/assert.h"

#include <rysq/core.hpp>
#include "cuda/detail.hpp"
#include "cuda/kernel/quartet.hpp"



#include "externals/cxx/namespace.hpp"

template<size_t N, size_t I = 0>
struct repeat {
    template<class F, typename T0, typename T1, typename T2, typename T3>
    static void apply(F f, T0 &t0, T1 &t1, T2 &t2, T3 &t3) {
	f(t0, t1, t2, t3, I);
	repeat<N-1,I+1>::apply(f, t0, t1, t2, t3);
    }
};

template<size_t I>
struct repeat<0,I> {
    template<class F, typename T0, typename T1, typename T2, typename T3>
    static void apply(F f, T0 &t0, T1 &t1, T2 &t2, T3 &t3) {}
};

template<size_t A0, size_t A1, size_t A2, size_t A3>
struct mask {
    template<size_t B0, size_t B1, size_t B2, size_t B3>
    struct compare {
	static const bool value = (A0 == B0 && A1 == B1 && A2 == B2 && A3 == B3);
    };
};						   


template<size_t begin, size_t end>
struct range {
    template<typename T, size_t N>
    __host__ __device__
    static T multiply(const T (&vector)[N]) {
	return vector[begin]*range<begin+1,end>::multiply(vector);
    }
};

template<size_t begin>
struct range<begin, begin> {
    template<typename T, size_t N>
    __host__ __device__
    static T multiply(const T (&vector)[N]) { return 1; }
};

BEGIN_NAMESPACE(rysq, cuda, kernel)

namespace fock {

    
    template<size_t _1,size_t _2, bool consecutive = (_1 == _2 - 1)>
    struct integral_index_ {
	template<typename T, typename U>
	__device__
	static T eval(const T (&N)[4], const U &index) {
	    // int shift;
	    // if (N[_1] == 1) shift = (1<<17)/1;
	    // else if (N[_1] == 3) shift = (1<<17)/3;
	    // else if (N[_1] == 6) shift = (1<<17)/6;
	    // // else (N[_1] == 10) shift = (1<<17)/10;
	    // else shift = (1<<17)/10;
	    // T j = (index*shift)>>17;
	    // j += ((j + 1)*N[_1]-1 < index);

	    T j = index/N[_1];
	    return ((index - j*N[_1])*range<0,_1>::multiply(N) +
		    j*range<0,_2>::multiply(N));
	}
    };

    template<size_t _1,size_t _2>
    struct integral_index_<_1, _2, true> {
	template<typename T, typename U>
	__device__
	static T eval(const T (&N)[4], const U &index) {
	    return index*range<0,_1>::multiply(N);
	}
    };

    template<size_t _1,size_t _2, typename T, typename U>
    __device__
    T integral_index(const T (&N)[4], const U &index) {
	return integral_index_<_1,_2>::eval(N, index);
    }


    template<size_t a, size_t b, size_t c, size_t d, class enable = void>
    struct transform {
	template< typename T, typename U>
	__device__
	static void apply(const T (&N)[4], const double *D, double *F,
			  double scale, const double *Q,
			  const U &rank, double *shmem = NULL) {

	    __syncthreads();
	    if (rank < N[c]*N[d]) shmem[rank] = D[rank];
	    __syncthreads();

	    if (rank >= N[a]*N[b]) return;
	    double f = 0;
	    Q += integral_index<a,b>(N, rank);
	    for (T j = 0; j < N[d]; ++j) {
		for (T i = 0; i < N[c]; ++i) {
		    f += (*shmem)*(*Q);
		    ++shmem;
		    Q += range<0,c>::multiply(N);
		}
		Q += range<0,d>::multiply(N) - N[c]*range<0,c>::multiply(N);
	    }
	    F[rank] = scale*f;
	}
    };

    // template<>
    // struct rysq::cuda::kernel::fock::transform<1,2,0,3> {
    //     template<typename T, typename U>
    //     __device__
    //     static void apply(const T (&N)[4], const double *D, double *F,
    // 		      double scale, const double *Q,
    // 		      const U &rank, double *shmem = NULL) {
    // 	if (rank >= N[1]*N[2]) return;
    // 	double f = 0;
    // 	Q += integral_index<1,2>(N, rank);//rank*N[0];
    // 	for (int j = 0; j < N[3]; ++j) {
    // 	    for (int i = 0; i < N[0]; ++i) {
    // 		f += (*D)*Q[i];
    // 		++D;
    // 	    }
    // 	    Q += N[0]*N[1]*N[2];
    // 	}
    // 	F[rank] = scale*f;
    //     }
    // };


    template<>
    struct transform<0,3,1,2> {
	template<typename T, typename U>
	__device__
	static void apply(const T (&N)[4], const double *D, double *F,
			  double scale, const double *Q,
			  const U &rank, double *shmem = NULL) {

	    __syncthreads();
	    if (rank < N[1]*N[2]) shmem[rank] = D[rank];
	    __syncthreads();
	    
	    if (rank >= N[0]*N[3]) return;
	    double f = 0;
	    Q += integral_index<0,3>(N, rank);
	    for (T i = 0; i < N[1]*N[2]; ++i) {
	     	f += shmem[i]*(*Q);
	     	Q += N[0];
	    }
	    F[rank] = scale*f;
	}
    };

    template<size_t a, size_t b, size_t c, size_t d>
    struct transform <a,b,c,d, 
		      typename boost::enable_if_c<
			  (mask<a,b,c,d>::template compare<0,1,2,3>::value ||
			   mask<a,b,c,d>::template compare<2,3,0,1>::value)
			  >::type> {
	template<typename T, typename U>
	__device__
	static void apply(const T (&N)[4], const double *D, double *F,
			  double scale, const double *Q,
			  const U &rank, double *shmem = NULL) {
	    const bool transposed = (mask<a,b,c,d>::template compare<2,3,0,1>::value);

	    __syncthreads();
	    if (rank < N[c]*N[d]) shmem[rank] = D[rank];
	    __syncthreads();

	    if (rank >= N[a]* N[b]) return;
	    double f = 0;
	    Q += rank*((transposed) ? N[c]*N[d] : 1);
	    for (T ab = 0; ab < N[c]*N[d]; ++ab) {
		f += shmem[ab]*(*Q);
		Q += ((transposed) ? 1 :  N[a]* N[b]);
	    }
	    F[rank] = scale*f;
	}
    };



    template<class density_type, class fock_type>
    struct Transform {
	const density_type D;
	fock_type F;
	Transform(const density_type D, fock_type F) : D(D), F(F) {}
	template<typename T, size_t M>
	__device__
	void operator()(const T (&idx)[4], const double (&Q)[M], size_t size,
			const thread_block &block, double *shmem) {
	    apply(idx, Q, size, D, F, block, shmem);
	}
	template<typename T, size_t M>
	__device__
	static void apply(const T (&idx)[4], const double (&Q)[M], size_t size,
			  const density_type D, fock_type F,
			  const thread_block &block, double *shmem) {
	    __syncthreads();
	    eri::Transform::apply(Q, size, shmem, block);
	    __syncthreads();
	     apply(idx, F.blocks(), D, F, shmem, block.rank, shmem + size);
	}
	template<typename T0, typename T1, typename T2>
	__device__
	static void apply(const T0 (&idx)[4], const T1 (&N)[4],
			  const density_type D, fock_type F, const double *eri, 
			  const T2 &rank, double *shmem) {
#define APPLY(D, F) apply(N, (D), (F), scale, eri, rank, shmem)
	    double scale = 1.0/(rysq::Quartet<rysq::Shell>::symmetry(QUARTET(idx)));
	    transform<0,1,2,3>::APPLY(D.tile(1, idx), F.tile(0, idx));
	    transform<2,3,0,1>::APPLY(D.tile(0, idx), F.tile(1, idx));
	    transform<0,2,1,3>::APPLY(D.tile(5, idx), F.tile(2, idx));
	    transform<0,3,1,2>::APPLY(D.tile(4, idx), F.tile(3, idx));
	    transform<1,2,0,3>::APPLY(D.tile(3, idx), F.tile(4, idx));
	    transform<1,3,0,2>::APPLY(D.tile(2, idx), F.tile(5, idx));
#undef APPLY
	}
    };

    // template< class Transform, class density_type, class fock_type>
    // __global__
    // static void kernel(const rysq::Int4 *quartets,
    // 		       const density_type D, fock_type F,
    // 		       const double * eri) {
    // 	// using namespace ::cuda::device;
    // 	// size_t size = range<0,4>::multiply(F.blocks());
    // 	// // copy(eri + grid::rank()*size, shmem, size,  Block());
    // 	// const ushort rank =  block::rank();
    // 	// __shared__ rysq::Int4 index;
    // 	// if (rank == 0) index = quartets[grid::rank()];
    // 	// __syncthreads();
    // 	// // Transform:: apply(index.elems,eri, size, NULL, Block (), shmem + size, D, F);
    //    }
 
}



END_NAMESPACE(rysq, cuda, kernel)

BEGIN_NAMESPACE(rysq, cuda)

#define IMPL_TYPE kernel::Eri<kernel::fock::Transform<density_set,mapped_fock_set> >

void detail::Fock::operator()(const detail::Centers &centers,
			      const detail::Quartets &quartets,
			      const density_set D, mapped_fock_set F,
			      const Parameters &p) {
    ushort threads = 0;
    for (int i = 0; i < 6; ++i) {
	threads = std::max<size_t>(threads, F.tile_size(i));
    }
    // std::cout <<  shared<<std::endl;
    size_t shared = (threads + range<0,4>::multiply(F.blocks()))*sizeof(double);

    typedef kernel::fock::Transform<density_set,mapped_fock_set> T;
    IMPL_TYPE *impl = static_cast<IMPL_TYPE*>(impl_);
    (*impl)(centers, quartets,  p, T( D, F), threads, shared);

    cudaThreadSynchronize();
    cuda_assert( );
}


detail::Fock::Fock(const rysq::Quartet<rysq::Shell> &quartet,
		   const Transpose &transpose) {
    //std::cout << quartet << std::endl;
    this->impl_ = IMPL_TYPE::new_(quartet, transpose);
    if (!this->impl_) throw std::exception();
}

detail::Fock::~Fock() {
    delete static_cast<IMPL_TYPE*>(this->impl_);
}

END_NAMESPACE(rysq, cuda)
