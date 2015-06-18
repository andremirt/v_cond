#include <boost/utility/binary.hpp>
#include <boost/utility/binary.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/first_n.hpp>
#include <boost/preprocessor/seq/enum.hpp>

#include "externals/cuda/exception.hpp"
#include "externals/cuda/types.hpp"
#include "externals/cuda/device/device.hpp"

#include "externals/cxx/utility.hpp"

#include "transpose.hpp"
#include "cuda/detail.hpp"
#include "cuda/kernel/quartet.hpp"
#include "cuda/kernel/recurrence.hpp"
#include "cuda/kernel/transfer.hpp"



#include "roots/rysq_roots.h"

#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/for_each.hpp>

// using namespace ::cuda;
// using namespace rysq::cuda;


#define Q0000 BOOST_BINARY(0000)
#define Q1000 BOOST_BINARY(0001)
#define Q1010 BOOST_BINARY(0101)
#define Q1100 BOOST_BINARY(0011)
#define Q1110 BOOST_BINARY(0111)
#define Q1011 BOOST_BINARY(1101)
#define Q1111 BOOST_BINARY(1111)

#define Q0101 BOOST_BINARY(1010)
#define Q0011 BOOST_BINARY(1100)
#define Q0001 BOOST_BINARY(1000)
#define Q0100 BOOST_BINARY(0010)

#include "externals/cxx/namespace.hpp"

BEGIN_NAMESPACE(rysq, cuda, kernel)


template<class Transform>
struct Eri {
    static Eri* new_(const rysq::Quartet<rysq::Shell> &quartet,
		     const rysq::Transpose &transpose);
    virtual ~Eri() {}
    virtual void operator()(const detail::Centers &centers,
			    const detail::Quartets &quartets,
			    const Parameters &p,
			    Transform transform,
			    ushort threads = 0, size_t shared = 0) = 0;
};

namespace eri {

    extern __shared__ double dmem[];

    struct configuration;

    struct Transform;

    template<int mask, size_t N, class Transform, size_t M = 0>
    struct Kernel;

    template<ushort M, int N, int mask, class Transform>
    __global__ 
    static void kernel(kernel::Quartet quartet,
		       const rysq::Int4 *quartets,
		       float cutoff,
		       Transform transform);


}

END_NAMESPACE(rysq, cuda, kernel)

__constant__ double2 eri_quartet_data[36*36 + 4*36];
__constant__ ushort3 eri_quartet_index2d[2048];

BEGIN_NAMESPACE(rysq, cuda, kernel, eri)
template<ushort M, int N, int mask, class Transform>
__global__ 
static void kernel(kernel::Quartet quartet, const rysq::Int4 *quartets,
		   float cutoff, Transform transform) {

    // typedef Quartet<N,mask> Quartet;

    using namespace ::cuda::device;
    
    __shared__ kernel::Bra bra;
    __shared__ kernel::Ket ket;

    const Block block;

    __shared__ ushort stride;

    __shared__ ptr3<double> I;
    __shared__ double *C, *Cps[4], *tmp;

    __shared__ double *t2W;

    __shared__  rysq::Int4 index;
    {

	if (block.rank == 0)  {
	    quartet.initialize(bra,	ket);

	    ushort size2d = bra.mij()*ket.nkl();
	    stride = size2d*3;

	    size_t offset = 0;
	    offset += bra.data(dmem);
	    // printf("%u\n",  (int)offset);
	    offset += ket.data(dmem + offset);
	    // printf("%u\n",  (int)offset);

	    C = dmem + offset;
	    offset += quartet.K();

	    for (int i = 0; i < 4; ++i) {
		Cps[i] = dmem + offset;
		offset += ((i < 2) ? bra.K : ket.K)*quartet.hybrid(i);
	    }

	    // printf("%u\n",  (int)offset);
	    t2W = dmem + offset;
	    offset += 2*N*quartet.K();
	    // printf("%u\n",  (int)offset);
	    I.x = dmem + offset + 0*size2d;
	    I.y = dmem + offset + 1*size2d;
	    I.z = dmem + offset + 2*size2d;
	    tmp = dmem + offset + size2d*3*N;
	    // printf("%p %i %p %p %p %p\n", dmem, size2d, I.x,I.y,I.z,tmp);
	
	    // index = quartets.get(grid::rank(), dmem);
	    index = quartets[grid::rank()];
	}
	int *pindex = (int*)dmem;
	if (block.rank < 4) {
	    pindex[block.rank] = 
		cxx::utility::permute(block.rank, index.elems, quartet.permutation());
	}
	quartet.centers(pindex, bra, ket, block);
    }

    __syncthreads();

    {
    	double2 *cAij = eri_quartet_data;
    	double2 *cAkl = cAij + bra.K;
    	bra.primitives(cAij, block);
    	ket.primitives(cAkl, block);
    }

    double2 *cAB = eri_quartet_data + (bra.K + ket.K);

    {
	ushort size = 0;
#pragma unroll
	for (int i = 0; i < 4; ++i) {
	    size += quartet.hybrid(i)*quartet.K(i);
	}
	double *tmp = t2W; // use as temporary
	copy((const double*)(cAB + quartet.K()), tmp, size, block);
	__syncthreads();
	
#define COEFFICIENTS(index)						\
	if (quartet.hybrid(index)) {					\
	    ushort K1 = quartet.K(index + !(index%2));			\
	    for (ushort j =  block.rank; j < K1;  j += block.size) {	\
		char K0 = quartet.K(index - (index%2));			\
		for (char i = 0; i < K0;  ++i){				\
		    Cps[index][i+j*K0] = tmp[((index%2 == 0) ? i : j)];	\
		}							\
	    }								\
	    tmp += quartet.K(index);					\
	}
	
	COEFFICIENTS(0);
	COEFFICIENTS(1);
	COEFFICIENTS(2);
	COEFFICIENTS(3);
#undef COEFFICIENTS
    }

    // known at compile time, correct up to eight roots
    const ushort N2 = (N < 5) ? N + N%2 : 8;

    __syncthreads();

    // partition block into N2 parts
    for (ushort K = block.rank/N2; K < bra.K*ket.K; K += block.size/N2) {
    	ushort Kkl = K/bra.K;
    	ushort Kij = K - Kkl*bra.K;

    	ushort pRank = block.rank%N2;
    	double2 AB = cAB[K];
    	if (pRank == 0) C[K] = AB.x*bra.e[Kij]*ket.e[Kkl];
    	//if(fabs(C[Kij]) < cutoff) continue;
    	double AB1 = AB.y;

    	double *t2 = t2W + 2*N*K;
    	double *W = t2 + N;

    	double rho = (bra.A[Kij])*(ket.B[Kkl])*AB1;
    	double X = rho*distance2(&bra.rA[Kij*3], &ket.rB[Kkl*3]);
    	rysq::roots<N>(X, t2, W, pRank);
    	if(pRank < N) t2[pRank] *= AB1;
    }


    ushort3 index2d[M];
    double Q[M];
    unsigned int hybrid = 0;
#pragma unroll
    for (int i = 0; i < M; ++i) {
    	ushort j = block.rank + i*block.size;
    	index2d[i] = (j < quartet.size()) ?
	    (eri_quartet_index2d[j]) : (make_ushort3(0,0,0));
    	Q[i] = 0.0;
	ushort size = 1;
#pragma unroll
	for (int k = 0; k < 4; ++k) {
	    if (quartet.hybrid(k)) hybrid |= ((j/size)%4 == 0) << (i*4+k);
	    size *=quartet.size(k);
	    // j /= quartet.size(k);
	}
    }

    // printf("my rank %i %i\n", block.rank, block.size);
    // if (block.rank == 0) {
    // 	for (int i = 0; i < quartet.size(); ++i) {
    // 	    // printf("xxx %i", (int)quartet.size());
    // 	    // cuda::println(quartet_index2d[i]);
    // 	}
    // }

    __syncthreads();

    // ket contractions
    for (ushort Kkl = 0; Kkl < ket.K; ++Kkl) {
    	for (ushort Kij = 0; Kij <  bra.K; ++Kij) {
    	    ushort K = Kij + Kkl*bra.K;

	    if (fabs(C[K]) < cutoff/10000) continue;

	    double *t2 = t2W + 2*N*K;
	    double *W = t2 + N;

	    __shared__ double rAB[3];

	    if (block.rank < 3)
		rAB[block.rank] = bra.rA[block.rank+Kij*3] - ket.rB[block.rank+Kkl*3];

	    double *G = (!(mask & Q0101) or (mask == Q1111)) ? I.x : tmp;
	    const ushort n = (mask & Q0011) ? ket.n : 1;

	    __syncthreads();
	    recurrence<N>(bra.m, n,
			  bra.A[Kij], ket.B[Kkl], bra.A1[Kij], ket.B1[Kkl],
			  rAB, &bra.rAi[Kij*3], &ket.rBi[Kkl*3], t2, W, G,
			  block.size, block.rank);
	    __syncthreads();


	    if (!(mask & Q0011)) { // xx|ss
	    	transfer1<N>(bra.mi, bra.mj, 1, bra.dr, G, I.x);
	    }
	    else if (!(mask & Q0001)) { // xx|xs
	    	transfer1<N>(bra.mi, bra.mj, ket.nk, bra.dr, G, I.x);
	    }
	    else if (!(mask & Q0100)) { // xs|xx
	    	transfer2<N>(bra.mi, ket.nk, ket.nl, ket.dr, G, I.x);
	    }
	    else if (mask == Q1111) { // xx|xx
		// printf("%i %i %i %p %p\n", bra.mi, bra.mj, ket.n,  G, tmp);
	    	transfer1<N>(bra.mi, bra.mj, ket.n, bra.dr, G, tmp);
	    	__syncthreads();
	    	transfer2<N>(bra.mi*bra.mj, ket.nk, ket.nl, ket.dr, tmp, I.x);
	    }
	    __syncthreads();

#pragma unroll
	    for (int i = 0; i < M; ++i) {
		double c = C[K];
		if ((hybrid >> i*4) & 0x1) c *= Cps[0][Kij];
		if ((hybrid >> i*4) & 0x2) c *= Cps[1][Kij];
		if ((hybrid >> i*4) & 0x4) c *= Cps[2][Kkl];
		if ((hybrid >> i*4) & 0x8) c *= Cps[3][Kkl];

		const ushort3 &idx = index2d[i];
		double q = 0.0;
		// q = integral_sum<N>(I.x + idx.x, I.y + idx.y, I.z + idx.z, stride);
#pragma unroll
		for (ushort a = 0, off = 0; a < N; ++a, off += stride) {
		    q += I.x[idx.x+off]*I.y[idx.y+off]*I.z[idx.z+off];  
		}
		Q[i] += c*q;
	    }

	}
    }    

    transform(index.elems, Q, quartet.size(), block, dmem);

}


struct Transform {
    double *eri;
    explicit Transform(double *eri): eri(eri) {}
    template<typename T, size_t  M>
    __device__    
    void operator()(const T (&index)[4], const double (&Q)[M], size_t size,
		    const thread_block &block, double *shmem) {
	apply(Q, size, eri + ::cuda::device::grid::rank()*size, block);
    }
    template<size_t  M>
    __device__    
    static void apply(const double (&Q)[M], size_t size,
		      double *eri,
		      const thread_block &block) {
	eri += block.rank;
#pragma unroll
	for (int i = 0; i < M; ++i) {
	    if (block.rank + i*block.size < size) *eri = rysq::SQRT_4PI5*Q[i];
	    eri += block.size;
	}
    }
    __device__    
    static void apply(const double (&Q)[1], size_t size, double *eri,
		      const thread_block &block) {
	if (block.rank < size) eri[block.rank] = rysq::SQRT_4PI5*Q[0];
    }
};


template<int mask, size_t N, class Transform, size_t M_>
struct Kernel : Eri <Transform> {

    static const size_t M = M_;


    static size_t get_occupancy(const rysq::Quartet<rysq::Shell> &quartet) {
	try {
	    size_t threads = std::max<size_t>(multiply(block(quartet)), 32);
	    return std::min(16384/(threads*60), 16384/shared(quartet));
	}
	catch (std::exception &e) {
	    return 0;
	}	
    }

    static dim3 block(const rysq::Quartet<rysq::Shell> &quartet) {
	using namespace cxx::utility;
	size_t m = (quartet[0].L + quartet[1].L + 1);
	size_t n = (quartet[2].L + quartet[3].L + 1);
	dim3 block;
	block.x = ceiling2(max(m, n, N));
	block.y = 3;
	block.z = max(N, qceiling(size_t(quartet.size()), (M*block.x*block.y)));

	if (multiply(block)*60 > 16384)
	    throw std::range_error("to many registers");

	return block;
    }

    static size_t shared(const rysq::Quartet<rysq::Shell> &quartet) {
	kernel::Bra bra;
	kernel::Ket ket;
	(kernel::Quartet(quartet)).initialize(bra, ket);

	size_t	shared = bra.data(NULL) + ket.data(NULL);
	shared += quartet.K(); // coefficients

	// hybrid s/p coefficients
	shared += quartet[0].is_hybrid()*bra.K;
	shared += quartet[1].is_hybrid()*bra.K;
	shared += quartet[2].is_hybrid()*ket.K;
	shared += quartet[3].is_hybrid()*ket.K;
	
	shared += 2*N*quartet.K(); // roots and weights
	shared += bra.mij()*ket.nkl()*(3*N); // 2-D integrals

	// temporary transfer memory
	size_t tmp = 0;
	if (mask == Q1000 || mask == Q1010) tmp = 0;
	else if (mask == Q1111) tmp = (bra.mi*bra.mj)*(ket.n)*(3*N);
	else tmp = (bra.m)*(ket.n)*(3*N);
	shared += tmp;

	shared *= sizeof(double);
	// shared = 15000;

	if (shared > (16384-1024))
	    throw std::range_error("resources exceeded");

	return shared;
    }

    Kernel() {}

    Kernel(const rysq::Quartet<rysq::Shell> &quartet)
	:  block_(block(quartet)), shared_(shared(quartet)),
	   quartet_(quartet)
    {
	// std::cout <<  "new kernel " << M << std::endl;
    }

    void operator()(const detail::Centers &centers,
		    const detail::Quartets &quartets,
		    const Parameters &p,
		    Transform transform,
		    ushort threads, size_t shared) {

	copy(quartet_.data(), "eri_quartet_data");
	::cuda::check_status();

	copy(quartet_.index2d(), "eri_quartet_index2d");
	::cuda::check_status();

	ushort permutation = transpose_permutation(transpose_.value);
	kernel::Quartet quartet(quartet_, centers, permutation);

	using std::max;
	using cxx::utility::qceiling;

	dim3  block = block_;
	block.z = max(qceiling<typeof(block.z)>(threads, block.x*block.y), block.z);
	if (block.x* block.y* block.z < threads)
	    throw std::runtime_error("");

	shared = max(this->shared_, shared);

	// std::cout <<  quartet. permutation() << std::endl;

	for (int i = 0; i <  quartets.size(); i += 65535) {
	    dim3 grid = dim3(std::min<size_t>(quartets.size() - i, 65535));
	    // std::cout <<  grid << block << shared << std::endl;
	    eri::kernel<M, N, mask, Transform> <<< grid, block, shared >>>
		(quartet, quartets.data() + i, p.cutoff, transform);

	    try {
		::cuda::check_status();
		cudaThreadSynchronize();
		::cuda::check_status();
	    }
	    catch (::cuda::configuration_error&) {
		std::cout << __FILE__ << ":" << __LINE__
			  << " invalid configuration: "
			  << grid << " " << block << " " << shared << std::endl;
		throw;
	    }
	}

    }

private:
    dim3 block_;
    size_t shared_;
    detail::Quartet quartet_;
    Transpose transpose_;
};




struct configuration {

    struct data {
	size_t M;
	size_t occupancy;
	const rysq::Quartet<rysq::Shell> &quartet;
	data(const rysq::Quartet<rysq::Shell> &quartet)
	    : M(0), occupancy(0), quartet(quartet) {}
    };

    struct find {
	data &data_;
	explicit find(configuration::data &data) : data_(data) {}
	template<class K>
	void operator()(const K &kernel) {
	    int o = kernel.get_occupancy(data_.quartet);
	    if (o > data_.occupancy) {
		// std::cout << K::M << std::endl;
		data_.occupancy = o;
		data_.M = K::M;
	    }
	}
    };

    struct create {
	typedef void* pointer;
	const configuration::data &data;
	pointer &kernel;
	explicit create(const configuration::data &data, pointer &kernel)
	    : data(data), kernel(kernel) {}
	template<class K>
	void operator()(const K &kernel) {
	    if (data.M == K::M) {
		if (this->kernel)
		    throw std::runtime_error("runtime_error");
		this->kernel = new K(data.quartet);
	    }
	}
	template<class T>
	operator Eri<T>*() {
	    return static_cast<Eri<T>*>(this->kernel);
	}
    };

};

namespace mpl = boost::mpl;

template<int mask, size_t N, class Transform>
struct Kernel< mask, N, Transform, 0> {
    typedef mpl::vector_c<int,1,3,6> M;

    template<class M>
    struct kernel_type {
	typedef Kernel<mask, N, Transform,  M::value> type;
    };

    typedef typename mpl::transform<M, kernel_type<mpl::_1>
    				    >::type kernel_types;

static Eri<Transform>* new_(const rysq::Quartet<rysq::Shell> &quartet,
			    const rysq::Transpose &transpose) {
    // std::cout <<  "new_" << std::endl;
    configuration::data data(quartet);
    mpl::for_each<kernel_types>(configuration::find(data));
    void *kernel = NULL;
    mpl::for_each<kernel_types>(configuration::create(data, kernel));
    return static_cast<Eri<Transform>*>(kernel);
}

};

END_NAMESPACE(rysq, cuda, kernel, eri)

BEGIN_NAMESPACE(rysq, cuda, kernel)

template<class Transform>
kernel::Eri<Transform>*
kernel::Eri<Transform>::new_(const rysq::Quartet<rysq::Shell> &quartet,
			     const rysq::Transpose &transpose) {

    using eri::Kernel;

    //    if (quartet[3].is_hybrid() || quartet[3].is_hybrid() || quartet[3].is_hybrid())
    // return NULL;

    if(quartet.size() > 1600) return NULL;
    if(quartet.size() < 32) return NULL;

    int N = quartet.L()/2 + 1;
    int mask = cxx::utility::bitmask(quartet[0].L, quartet[1].L,
				     quartet[2].L, quartet[3].L);

#define ELEM BOOST_PP_SEQ_ELEM
#define ENUM BOOST_PP_SEQ_ENUM
#define FOR_EACH_PRODUCT BOOST_PP_SEQ_FOR_EACH_PRODUCT

#define SEQ_MASK	(Q1111)(Q1110)(Q1011)(Q1100)(Q1010)

#define KERNEL(r, params)						\
    if ((mask == ELEM(0, params)) && (N == ELEM(1, params))) {		\
	try {								\
	    return Kernel<ENUM(params), Transform>::new_(quartet, transpose); \
	}								\
	catch (std::exception&) { return NULL; }			\
    } else

#define SEQ_N		(2)(3)(4)(5)
    FOR_EACH_PRODUCT(KERNEL, (SEQ_MASK)(SEQ_N)) {
	return NULL;
    }
#undef SEQ_N

#undef KERNEL
#undef SEQ_MASK

#undef ELEM
#undef ENUM
#undef FOR_EACH_PRODUCT

}

END_NAMESPACE(rysq, cuda, kernel)
