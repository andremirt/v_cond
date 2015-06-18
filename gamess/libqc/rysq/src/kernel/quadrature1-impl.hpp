/**  
 @file 
 @warning Automatically Generated
*/
/**  
 @warning AUTOMATICALLY GENERATED
*/



#ifndef _RYSQ_KERNEL_QUADRATURE1_IMPL_HPP_
#define _RYSQ_KERNEL_QUADRATURE1_IMPL_HPP_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// #include <pmmintrin.h>

#ifdef RYSQ_WITH_SSE
#warning "Using SSE instructions"
#include <pmmintrin.h>
#include <xmmintrin.h>


#define D128 __m128d
#define ZERO _mm_setzero_pd()
#define SET1(v) _mm_set1_pd((v))

#define LOAD(m) _mm_load_pd((m))
#define LOADU(m) _mm_loadu_pd((m))
#define LOAD1(m) _mm_load_sd((m))
#define	LOADDUP(m) _mm_loaddup_pd((m))

#define STORE(m,r) _mm_store_pd((m), (r))
#define STOREU(m,r) _mm_storeu_pd((m), (r))
#define STORE1(m,r) _mm_store_sd((m), (r))

#define MUL(x,y) _mm_mul_pd((x), (y))
#define ADD(a,b) _mm_add_pd((a), (b))
#define HADD(a,b) _mm_hadd_pd((a), (b))

#define MUL1(x,y) _mm_mul_sd((x), (y))
#define ADD1(a,b) _mm_add_sd((a), (b))

#endif

#include <math.h>
#include "meta.hpp"
#include "cxx/sugar/dimension.hpp"


BEGIN_NAMESPACE(rysq, kernel, quadrature)


    //unrolled bras

#define Ix(a,i,j) DIMENSION(Ix, (NT,Li1,*), (a,i,j))
#define Iy(a,i,j) DIMENSION(Iy, (NT,Li1,*), (a,i,j))
#define Iz(a,i,j) DIMENSION(Iz, (NT,Li1,*), (a,i,j))

/** 
    @brief <ss| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::S,rysq::S> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::S,rysq::S> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 1;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
#else
     //double qK0 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x00, y00), z00));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x00, y00), z00));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD1(qK0, MUL1(C00, HADD(q0, q0)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,0,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	num += 1; //num += (fabs(I[0]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	}
	else {
	}
	STORE1(&I[0], ADD1(qK0, LOAD1(&I[0])));
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[0]*NORMALIZE[0]*qK0;
	// num += (fabs(I[0]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <ps| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::P,rysq::S> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::P,rysq::S> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 2;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x10, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x00, y10), z00));
	    q2 = ADD(q2, MUL(MUL(x00, y00), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x10, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x00, y10), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x00, y00), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD1(qK2, MUL1(C00, HADD(q2, q2)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,0,0);
	    q2 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	num += 3; //num += (fabs(I[2]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	}
	STORE1(&I[2], ADD1(qK2, LOAD1(&I[2])));
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[1]*NORMALIZE[0]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[2]*NORMALIZE[0]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[3]*NORMALIZE[0]*qK2;
	// num += (fabs(I[2]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <ds| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::D,rysq::S> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::D,rysq::S> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 3;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x20, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x00, y20), z00));
	    q2 = ADD(q2, MUL(MUL(x00, y00), z20));
	    q3 = ADD(q3, MUL(MUL(x10, y10), z00));
	    q4 = ADD(q4, MUL(MUL(x10, y00), z10));
	    q5 = ADD(q5, MUL(MUL(x00, y10), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x20, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x00, y20), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x00, y00), z20));
	    q3 = ADD1(q3, MUL1(MUL1(x10, y10), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x10, y00), z10));
	    q5 = ADD1(q5, MUL1(MUL1(x00, y10), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD(qK2, MUL(C00, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C00, HADD(q4, q5)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,2,0)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,0)*Iy(a,2,0)*Iz(a,0,0);
	    q2 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,2,0);
	    q3 += Ix(a,1,0)*Iy(a,1,0)*Iz(a,0,0);
	    q4 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,1,0);
	    q5 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+0];
	//I[3] += q3*C[k+0];
	I[3] += q3*C_[0];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+0];
	//I[5] += q5*C[k+0];
	I[5] += q5*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	num += 6; //num += (fabs(I[4]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[4]*NORMALIZE[0]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[5]*NORMALIZE[0]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[6]*NORMALIZE[0]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[7]*NORMALIZE[0]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[8]*NORMALIZE[0]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[9]*NORMALIZE[0]*qK5;
	// num += (fabs(I[5]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <fs| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::F,rysq::S> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::F,rysq::S> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 4;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x30 = LOAD(&Ix(a,3,0));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y30 = LOAD(&Iy(a,3,0));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z30 = LOAD(&Iz(a,3,0));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x30, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x00, y30), z00));
	    q2 = ADD(q2, MUL(MUL(x00, y00), z30));
	    q3 = ADD(q3, MUL(MUL(x20, y10), z00));
	    q4 = ADD(q4, MUL(MUL(x20, y00), z10));
	    q5 = ADD(q5, MUL(MUL(x10, y20), z00));
	    q6 = ADD(q6, MUL(MUL(x00, y20), z10));
	    q7 = ADD(q7, MUL(MUL(x10, y00), z20));
	    q8 = ADD(q8, MUL(MUL(x00, y10), z20));
	    q9 = ADD(q9, MUL(MUL(x10, y10), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x30 = LOAD1(&Ix(a,3,0));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y30 = LOAD1(&Iy(a,3,0));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z30 = LOAD1(&Iz(a,3,0));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x30, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x00, y30), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x00, y00), z30));
	    q3 = ADD1(q3, MUL1(MUL1(x20, y10), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x20, y00), z10));
	    q5 = ADD1(q5, MUL1(MUL1(x10, y20), z00));
	    q6 = ADD1(q6, MUL1(MUL1(x00, y20), z10));
	    q7 = ADD1(q7, MUL1(MUL1(x10, y00), z20));
	    q8 = ADD1(q8, MUL1(MUL1(x00, y10), z20));
	    q9 = ADD1(q9, MUL1(MUL1(x10, y10), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD(qK2, MUL(C00, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C00, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C00, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C00, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,3,0)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,0)*Iy(a,3,0)*Iz(a,0,0);
	    q2 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,3,0);
	    q3 += Ix(a,2,0)*Iy(a,1,0)*Iz(a,0,0);
	    q4 += Ix(a,2,0)*Iy(a,0,0)*Iz(a,1,0);
	    q5 += Ix(a,1,0)*Iy(a,2,0)*Iz(a,0,0);
	    q6 += Ix(a,0,0)*Iy(a,2,0)*Iz(a,1,0);
	    q7 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,2,0);
	    q8 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,2,0);
	    q9 += Ix(a,1,0)*Iy(a,1,0)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+0];
	//I[3] += q3*C[k+0];
	I[3] += q3*C_[0];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+0];
	//I[5] += q5*C[k+0];
	I[5] += q5*C_[0];
	//qK6 += q6*C[k+0];
	//I[6] += q6*C[k+0];
	I[6] += q6*C_[0];
	//qK7 += q7*C[k+0];
	//I[7] += q7*C[k+0];
	I[7] += q7*C_[0];
	//qK8 += q8*C[k+0];
	//I[8] += q8*C[k+0];
	I[8] += q8*C_[0];
	//qK9 += q9*C[k+0];
	//I[9] += q9*C[k+0];
	I[9] += q9*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[10]*NORMALIZE[0]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[11]*NORMALIZE[0]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[12]*NORMALIZE[0]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[13]*NORMALIZE[0]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[14]*NORMALIZE[0]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[15]*NORMALIZE[0]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[16]*NORMALIZE[0]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[17]*NORMALIZE[0]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[18]*NORMALIZE[0]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[19]*NORMALIZE[0]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <sps| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::SP,rysq::S> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::SP,rysq::S> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 2;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x00, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x10, y00), z00));
	    q2 = ADD(q2, MUL(MUL(x00, y10), z00));
	    q3 = ADD(q3, MUL(MUL(x00, y00), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x00, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x10, y00), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x00, y10), z00));
	    q3 = ADD1(q3, MUL1(MUL1(x00, y00), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C01 = LOADU(&C[k+0]);
	qK0 = ADD(qK0, MUL(C01, HADD(q0, q1)));
	D128 C11 = LOADDUP(&C[k+1]);
	qK2 = ADD(qK2, MUL(C11, HADD(q2, q3)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,0,0);
	    q2 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,0,0);
	    q3 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+1];
	//I[1] += q1*C[k+1];
	I[1] += q1*C_[1];
	//qK2 += q2*C[k+1];
	//I[2] += q2*C[k+1];
	I[2] += q2*C_[1];
	//qK3 += q3*C[k+1];
	//I[3] += q3*C[k+1];
	I[3] += q3*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	num += 4; //num += (fabs(I[2]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[0]*NORMALIZE[0]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[1]*NORMALIZE[0]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[2]*NORMALIZE[0]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[3]*NORMALIZE[0]*qK3;
	// num += (fabs(I[3]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <sp| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::S,rysq::P> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::S,rysq::P> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 1;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x01, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x00, y01), z00));
	    q2 = ADD(q2, MUL(MUL(x00, y00), z01));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x01, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x00, y01), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x00, y00), z01));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD1(qK2, MUL1(C00, HADD(q2, q2)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,0,0);
	    q2 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,0,1);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	num += 3; //num += (fabs(I[2]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	}
	STORE1(&I[2], ADD1(qK2, LOAD1(&I[2])));
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[0]*NORMALIZE[1]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[0]*NORMALIZE[2]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[0]*NORMALIZE[3]*qK2;
	// num += (fabs(I[2]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <pp| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::P,rysq::P> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::P,rysq::P> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 2;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    q0 = ADD(q0, MUL(MUL(x11, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x01, y10), z00));
	    q2 = ADD(q2, MUL(MUL(x01, y00), z10));
	    q3 = ADD(q3, MUL(MUL(x10, y01), z00));
	    q4 = ADD(q4, MUL(MUL(x00, y11), z00));
	    q5 = ADD(q5, MUL(MUL(x00, y01), z10));
	    q6 = ADD(q6, MUL(MUL(x10, y00), z01));
	    q7 = ADD(q7, MUL(MUL(x00, y10), z01));
	    q8 = ADD(q8, MUL(MUL(x00, y00), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    q0 = ADD1(q0, MUL1(MUL1(x11, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x01, y10), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x01, y00), z10));
	    q3 = ADD1(q3, MUL1(MUL1(x10, y01), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x00, y11), z00));
	    q5 = ADD1(q5, MUL1(MUL1(x00, y01), z10));
	    q6 = ADD1(q6, MUL1(MUL1(x10, y00), z01));
	    q7 = ADD1(q7, MUL1(MUL1(x00, y10), z01));
	    q8 = ADD1(q8, MUL1(MUL1(x00, y00), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD(qK2, MUL(C00, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C00, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C00, HADD(q6, q7)));
	qK8 = ADD1(qK8, MUL1(C00, HADD(q8, q8)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,0,0);
	    q2 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,1,0);
	    q3 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,0,0);
	    q4 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,0,0);
	    q5 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,1,0);
	    q6 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,0,1);
	    q7 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,0,1);
	    q8 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+0];
	//I[3] += q3*C[k+0];
	I[3] += q3*C_[0];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+0];
	//I[5] += q5*C[k+0];
	I[5] += q5*C_[0];
	//qK6 += q6*C[k+0];
	//I[6] += q6*C[k+0];
	I[6] += q6*C_[0];
	//qK7 += q7*C[k+0];
	//I[7] += q7*C[k+0];
	I[7] += q7*C_[0];
	//qK8 += q8*C[k+0];
	//I[8] += q8*C[k+0];
	I[8] += q8*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 9; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	}
	STORE1(&I[8], ADD1(qK8, LOAD1(&I[8])));
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[1]*NORMALIZE[1]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[2]*NORMALIZE[1]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[3]*NORMALIZE[1]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[1]*NORMALIZE[2]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[2]*NORMALIZE[2]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[3]*NORMALIZE[2]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[1]*NORMALIZE[3]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[2]*NORMALIZE[3]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[3]*NORMALIZE[3]*qK8;
	// num += (fabs(I[8]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <dp| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::D,rysq::P> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::D,rysq::P> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 3;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x21 = LOAD(&Ix(a,2,1));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y21 = LOAD(&Iy(a,2,1));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x21, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x01, y20), z00));
	    q2 = ADD(q2, MUL(MUL(x01, y00), z20));
	    q3 = ADD(q3, MUL(MUL(x11, y10), z00));
	    q4 = ADD(q4, MUL(MUL(x11, y00), z10));
	    q5 = ADD(q5, MUL(MUL(x01, y10), z10));
	    q6 = ADD(q6, MUL(MUL(x20, y01), z00));
	    q7 = ADD(q7, MUL(MUL(x00, y21), z00));
	    q8 = ADD(q8, MUL(MUL(x00, y01), z20));
	    q9 = ADD(q9, MUL(MUL(x10, y11), z00));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x21 = LOAD1(&Ix(a,2,1));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y21 = LOAD1(&Iy(a,2,1));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x21, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x01, y20), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x01, y00), z20));
	    q3 = ADD1(q3, MUL1(MUL1(x11, y10), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x11, y00), z10));
	    q5 = ADD1(q5, MUL1(MUL1(x01, y10), z10));
	    q6 = ADD1(q6, MUL1(MUL1(x20, y01), z00));
	    q7 = ADD1(q7, MUL1(MUL1(x00, y21), z00));
	    q8 = ADD1(q8, MUL1(MUL1(x00, y01), z20));
	    q9 = ADD1(q9, MUL1(MUL1(x10, y11), z00));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD(qK2, MUL(C00, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C00, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C00, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C00, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,2,1)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,1)*Iy(a,2,0)*Iz(a,0,0);
	    q2 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,2,0);
	    q3 += Ix(a,1,1)*Iy(a,1,0)*Iz(a,0,0);
	    q4 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,1,0);
	    q5 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,1,0);
	    q6 += Ix(a,2,0)*Iy(a,0,1)*Iz(a,0,0);
	    q7 += Ix(a,0,0)*Iy(a,2,1)*Iz(a,0,0);
	    q8 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,2,0);
	    q9 += Ix(a,1,0)*Iy(a,1,1)*Iz(a,0,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+0];
	//I[3] += q3*C[k+0];
	I[3] += q3*C_[0];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+0];
	//I[5] += q5*C[k+0];
	I[5] += q5*C_[0];
	//qK6 += q6*C[k+0];
	//I[6] += q6*C[k+0];
	I[6] += q6*C_[0];
	//qK7 += q7*C[k+0];
	//I[7] += q7*C[k+0];
	I[7] += q7*C_[0];
	//qK8 += q8*C[k+0];
	//I[8] += q8*C[k+0];
	I[8] += q8*C_[0];
	//qK9 += q9*C[k+0];
	//I[9] += q9*C[k+0];
	I[9] += q9*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[4]*NORMALIZE[1]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[5]*NORMALIZE[1]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[6]*NORMALIZE[1]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[7]*NORMALIZE[1]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[8]*NORMALIZE[1]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[9]*NORMALIZE[1]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[4]*NORMALIZE[2]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[5]*NORMALIZE[2]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[6]*NORMALIZE[2]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[7]*NORMALIZE[2]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK10 = ZERO;
     D128 qK12 = ZERO;
     D128 qK14 = ZERO;
     D128 qK16 = ZERO;
#else
     //double qK10 = 0.0;
     //double qK11 = 0.0;
     //double qK12 = 0.0;
     //double qK13 = 0.0;
     //double qK14 = 0.0;
     //double qK15 = 0.0;
     //double qK16 = 0.0;
     //double qK17 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q10 = ZERO;
	D128 q11 = ZERO;
	D128 q12 = ZERO;
	D128 q13 = ZERO;
	D128 q14 = ZERO;
	D128 q15 = ZERO;
	D128 q16 = ZERO;
	D128 q17 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    D128 z21 = LOAD(&Iz(a,2,1));
	    q10 = ADD(q10, MUL(MUL(x10, y01), z10));
	    q11 = ADD(q11, MUL(MUL(x00, y11), z10));
	    q12 = ADD(q12, MUL(MUL(x20, y00), z01));
	    q13 = ADD(q13, MUL(MUL(x00, y20), z01));
	    q14 = ADD(q14, MUL(MUL(x00, y00), z21));
	    q15 = ADD(q15, MUL(MUL(x10, y10), z01));
	    q16 = ADD(q16, MUL(MUL(x10, y00), z11));
	    q17 = ADD(q17, MUL(MUL(x00, y10), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    D128 z21 = LOAD1(&Iz(a,2,1));
	    q10 = ADD1(q10, MUL1(MUL1(x10, y01), z10));
	    q11 = ADD1(q11, MUL1(MUL1(x00, y11), z10));
	    q12 = ADD1(q12, MUL1(MUL1(x20, y00), z01));
	    q13 = ADD1(q13, MUL1(MUL1(x00, y20), z01));
	    q14 = ADD1(q14, MUL1(MUL1(x00, y00), z21));
	    q15 = ADD1(q15, MUL1(MUL1(x10, y10), z01));
	    q16 = ADD1(q16, MUL1(MUL1(x10, y00), z11));
	    q17 = ADD1(q17, MUL1(MUL1(x00, y10), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK10 = ADD(qK10, MUL(C00, HADD(q10, q11)));
	qK12 = ADD(qK12, MUL(C00, HADD(q12, q13)));
	qK14 = ADD(qK14, MUL(C00, HADD(q14, q15)));
	qK16 = ADD(qK16, MUL(C00, HADD(q16, q17)));

#else // SSE
	    
	// function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
	T q16 = 0.0;
	T q17 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,1,0);
	    q11 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,1,0);
	    q12 += Ix(a,2,0)*Iy(a,0,0)*Iz(a,0,1);
	    q13 += Ix(a,0,0)*Iy(a,2,0)*Iz(a,0,1);
	    q14 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,2,1);
	    q15 += Ix(a,1,0)*Iy(a,1,0)*Iz(a,0,1);
	    q16 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,1,1);
	    q17 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK10 += q10*C[k+0];
	//I[10] += q10*C[k+0];
	I[10] += q10*C_[0];
	//qK11 += q11*C[k+0];
	//I[11] += q11*C[k+0];
	I[11] += q11*C_[0];
	//qK12 += q12*C[k+0];
	//I[12] += q12*C[k+0];
	I[12] += q12*C_[0];
	//qK13 += q13*C[k+0];
	//I[13] += q13*C[k+0];
	I[13] += q13*C_[0];
	//qK14 += q14*C[k+0];
	//I[14] += q14*C[k+0];
	I[14] += q14*C_[0];
	//qK15 += q15*C[k+0];
	//I[15] += q15*C[k+0];
	I[15] += q15*C_[0];
	//qK16 += q16*C[k+0];
	//I[16] += q16*C[k+0];
	I[16] += q16*C_[0];
	//qK17 += q17*C[k+0];
	//I[17] += q17*C[k+0];
	I[17] += q17*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK10 = MUL(q, qK10);
	qK12 = MUL(q, qK12);
	qK14 = MUL(q, qK14);
	qK16 = MUL(q, qK16);
	num += 8; //num += (fabs(I[16]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[10]) & 0xF) {
	    // 10
	    STOREU(&I[10], ADD(qK10, LOADU(&I[10])));
	    STOREU(&I[12], ADD(qK12, LOADU(&I[12])));
	    STOREU(&I[14], ADD(qK14, LOADU(&I[14])));
	    STOREU(&I[16], ADD(qK16, LOADU(&I[16])));
	}
	else {
	    STORE(&I[10], ADD(qK10, LOADU(&I[10])));
	    STORE(&I[12], ADD(qK12, LOADU(&I[12])));
	    STORE(&I[14], ADD(qK14, LOADU(&I[14])));
	    STORE(&I[16], ADD(qK16, LOADU(&I[16])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[10] += scale*NORMALIZE[8]*NORMALIZE[2]*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*NORMALIZE[9]*NORMALIZE[2]*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*NORMALIZE[4]*NORMALIZE[3]*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*NORMALIZE[5]*NORMALIZE[3]*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*NORMALIZE[6]*NORMALIZE[3]*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*NORMALIZE[7]*NORMALIZE[3]*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*NORMALIZE[8]*NORMALIZE[3]*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*NORMALIZE[9]*NORMALIZE[3]*qK17;
	// num += (fabs(I[17]) >= tol);
    }
    else {
	// I[10] += scale*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*qK17;
	// num += (fabs(I[17]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <fp| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::F,rysq::P> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::F,rysq::P> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 4;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x31 = LOAD(&Ix(a,3,1));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 x21 = LOAD(&Ix(a,2,1));
	    D128 y30 = LOAD(&Iy(a,3,0));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z30 = LOAD(&Iz(a,3,0));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x31, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x01, y30), z00));
	    q2 = ADD(q2, MUL(MUL(x01, y00), z30));
	    q3 = ADD(q3, MUL(MUL(x21, y10), z00));
	    q4 = ADD(q4, MUL(MUL(x21, y00), z10));
	    q5 = ADD(q5, MUL(MUL(x11, y20), z00));
	    q6 = ADD(q6, MUL(MUL(x01, y20), z10));
	    q7 = ADD(q7, MUL(MUL(x11, y00), z20));
	    q8 = ADD(q8, MUL(MUL(x01, y10), z20));
	    q9 = ADD(q9, MUL(MUL(x11, y10), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x31 = LOAD1(&Ix(a,3,1));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 x21 = LOAD1(&Ix(a,2,1));
	    D128 y30 = LOAD1(&Iy(a,3,0));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z30 = LOAD1(&Iz(a,3,0));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x31, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x01, y30), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x01, y00), z30));
	    q3 = ADD1(q3, MUL1(MUL1(x21, y10), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x21, y00), z10));
	    q5 = ADD1(q5, MUL1(MUL1(x11, y20), z00));
	    q6 = ADD1(q6, MUL1(MUL1(x01, y20), z10));
	    q7 = ADD1(q7, MUL1(MUL1(x11, y00), z20));
	    q8 = ADD1(q8, MUL1(MUL1(x01, y10), z20));
	    q9 = ADD1(q9, MUL1(MUL1(x11, y10), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD(qK2, MUL(C00, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C00, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C00, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C00, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,3,1)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,1)*Iy(a,3,0)*Iz(a,0,0);
	    q2 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,3,0);
	    q3 += Ix(a,2,1)*Iy(a,1,0)*Iz(a,0,0);
	    q4 += Ix(a,2,1)*Iy(a,0,0)*Iz(a,1,0);
	    q5 += Ix(a,1,1)*Iy(a,2,0)*Iz(a,0,0);
	    q6 += Ix(a,0,1)*Iy(a,2,0)*Iz(a,1,0);
	    q7 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,2,0);
	    q8 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,2,0);
	    q9 += Ix(a,1,1)*Iy(a,1,0)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+0];
	//I[3] += q3*C[k+0];
	I[3] += q3*C_[0];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+0];
	//I[5] += q5*C[k+0];
	I[5] += q5*C_[0];
	//qK6 += q6*C[k+0];
	//I[6] += q6*C[k+0];
	I[6] += q6*C_[0];
	//qK7 += q7*C[k+0];
	//I[7] += q7*C[k+0];
	I[7] += q7*C_[0];
	//qK8 += q8*C[k+0];
	//I[8] += q8*C[k+0];
	I[8] += q8*C_[0];
	//qK9 += q9*C[k+0];
	//I[9] += q9*C[k+0];
	I[9] += q9*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[10]*NORMALIZE[1]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[11]*NORMALIZE[1]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[12]*NORMALIZE[1]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[13]*NORMALIZE[1]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[14]*NORMALIZE[1]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[15]*NORMALIZE[1]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[16]*NORMALIZE[1]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[17]*NORMALIZE[1]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[18]*NORMALIZE[1]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[19]*NORMALIZE[1]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK10 = ZERO;
     D128 qK12 = ZERO;
     D128 qK14 = ZERO;
     D128 qK16 = ZERO;
     D128 qK18 = ZERO;
#else
     //double qK10 = 0.0;
     //double qK11 = 0.0;
     //double qK12 = 0.0;
     //double qK13 = 0.0;
     //double qK14 = 0.0;
     //double qK15 = 0.0;
     //double qK16 = 0.0;
     //double qK17 = 0.0;
     //double qK18 = 0.0;
     //double qK19 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q10 = ZERO;
	D128 q11 = ZERO;
	D128 q12 = ZERO;
	D128 q13 = ZERO;
	D128 q14 = ZERO;
	D128 q15 = ZERO;
	D128 q16 = ZERO;
	D128 q17 = ZERO;
	D128 q18 = ZERO;
	D128 q19 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x30 = LOAD(&Ix(a,3,0));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y31 = LOAD(&Iy(a,3,1));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 y21 = LOAD(&Iy(a,2,1));
	    D128 z30 = LOAD(&Iz(a,3,0));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q10 = ADD(q10, MUL(MUL(x30, y01), z00));
	    q11 = ADD(q11, MUL(MUL(x00, y31), z00));
	    q12 = ADD(q12, MUL(MUL(x00, y01), z30));
	    q13 = ADD(q13, MUL(MUL(x20, y11), z00));
	    q14 = ADD(q14, MUL(MUL(x20, y01), z10));
	    q15 = ADD(q15, MUL(MUL(x10, y21), z00));
	    q16 = ADD(q16, MUL(MUL(x00, y21), z10));
	    q17 = ADD(q17, MUL(MUL(x10, y01), z20));
	    q18 = ADD(q18, MUL(MUL(x00, y11), z20));
	    q19 = ADD(q19, MUL(MUL(x10, y11), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x30 = LOAD1(&Ix(a,3,0));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y31 = LOAD1(&Iy(a,3,1));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 y21 = LOAD1(&Iy(a,2,1));
	    D128 z30 = LOAD1(&Iz(a,3,0));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q10 = ADD1(q10, MUL1(MUL1(x30, y01), z00));
	    q11 = ADD1(q11, MUL1(MUL1(x00, y31), z00));
	    q12 = ADD1(q12, MUL1(MUL1(x00, y01), z30));
	    q13 = ADD1(q13, MUL1(MUL1(x20, y11), z00));
	    q14 = ADD1(q14, MUL1(MUL1(x20, y01), z10));
	    q15 = ADD1(q15, MUL1(MUL1(x10, y21), z00));
	    q16 = ADD1(q16, MUL1(MUL1(x00, y21), z10));
	    q17 = ADD1(q17, MUL1(MUL1(x10, y01), z20));
	    q18 = ADD1(q18, MUL1(MUL1(x00, y11), z20));
	    q19 = ADD1(q19, MUL1(MUL1(x10, y11), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK10 = ADD(qK10, MUL(C00, HADD(q10, q11)));
	qK12 = ADD(qK12, MUL(C00, HADD(q12, q13)));
	qK14 = ADD(qK14, MUL(C00, HADD(q14, q15)));
	qK16 = ADD(qK16, MUL(C00, HADD(q16, q17)));
	qK18 = ADD(qK18, MUL(C00, HADD(q18, q19)));

#else // SSE
	    
	// function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
	T q16 = 0.0;
	T q17 = 0.0;
	T q18 = 0.0;
	T q19 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += Ix(a,3,0)*Iy(a,0,1)*Iz(a,0,0);
	    q11 += Ix(a,0,0)*Iy(a,3,1)*Iz(a,0,0);
	    q12 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,3,0);
	    q13 += Ix(a,2,0)*Iy(a,1,1)*Iz(a,0,0);
	    q14 += Ix(a,2,0)*Iy(a,0,1)*Iz(a,1,0);
	    q15 += Ix(a,1,0)*Iy(a,2,1)*Iz(a,0,0);
	    q16 += Ix(a,0,0)*Iy(a,2,1)*Iz(a,1,0);
	    q17 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,2,0);
	    q18 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,2,0);
	    q19 += Ix(a,1,0)*Iy(a,1,1)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK10 += q10*C[k+0];
	//I[10] += q10*C[k+0];
	I[10] += q10*C_[0];
	//qK11 += q11*C[k+0];
	//I[11] += q11*C[k+0];
	I[11] += q11*C_[0];
	//qK12 += q12*C[k+0];
	//I[12] += q12*C[k+0];
	I[12] += q12*C_[0];
	//qK13 += q13*C[k+0];
	//I[13] += q13*C[k+0];
	I[13] += q13*C_[0];
	//qK14 += q14*C[k+0];
	//I[14] += q14*C[k+0];
	I[14] += q14*C_[0];
	//qK15 += q15*C[k+0];
	//I[15] += q15*C[k+0];
	I[15] += q15*C_[0];
	//qK16 += q16*C[k+0];
	//I[16] += q16*C[k+0];
	I[16] += q16*C_[0];
	//qK17 += q17*C[k+0];
	//I[17] += q17*C[k+0];
	I[17] += q17*C_[0];
	//qK18 += q18*C[k+0];
	//I[18] += q18*C[k+0];
	I[18] += q18*C_[0];
	//qK19 += q19*C[k+0];
	//I[19] += q19*C[k+0];
	I[19] += q19*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK10 = MUL(q, qK10);
	qK12 = MUL(q, qK12);
	qK14 = MUL(q, qK14);
	qK16 = MUL(q, qK16);
	qK18 = MUL(q, qK18);
	num += 10; //num += (fabs(I[18]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[10]) & 0xF) {
	    // 10
	    STOREU(&I[10], ADD(qK10, LOADU(&I[10])));
	    STOREU(&I[12], ADD(qK12, LOADU(&I[12])));
	    STOREU(&I[14], ADD(qK14, LOADU(&I[14])));
	    STOREU(&I[16], ADD(qK16, LOADU(&I[16])));
	    STOREU(&I[18], ADD(qK18, LOADU(&I[18])));
	}
	else {
	    STORE(&I[10], ADD(qK10, LOADU(&I[10])));
	    STORE(&I[12], ADD(qK12, LOADU(&I[12])));
	    STORE(&I[14], ADD(qK14, LOADU(&I[14])));
	    STORE(&I[16], ADD(qK16, LOADU(&I[16])));
	    STORE(&I[18], ADD(qK18, LOADU(&I[18])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[10] += scale*NORMALIZE[10]*NORMALIZE[2]*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*NORMALIZE[11]*NORMALIZE[2]*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*NORMALIZE[12]*NORMALIZE[2]*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*NORMALIZE[13]*NORMALIZE[2]*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*NORMALIZE[14]*NORMALIZE[2]*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*NORMALIZE[15]*NORMALIZE[2]*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*NORMALIZE[16]*NORMALIZE[2]*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*NORMALIZE[17]*NORMALIZE[2]*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*NORMALIZE[18]*NORMALIZE[2]*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*NORMALIZE[19]*NORMALIZE[2]*qK19;
	// num += (fabs(I[19]) >= tol);
    }
    else {
	// I[10] += scale*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*qK19;
	// num += (fabs(I[19]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK20 = ZERO;
     D128 qK22 = ZERO;
     D128 qK24 = ZERO;
     D128 qK26 = ZERO;
     D128 qK28 = ZERO;
#else
     //double qK20 = 0.0;
     //double qK21 = 0.0;
     //double qK22 = 0.0;
     //double qK23 = 0.0;
     //double qK24 = 0.0;
     //double qK25 = 0.0;
     //double qK26 = 0.0;
     //double qK27 = 0.0;
     //double qK28 = 0.0;
     //double qK29 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q20 = ZERO;
	D128 q21 = ZERO;
	D128 q22 = ZERO;
	D128 q23 = ZERO;
	D128 q24 = ZERO;
	D128 q25 = ZERO;
	D128 q26 = ZERO;
	D128 q27 = ZERO;
	D128 q28 = ZERO;
	D128 q29 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x30 = LOAD(&Ix(a,3,0));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y30 = LOAD(&Iy(a,3,0));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z31 = LOAD(&Iz(a,3,1));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    D128 z21 = LOAD(&Iz(a,2,1));
	    q20 = ADD(q20, MUL(MUL(x30, y00), z01));
	    q21 = ADD(q21, MUL(MUL(x00, y30), z01));
	    q22 = ADD(q22, MUL(MUL(x00, y00), z31));
	    q23 = ADD(q23, MUL(MUL(x20, y10), z01));
	    q24 = ADD(q24, MUL(MUL(x20, y00), z11));
	    q25 = ADD(q25, MUL(MUL(x10, y20), z01));
	    q26 = ADD(q26, MUL(MUL(x00, y20), z11));
	    q27 = ADD(q27, MUL(MUL(x10, y00), z21));
	    q28 = ADD(q28, MUL(MUL(x00, y10), z21));
	    q29 = ADD(q29, MUL(MUL(x10, y10), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x30 = LOAD1(&Ix(a,3,0));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y30 = LOAD1(&Iy(a,3,0));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z31 = LOAD1(&Iz(a,3,1));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    D128 z21 = LOAD1(&Iz(a,2,1));
	    q20 = ADD1(q20, MUL1(MUL1(x30, y00), z01));
	    q21 = ADD1(q21, MUL1(MUL1(x00, y30), z01));
	    q22 = ADD1(q22, MUL1(MUL1(x00, y00), z31));
	    q23 = ADD1(q23, MUL1(MUL1(x20, y10), z01));
	    q24 = ADD1(q24, MUL1(MUL1(x20, y00), z11));
	    q25 = ADD1(q25, MUL1(MUL1(x10, y20), z01));
	    q26 = ADD1(q26, MUL1(MUL1(x00, y20), z11));
	    q27 = ADD1(q27, MUL1(MUL1(x10, y00), z21));
	    q28 = ADD1(q28, MUL1(MUL1(x00, y10), z21));
	    q29 = ADD1(q29, MUL1(MUL1(x10, y10), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK20 = ADD(qK20, MUL(C00, HADD(q20, q21)));
	qK22 = ADD(qK22, MUL(C00, HADD(q22, q23)));
	qK24 = ADD(qK24, MUL(C00, HADD(q24, q25)));
	qK26 = ADD(qK26, MUL(C00, HADD(q26, q27)));
	qK28 = ADD(qK28, MUL(C00, HADD(q28, q29)));

#else // SSE
	    
	// function registers
	T q20 = 0.0;
	T q21 = 0.0;
	T q22 = 0.0;
	T q23 = 0.0;
	T q24 = 0.0;
	T q25 = 0.0;
	T q26 = 0.0;
	T q27 = 0.0;
	T q28 = 0.0;
	T q29 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q20 += Ix(a,3,0)*Iy(a,0,0)*Iz(a,0,1);
	    q21 += Ix(a,0,0)*Iy(a,3,0)*Iz(a,0,1);
	    q22 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,3,1);
	    q23 += Ix(a,2,0)*Iy(a,1,0)*Iz(a,0,1);
	    q24 += Ix(a,2,0)*Iy(a,0,0)*Iz(a,1,1);
	    q25 += Ix(a,1,0)*Iy(a,2,0)*Iz(a,0,1);
	    q26 += Ix(a,0,0)*Iy(a,2,0)*Iz(a,1,1);
	    q27 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,2,1);
	    q28 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,2,1);
	    q29 += Ix(a,1,0)*Iy(a,1,0)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK20 += q20*C[k+0];
	//I[20] += q20*C[k+0];
	I[20] += q20*C_[0];
	//qK21 += q21*C[k+0];
	//I[21] += q21*C[k+0];
	I[21] += q21*C_[0];
	//qK22 += q22*C[k+0];
	//I[22] += q22*C[k+0];
	I[22] += q22*C_[0];
	//qK23 += q23*C[k+0];
	//I[23] += q23*C[k+0];
	I[23] += q23*C_[0];
	//qK24 += q24*C[k+0];
	//I[24] += q24*C[k+0];
	I[24] += q24*C_[0];
	//qK25 += q25*C[k+0];
	//I[25] += q25*C[k+0];
	I[25] += q25*C_[0];
	//qK26 += q26*C[k+0];
	//I[26] += q26*C[k+0];
	I[26] += q26*C_[0];
	//qK27 += q27*C[k+0];
	//I[27] += q27*C[k+0];
	I[27] += q27*C_[0];
	//qK28 += q28*C[k+0];
	//I[28] += q28*C[k+0];
	I[28] += q28*C_[0];
	//qK29 += q29*C[k+0];
	//I[29] += q29*C[k+0];
	I[29] += q29*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK20 = MUL(q, qK20);
	qK22 = MUL(q, qK22);
	qK24 = MUL(q, qK24);
	qK26 = MUL(q, qK26);
	qK28 = MUL(q, qK28);
	num += 10; //num += (fabs(I[28]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[20]) & 0xF) {
	    // 20
	    STOREU(&I[20], ADD(qK20, LOADU(&I[20])));
	    STOREU(&I[22], ADD(qK22, LOADU(&I[22])));
	    STOREU(&I[24], ADD(qK24, LOADU(&I[24])));
	    STOREU(&I[26], ADD(qK26, LOADU(&I[26])));
	    STOREU(&I[28], ADD(qK28, LOADU(&I[28])));
	}
	else {
	    STORE(&I[20], ADD(qK20, LOADU(&I[20])));
	    STORE(&I[22], ADD(qK22, LOADU(&I[22])));
	    STORE(&I[24], ADD(qK24, LOADU(&I[24])));
	    STORE(&I[26], ADD(qK26, LOADU(&I[26])));
	    STORE(&I[28], ADD(qK28, LOADU(&I[28])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[20] += scale*NORMALIZE[10]*NORMALIZE[3]*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*NORMALIZE[11]*NORMALIZE[3]*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*NORMALIZE[12]*NORMALIZE[3]*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*NORMALIZE[13]*NORMALIZE[3]*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*NORMALIZE[14]*NORMALIZE[3]*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*NORMALIZE[15]*NORMALIZE[3]*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*NORMALIZE[16]*NORMALIZE[3]*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*NORMALIZE[17]*NORMALIZE[3]*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*NORMALIZE[18]*NORMALIZE[3]*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*NORMALIZE[19]*NORMALIZE[3]*qK29;
	// num += (fabs(I[29]) >= tol);
    }
    else {
	// I[20] += scale*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*qK29;
	// num += (fabs(I[29]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <spp| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::SP,rysq::P> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::SP,rysq::P> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 2;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x01, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x11, y00), z00));
	    q2 = ADD(q2, MUL(MUL(x01, y10), z00));
	    q3 = ADD(q3, MUL(MUL(x01, y00), z10));
	    q4 = ADD(q4, MUL(MUL(x00, y01), z00));
	    q5 = ADD(q5, MUL(MUL(x10, y01), z00));
	    q6 = ADD(q6, MUL(MUL(x00, y11), z00));
	    q7 = ADD(q7, MUL(MUL(x00, y01), z10));
	    q8 = ADD(q8, MUL(MUL(x00, y00), z01));
	    q9 = ADD(q9, MUL(MUL(x10, y00), z01));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x01, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x11, y00), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x01, y10), z00));
	    q3 = ADD1(q3, MUL1(MUL1(x01, y00), z10));
	    q4 = ADD1(q4, MUL1(MUL1(x00, y01), z00));
	    q5 = ADD1(q5, MUL1(MUL1(x10, y01), z00));
	    q6 = ADD1(q6, MUL1(MUL1(x00, y11), z00));
	    q7 = ADD1(q7, MUL1(MUL1(x00, y01), z10));
	    q8 = ADD1(q8, MUL1(MUL1(x00, y00), z01));
	    q9 = ADD1(q9, MUL1(MUL1(x10, y00), z01));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C01 = LOADU(&C[k+0]);
	qK0 = ADD(qK0, MUL(C01, HADD(q0, q1)));
	D128 C11 = LOADDUP(&C[k+1]);
	qK2 = ADD(qK2, MUL(C11, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C01, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C11, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C01, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,0,0);
	    q2 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,0,0);
	    q3 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,1,0);
	    q4 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,0,0);
	    q5 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,0,0);
	    q6 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,0,0);
	    q7 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,1,0);
	    q8 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,0,1);
	    q9 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,0,1);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+1];
	//I[1] += q1*C[k+1];
	I[1] += q1*C_[1];
	//qK2 += q2*C[k+1];
	//I[2] += q2*C[k+1];
	I[2] += q2*C_[1];
	//qK3 += q3*C[k+1];
	//I[3] += q3*C[k+1];
	I[3] += q3*C_[1];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+1];
	//I[5] += q5*C[k+1];
	I[5] += q5*C_[1];
	//qK6 += q6*C[k+1];
	//I[6] += q6*C[k+1];
	I[6] += q6*C_[1];
	//qK7 += q7*C[k+1];
	//I[7] += q7*C[k+1];
	I[7] += q7*C_[1];
	//qK8 += q8*C[k+0];
	//I[8] += q8*C[k+0];
	I[8] += q8*C_[0];
	//qK9 += q9*C[k+1];
	//I[9] += q9*C[k+1];
	I[9] += q9*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[0]*NORMALIZE[1]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[1]*NORMALIZE[1]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[2]*NORMALIZE[1]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[3]*NORMALIZE[1]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[0]*NORMALIZE[2]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[1]*NORMALIZE[2]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[2]*NORMALIZE[2]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[3]*NORMALIZE[2]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[0]*NORMALIZE[3]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[1]*NORMALIZE[3]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK10 = ZERO;
#else
     //double qK10 = 0.0;
     //double qK11 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q10 = ZERO;
	D128 q11 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    q10 = ADD(q10, MUL(MUL(x00, y10), z01));
	    q11 = ADD(q11, MUL(MUL(x00, y00), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    q10 = ADD1(q10, MUL1(MUL1(x00, y10), z01));
	    q11 = ADD1(q11, MUL1(MUL1(x00, y00), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C11 = LOADDUP(&C[k+1]);
	qK10 = ADD(qK10, MUL(C11, HADD(q10, q11)));

#else // SSE
	    
	// function registers
	T q10 = 0.0;
	T q11 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,0,1);
	    q11 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK10 += q10*C[k+1];
	//I[10] += q10*C[k+1];
	I[10] += q10*C_[1];
	//qK11 += q11*C[k+1];
	//I[11] += q11*C[k+1];
	I[11] += q11*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK10 = MUL(q, qK10);
	num += 2; //num += (fabs(I[10]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[10]) & 0xF) {
	    // 10
	    STOREU(&I[10], ADD(qK10, LOADU(&I[10])));
	}
	else {
	    STORE(&I[10], ADD(qK10, LOADU(&I[10])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[10] += scale*NORMALIZE[2]*NORMALIZE[3]*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*NORMALIZE[3]*NORMALIZE[3]*qK11;
	// num += (fabs(I[11]) >= tol);
    }
    else {
	// I[10] += scale*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*qK11;
	// num += (fabs(I[11]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <sd| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::S,rysq::D> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::S,rysq::D> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 1;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x02 = LOAD(&Ix(a,0,2));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y02 = LOAD(&Iy(a,0,2));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z02 = LOAD(&Iz(a,0,2));
	    q0 = ADD(q0, MUL(MUL(x02, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x00, y02), z00));
	    q2 = ADD(q2, MUL(MUL(x00, y00), z02));
	    q3 = ADD(q3, MUL(MUL(x01, y01), z00));
	    q4 = ADD(q4, MUL(MUL(x01, y00), z01));
	    q5 = ADD(q5, MUL(MUL(x00, y01), z01));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x02 = LOAD1(&Ix(a,0,2));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y02 = LOAD1(&Iy(a,0,2));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z02 = LOAD1(&Iz(a,0,2));
	    q0 = ADD1(q0, MUL1(MUL1(x02, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x00, y02), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x00, y00), z02));
	    q3 = ADD1(q3, MUL1(MUL1(x01, y01), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x01, y00), z01));
	    q5 = ADD1(q5, MUL1(MUL1(x00, y01), z01));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD(qK2, MUL(C00, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C00, HADD(q4, q5)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,0,2)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,0)*Iy(a,0,2)*Iz(a,0,0);
	    q2 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,0,2);
	    q3 += Ix(a,0,1)*Iy(a,0,1)*Iz(a,0,0);
	    q4 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,0,1);
	    q5 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,0,1);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+0];
	//I[3] += q3*C[k+0];
	I[3] += q3*C_[0];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+0];
	//I[5] += q5*C[k+0];
	I[5] += q5*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	num += 6; //num += (fabs(I[4]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[0]*NORMALIZE[4]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[0]*NORMALIZE[5]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[0]*NORMALIZE[6]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[0]*NORMALIZE[7]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[0]*NORMALIZE[8]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[0]*NORMALIZE[9]*qK5;
	// num += (fabs(I[5]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <pd| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::P,rysq::D> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::P,rysq::D> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 2;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x12 = LOAD(&Ix(a,1,2));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x02 = LOAD(&Ix(a,0,2));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 y12 = LOAD(&Iy(a,1,2));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y02 = LOAD(&Iy(a,0,2));
	    D128 z12 = LOAD(&Iz(a,1,2));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z02 = LOAD(&Iz(a,0,2));
	    q0 = ADD(q0, MUL(MUL(x12, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x02, y10), z00));
	    q2 = ADD(q2, MUL(MUL(x02, y00), z10));
	    q3 = ADD(q3, MUL(MUL(x10, y02), z00));
	    q4 = ADD(q4, MUL(MUL(x00, y12), z00));
	    q5 = ADD(q5, MUL(MUL(x00, y02), z10));
	    q6 = ADD(q6, MUL(MUL(x10, y00), z02));
	    q7 = ADD(q7, MUL(MUL(x00, y10), z02));
	    q8 = ADD(q8, MUL(MUL(x00, y00), z12));
	    q9 = ADD(q9, MUL(MUL(x11, y01), z00));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x12 = LOAD1(&Ix(a,1,2));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x02 = LOAD1(&Ix(a,0,2));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 y12 = LOAD1(&Iy(a,1,2));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y02 = LOAD1(&Iy(a,0,2));
	    D128 z12 = LOAD1(&Iz(a,1,2));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z02 = LOAD1(&Iz(a,0,2));
	    q0 = ADD1(q0, MUL1(MUL1(x12, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x02, y10), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x02, y00), z10));
	    q3 = ADD1(q3, MUL1(MUL1(x10, y02), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x00, y12), z00));
	    q5 = ADD1(q5, MUL1(MUL1(x00, y02), z10));
	    q6 = ADD1(q6, MUL1(MUL1(x10, y00), z02));
	    q7 = ADD1(q7, MUL1(MUL1(x00, y10), z02));
	    q8 = ADD1(q8, MUL1(MUL1(x00, y00), z12));
	    q9 = ADD1(q9, MUL1(MUL1(x11, y01), z00));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD(qK2, MUL(C00, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C00, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C00, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C00, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,1,2)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,2)*Iy(a,1,0)*Iz(a,0,0);
	    q2 += Ix(a,0,2)*Iy(a,0,0)*Iz(a,1,0);
	    q3 += Ix(a,1,0)*Iy(a,0,2)*Iz(a,0,0);
	    q4 += Ix(a,0,0)*Iy(a,1,2)*Iz(a,0,0);
	    q5 += Ix(a,0,0)*Iy(a,0,2)*Iz(a,1,0);
	    q6 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,0,2);
	    q7 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,0,2);
	    q8 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,1,2);
	    q9 += Ix(a,1,1)*Iy(a,0,1)*Iz(a,0,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+0];
	//I[3] += q3*C[k+0];
	I[3] += q3*C_[0];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+0];
	//I[5] += q5*C[k+0];
	I[5] += q5*C_[0];
	//qK6 += q6*C[k+0];
	//I[6] += q6*C[k+0];
	I[6] += q6*C_[0];
	//qK7 += q7*C[k+0];
	//I[7] += q7*C[k+0];
	I[7] += q7*C_[0];
	//qK8 += q8*C[k+0];
	//I[8] += q8*C[k+0];
	I[8] += q8*C_[0];
	//qK9 += q9*C[k+0];
	//I[9] += q9*C[k+0];
	I[9] += q9*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[1]*NORMALIZE[4]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[2]*NORMALIZE[4]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[3]*NORMALIZE[4]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[1]*NORMALIZE[5]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[2]*NORMALIZE[5]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[3]*NORMALIZE[5]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[1]*NORMALIZE[6]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[2]*NORMALIZE[6]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[3]*NORMALIZE[6]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[1]*NORMALIZE[7]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK10 = ZERO;
     D128 qK12 = ZERO;
     D128 qK14 = ZERO;
     D128 qK16 = ZERO;
#else
     //double qK10 = 0.0;
     //double qK11 = 0.0;
     //double qK12 = 0.0;
     //double qK13 = 0.0;
     //double qK14 = 0.0;
     //double qK15 = 0.0;
     //double qK16 = 0.0;
     //double qK17 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q10 = ZERO;
	D128 q11 = ZERO;
	D128 q12 = ZERO;
	D128 q13 = ZERO;
	D128 q14 = ZERO;
	D128 q15 = ZERO;
	D128 q16 = ZERO;
	D128 q17 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    q10 = ADD(q10, MUL(MUL(x01, y11), z00));
	    q11 = ADD(q11, MUL(MUL(x01, y01), z10));
	    q12 = ADD(q12, MUL(MUL(x11, y00), z01));
	    q13 = ADD(q13, MUL(MUL(x01, y10), z01));
	    q14 = ADD(q14, MUL(MUL(x01, y00), z11));
	    q15 = ADD(q15, MUL(MUL(x10, y01), z01));
	    q16 = ADD(q16, MUL(MUL(x00, y11), z01));
	    q17 = ADD(q17, MUL(MUL(x00, y01), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    q10 = ADD1(q10, MUL1(MUL1(x01, y11), z00));
	    q11 = ADD1(q11, MUL1(MUL1(x01, y01), z10));
	    q12 = ADD1(q12, MUL1(MUL1(x11, y00), z01));
	    q13 = ADD1(q13, MUL1(MUL1(x01, y10), z01));
	    q14 = ADD1(q14, MUL1(MUL1(x01, y00), z11));
	    q15 = ADD1(q15, MUL1(MUL1(x10, y01), z01));
	    q16 = ADD1(q16, MUL1(MUL1(x00, y11), z01));
	    q17 = ADD1(q17, MUL1(MUL1(x00, y01), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK10 = ADD(qK10, MUL(C00, HADD(q10, q11)));
	qK12 = ADD(qK12, MUL(C00, HADD(q12, q13)));
	qK14 = ADD(qK14, MUL(C00, HADD(q14, q15)));
	qK16 = ADD(qK16, MUL(C00, HADD(q16, q17)));

#else // SSE
	    
	// function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
	T q16 = 0.0;
	T q17 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += Ix(a,0,1)*Iy(a,1,1)*Iz(a,0,0);
	    q11 += Ix(a,0,1)*Iy(a,0,1)*Iz(a,1,0);
	    q12 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,0,1);
	    q13 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,0,1);
	    q14 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,1,1);
	    q15 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,0,1);
	    q16 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,0,1);
	    q17 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK10 += q10*C[k+0];
	//I[10] += q10*C[k+0];
	I[10] += q10*C_[0];
	//qK11 += q11*C[k+0];
	//I[11] += q11*C[k+0];
	I[11] += q11*C_[0];
	//qK12 += q12*C[k+0];
	//I[12] += q12*C[k+0];
	I[12] += q12*C_[0];
	//qK13 += q13*C[k+0];
	//I[13] += q13*C[k+0];
	I[13] += q13*C_[0];
	//qK14 += q14*C[k+0];
	//I[14] += q14*C[k+0];
	I[14] += q14*C_[0];
	//qK15 += q15*C[k+0];
	//I[15] += q15*C[k+0];
	I[15] += q15*C_[0];
	//qK16 += q16*C[k+0];
	//I[16] += q16*C[k+0];
	I[16] += q16*C_[0];
	//qK17 += q17*C[k+0];
	//I[17] += q17*C[k+0];
	I[17] += q17*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK10 = MUL(q, qK10);
	qK12 = MUL(q, qK12);
	qK14 = MUL(q, qK14);
	qK16 = MUL(q, qK16);
	num += 8; //num += (fabs(I[16]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[10]) & 0xF) {
	    // 10
	    STOREU(&I[10], ADD(qK10, LOADU(&I[10])));
	    STOREU(&I[12], ADD(qK12, LOADU(&I[12])));
	    STOREU(&I[14], ADD(qK14, LOADU(&I[14])));
	    STOREU(&I[16], ADD(qK16, LOADU(&I[16])));
	}
	else {
	    STORE(&I[10], ADD(qK10, LOADU(&I[10])));
	    STORE(&I[12], ADD(qK12, LOADU(&I[12])));
	    STORE(&I[14], ADD(qK14, LOADU(&I[14])));
	    STORE(&I[16], ADD(qK16, LOADU(&I[16])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[10] += scale*NORMALIZE[2]*NORMALIZE[7]*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*NORMALIZE[3]*NORMALIZE[7]*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*NORMALIZE[1]*NORMALIZE[8]*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*NORMALIZE[2]*NORMALIZE[8]*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*NORMALIZE[3]*NORMALIZE[8]*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*NORMALIZE[1]*NORMALIZE[9]*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*NORMALIZE[2]*NORMALIZE[9]*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*NORMALIZE[3]*NORMALIZE[9]*qK17;
	// num += (fabs(I[17]) >= tol);
    }
    else {
	// I[10] += scale*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*qK17;
	// num += (fabs(I[17]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <dd| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::D,rysq::D> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::D,rysq::D> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 3;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x12 = LOAD(&Ix(a,1,2));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x22 = LOAD(&Ix(a,2,2));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x02 = LOAD(&Ix(a,0,2));
	    D128 y12 = LOAD(&Iy(a,1,2));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y22 = LOAD(&Iy(a,2,2));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y02 = LOAD(&Iy(a,0,2));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x22, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x02, y20), z00));
	    q2 = ADD(q2, MUL(MUL(x02, y00), z20));
	    q3 = ADD(q3, MUL(MUL(x12, y10), z00));
	    q4 = ADD(q4, MUL(MUL(x12, y00), z10));
	    q5 = ADD(q5, MUL(MUL(x02, y10), z10));
	    q6 = ADD(q6, MUL(MUL(x20, y02), z00));
	    q7 = ADD(q7, MUL(MUL(x00, y22), z00));
	    q8 = ADD(q8, MUL(MUL(x00, y02), z20));
	    q9 = ADD(q9, MUL(MUL(x10, y12), z00));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x12 = LOAD1(&Ix(a,1,2));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x22 = LOAD1(&Ix(a,2,2));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x02 = LOAD1(&Ix(a,0,2));
	    D128 y12 = LOAD1(&Iy(a,1,2));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y22 = LOAD1(&Iy(a,2,2));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y02 = LOAD1(&Iy(a,0,2));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x22, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x02, y20), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x02, y00), z20));
	    q3 = ADD1(q3, MUL1(MUL1(x12, y10), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x12, y00), z10));
	    q5 = ADD1(q5, MUL1(MUL1(x02, y10), z10));
	    q6 = ADD1(q6, MUL1(MUL1(x20, y02), z00));
	    q7 = ADD1(q7, MUL1(MUL1(x00, y22), z00));
	    q8 = ADD1(q8, MUL1(MUL1(x00, y02), z20));
	    q9 = ADD1(q9, MUL1(MUL1(x10, y12), z00));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD(qK2, MUL(C00, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C00, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C00, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C00, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,2,2)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,2)*Iy(a,2,0)*Iz(a,0,0);
	    q2 += Ix(a,0,2)*Iy(a,0,0)*Iz(a,2,0);
	    q3 += Ix(a,1,2)*Iy(a,1,0)*Iz(a,0,0);
	    q4 += Ix(a,1,2)*Iy(a,0,0)*Iz(a,1,0);
	    q5 += Ix(a,0,2)*Iy(a,1,0)*Iz(a,1,0);
	    q6 += Ix(a,2,0)*Iy(a,0,2)*Iz(a,0,0);
	    q7 += Ix(a,0,0)*Iy(a,2,2)*Iz(a,0,0);
	    q8 += Ix(a,0,0)*Iy(a,0,2)*Iz(a,2,0);
	    q9 += Ix(a,1,0)*Iy(a,1,2)*Iz(a,0,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+0];
	//I[3] += q3*C[k+0];
	I[3] += q3*C_[0];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+0];
	//I[5] += q5*C[k+0];
	I[5] += q5*C_[0];
	//qK6 += q6*C[k+0];
	//I[6] += q6*C[k+0];
	I[6] += q6*C_[0];
	//qK7 += q7*C[k+0];
	//I[7] += q7*C[k+0];
	I[7] += q7*C_[0];
	//qK8 += q8*C[k+0];
	//I[8] += q8*C[k+0];
	I[8] += q8*C_[0];
	//qK9 += q9*C[k+0];
	//I[9] += q9*C[k+0];
	I[9] += q9*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[4]*NORMALIZE[4]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[5]*NORMALIZE[4]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[6]*NORMALIZE[4]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[7]*NORMALIZE[4]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[8]*NORMALIZE[4]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[9]*NORMALIZE[4]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[4]*NORMALIZE[5]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[5]*NORMALIZE[5]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[6]*NORMALIZE[5]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[7]*NORMALIZE[5]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK10 = ZERO;
     D128 qK12 = ZERO;
     D128 qK14 = ZERO;
     D128 qK16 = ZERO;
     D128 qK18 = ZERO;
#else
     //double qK10 = 0.0;
     //double qK11 = 0.0;
     //double qK12 = 0.0;
     //double qK13 = 0.0;
     //double qK14 = 0.0;
     //double qK15 = 0.0;
     //double qK16 = 0.0;
     //double qK17 = 0.0;
     //double qK18 = 0.0;
     //double qK19 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q10 = ZERO;
	D128 q11 = ZERO;
	D128 q12 = ZERO;
	D128 q13 = ZERO;
	D128 q14 = ZERO;
	D128 q15 = ZERO;
	D128 q16 = ZERO;
	D128 q17 = ZERO;
	D128 q18 = ZERO;
	D128 q19 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x21 = LOAD(&Ix(a,2,1));
	    D128 y12 = LOAD(&Iy(a,1,2));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y21 = LOAD(&Iy(a,2,1));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y02 = LOAD(&Iy(a,0,2));
	    D128 z12 = LOAD(&Iz(a,1,2));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z02 = LOAD(&Iz(a,0,2));
	    D128 z22 = LOAD(&Iz(a,2,2));
	    q10 = ADD(q10, MUL(MUL(x10, y02), z10));
	    q11 = ADD(q11, MUL(MUL(x00, y12), z10));
	    q12 = ADD(q12, MUL(MUL(x20, y00), z02));
	    q13 = ADD(q13, MUL(MUL(x00, y20), z02));
	    q14 = ADD(q14, MUL(MUL(x00, y00), z22));
	    q15 = ADD(q15, MUL(MUL(x10, y10), z02));
	    q16 = ADD(q16, MUL(MUL(x10, y00), z12));
	    q17 = ADD(q17, MUL(MUL(x00, y10), z12));
	    q18 = ADD(q18, MUL(MUL(x21, y01), z00));
	    q19 = ADD(q19, MUL(MUL(x01, y21), z00));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x21 = LOAD1(&Ix(a,2,1));
	    D128 y12 = LOAD1(&Iy(a,1,2));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y21 = LOAD1(&Iy(a,2,1));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y02 = LOAD1(&Iy(a,0,2));
	    D128 z12 = LOAD1(&Iz(a,1,2));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z02 = LOAD1(&Iz(a,0,2));
	    D128 z22 = LOAD1(&Iz(a,2,2));
	    q10 = ADD1(q10, MUL1(MUL1(x10, y02), z10));
	    q11 = ADD1(q11, MUL1(MUL1(x00, y12), z10));
	    q12 = ADD1(q12, MUL1(MUL1(x20, y00), z02));
	    q13 = ADD1(q13, MUL1(MUL1(x00, y20), z02));
	    q14 = ADD1(q14, MUL1(MUL1(x00, y00), z22));
	    q15 = ADD1(q15, MUL1(MUL1(x10, y10), z02));
	    q16 = ADD1(q16, MUL1(MUL1(x10, y00), z12));
	    q17 = ADD1(q17, MUL1(MUL1(x00, y10), z12));
	    q18 = ADD1(q18, MUL1(MUL1(x21, y01), z00));
	    q19 = ADD1(q19, MUL1(MUL1(x01, y21), z00));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK10 = ADD(qK10, MUL(C00, HADD(q10, q11)));
	qK12 = ADD(qK12, MUL(C00, HADD(q12, q13)));
	qK14 = ADD(qK14, MUL(C00, HADD(q14, q15)));
	qK16 = ADD(qK16, MUL(C00, HADD(q16, q17)));
	qK18 = ADD(qK18, MUL(C00, HADD(q18, q19)));

#else // SSE
	    
	// function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
	T q16 = 0.0;
	T q17 = 0.0;
	T q18 = 0.0;
	T q19 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += Ix(a,1,0)*Iy(a,0,2)*Iz(a,1,0);
	    q11 += Ix(a,0,0)*Iy(a,1,2)*Iz(a,1,0);
	    q12 += Ix(a,2,0)*Iy(a,0,0)*Iz(a,0,2);
	    q13 += Ix(a,0,0)*Iy(a,2,0)*Iz(a,0,2);
	    q14 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,2,2);
	    q15 += Ix(a,1,0)*Iy(a,1,0)*Iz(a,0,2);
	    q16 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,1,2);
	    q17 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,1,2);
	    q18 += Ix(a,2,1)*Iy(a,0,1)*Iz(a,0,0);
	    q19 += Ix(a,0,1)*Iy(a,2,1)*Iz(a,0,0);
	}
	    
	//contraction coefficients
	//qK10 += q10*C[k+0];
	//I[10] += q10*C[k+0];
	I[10] += q10*C_[0];
	//qK11 += q11*C[k+0];
	//I[11] += q11*C[k+0];
	I[11] += q11*C_[0];
	//qK12 += q12*C[k+0];
	//I[12] += q12*C[k+0];
	I[12] += q12*C_[0];
	//qK13 += q13*C[k+0];
	//I[13] += q13*C[k+0];
	I[13] += q13*C_[0];
	//qK14 += q14*C[k+0];
	//I[14] += q14*C[k+0];
	I[14] += q14*C_[0];
	//qK15 += q15*C[k+0];
	//I[15] += q15*C[k+0];
	I[15] += q15*C_[0];
	//qK16 += q16*C[k+0];
	//I[16] += q16*C[k+0];
	I[16] += q16*C_[0];
	//qK17 += q17*C[k+0];
	//I[17] += q17*C[k+0];
	I[17] += q17*C_[0];
	//qK18 += q18*C[k+0];
	//I[18] += q18*C[k+0];
	I[18] += q18*C_[0];
	//qK19 += q19*C[k+0];
	//I[19] += q19*C[k+0];
	I[19] += q19*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK10 = MUL(q, qK10);
	qK12 = MUL(q, qK12);
	qK14 = MUL(q, qK14);
	qK16 = MUL(q, qK16);
	qK18 = MUL(q, qK18);
	num += 10; //num += (fabs(I[18]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[10]) & 0xF) {
	    // 10
	    STOREU(&I[10], ADD(qK10, LOADU(&I[10])));
	    STOREU(&I[12], ADD(qK12, LOADU(&I[12])));
	    STOREU(&I[14], ADD(qK14, LOADU(&I[14])));
	    STOREU(&I[16], ADD(qK16, LOADU(&I[16])));
	    STOREU(&I[18], ADD(qK18, LOADU(&I[18])));
	}
	else {
	    STORE(&I[10], ADD(qK10, LOADU(&I[10])));
	    STORE(&I[12], ADD(qK12, LOADU(&I[12])));
	    STORE(&I[14], ADD(qK14, LOADU(&I[14])));
	    STORE(&I[16], ADD(qK16, LOADU(&I[16])));
	    STORE(&I[18], ADD(qK18, LOADU(&I[18])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[10] += scale*NORMALIZE[8]*NORMALIZE[5]*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*NORMALIZE[9]*NORMALIZE[5]*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*NORMALIZE[4]*NORMALIZE[6]*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*NORMALIZE[5]*NORMALIZE[6]*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*NORMALIZE[6]*NORMALIZE[6]*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*NORMALIZE[7]*NORMALIZE[6]*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*NORMALIZE[8]*NORMALIZE[6]*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*NORMALIZE[9]*NORMALIZE[6]*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*NORMALIZE[4]*NORMALIZE[7]*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*NORMALIZE[5]*NORMALIZE[7]*qK19;
	// num += (fabs(I[19]) >= tol);
    }
    else {
	// I[10] += scale*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*qK19;
	// num += (fabs(I[19]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK20 = ZERO;
     D128 qK22 = ZERO;
     D128 qK24 = ZERO;
     D128 qK26 = ZERO;
     D128 qK28 = ZERO;
#else
     //double qK20 = 0.0;
     //double qK21 = 0.0;
     //double qK22 = 0.0;
     //double qK23 = 0.0;
     //double qK24 = 0.0;
     //double qK25 = 0.0;
     //double qK26 = 0.0;
     //double qK27 = 0.0;
     //double qK28 = 0.0;
     //double qK29 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q20 = ZERO;
	D128 q21 = ZERO;
	D128 q22 = ZERO;
	D128 q23 = ZERO;
	D128 q24 = ZERO;
	D128 q25 = ZERO;
	D128 q26 = ZERO;
	D128 q27 = ZERO;
	D128 q28 = ZERO;
	D128 q29 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 x21 = LOAD(&Ix(a,2,1));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z21 = LOAD(&Iz(a,2,1));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    q20 = ADD(q20, MUL(MUL(x01, y01), z20));
	    q21 = ADD(q21, MUL(MUL(x11, y11), z00));
	    q22 = ADD(q22, MUL(MUL(x11, y01), z10));
	    q23 = ADD(q23, MUL(MUL(x01, y11), z10));
	    q24 = ADD(q24, MUL(MUL(x21, y00), z01));
	    q25 = ADD(q25, MUL(MUL(x01, y20), z01));
	    q26 = ADD(q26, MUL(MUL(x01, y00), z21));
	    q27 = ADD(q27, MUL(MUL(x11, y10), z01));
	    q28 = ADD(q28, MUL(MUL(x11, y00), z11));
	    q29 = ADD(q29, MUL(MUL(x01, y10), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 x21 = LOAD1(&Ix(a,2,1));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z21 = LOAD1(&Iz(a,2,1));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    q20 = ADD1(q20, MUL1(MUL1(x01, y01), z20));
	    q21 = ADD1(q21, MUL1(MUL1(x11, y11), z00));
	    q22 = ADD1(q22, MUL1(MUL1(x11, y01), z10));
	    q23 = ADD1(q23, MUL1(MUL1(x01, y11), z10));
	    q24 = ADD1(q24, MUL1(MUL1(x21, y00), z01));
	    q25 = ADD1(q25, MUL1(MUL1(x01, y20), z01));
	    q26 = ADD1(q26, MUL1(MUL1(x01, y00), z21));
	    q27 = ADD1(q27, MUL1(MUL1(x11, y10), z01));
	    q28 = ADD1(q28, MUL1(MUL1(x11, y00), z11));
	    q29 = ADD1(q29, MUL1(MUL1(x01, y10), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK20 = ADD(qK20, MUL(C00, HADD(q20, q21)));
	qK22 = ADD(qK22, MUL(C00, HADD(q22, q23)));
	qK24 = ADD(qK24, MUL(C00, HADD(q24, q25)));
	qK26 = ADD(qK26, MUL(C00, HADD(q26, q27)));
	qK28 = ADD(qK28, MUL(C00, HADD(q28, q29)));

#else // SSE
	    
	// function registers
	T q20 = 0.0;
	T q21 = 0.0;
	T q22 = 0.0;
	T q23 = 0.0;
	T q24 = 0.0;
	T q25 = 0.0;
	T q26 = 0.0;
	T q27 = 0.0;
	T q28 = 0.0;
	T q29 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q20 += Ix(a,0,1)*Iy(a,0,1)*Iz(a,2,0);
	    q21 += Ix(a,1,1)*Iy(a,1,1)*Iz(a,0,0);
	    q22 += Ix(a,1,1)*Iy(a,0,1)*Iz(a,1,0);
	    q23 += Ix(a,0,1)*Iy(a,1,1)*Iz(a,1,0);
	    q24 += Ix(a,2,1)*Iy(a,0,0)*Iz(a,0,1);
	    q25 += Ix(a,0,1)*Iy(a,2,0)*Iz(a,0,1);
	    q26 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,2,1);
	    q27 += Ix(a,1,1)*Iy(a,1,0)*Iz(a,0,1);
	    q28 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,1,1);
	    q29 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK20 += q20*C[k+0];
	//I[20] += q20*C[k+0];
	I[20] += q20*C_[0];
	//qK21 += q21*C[k+0];
	//I[21] += q21*C[k+0];
	I[21] += q21*C_[0];
	//qK22 += q22*C[k+0];
	//I[22] += q22*C[k+0];
	I[22] += q22*C_[0];
	//qK23 += q23*C[k+0];
	//I[23] += q23*C[k+0];
	I[23] += q23*C_[0];
	//qK24 += q24*C[k+0];
	//I[24] += q24*C[k+0];
	I[24] += q24*C_[0];
	//qK25 += q25*C[k+0];
	//I[25] += q25*C[k+0];
	I[25] += q25*C_[0];
	//qK26 += q26*C[k+0];
	//I[26] += q26*C[k+0];
	I[26] += q26*C_[0];
	//qK27 += q27*C[k+0];
	//I[27] += q27*C[k+0];
	I[27] += q27*C_[0];
	//qK28 += q28*C[k+0];
	//I[28] += q28*C[k+0];
	I[28] += q28*C_[0];
	//qK29 += q29*C[k+0];
	//I[29] += q29*C[k+0];
	I[29] += q29*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK20 = MUL(q, qK20);
	qK22 = MUL(q, qK22);
	qK24 = MUL(q, qK24);
	qK26 = MUL(q, qK26);
	qK28 = MUL(q, qK28);
	num += 10; //num += (fabs(I[28]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[20]) & 0xF) {
	    // 20
	    STOREU(&I[20], ADD(qK20, LOADU(&I[20])));
	    STOREU(&I[22], ADD(qK22, LOADU(&I[22])));
	    STOREU(&I[24], ADD(qK24, LOADU(&I[24])));
	    STOREU(&I[26], ADD(qK26, LOADU(&I[26])));
	    STOREU(&I[28], ADD(qK28, LOADU(&I[28])));
	}
	else {
	    STORE(&I[20], ADD(qK20, LOADU(&I[20])));
	    STORE(&I[22], ADD(qK22, LOADU(&I[22])));
	    STORE(&I[24], ADD(qK24, LOADU(&I[24])));
	    STORE(&I[26], ADD(qK26, LOADU(&I[26])));
	    STORE(&I[28], ADD(qK28, LOADU(&I[28])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[20] += scale*NORMALIZE[6]*NORMALIZE[7]*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*NORMALIZE[7]*NORMALIZE[7]*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*NORMALIZE[8]*NORMALIZE[7]*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*NORMALIZE[9]*NORMALIZE[7]*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*NORMALIZE[4]*NORMALIZE[8]*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*NORMALIZE[5]*NORMALIZE[8]*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*NORMALIZE[6]*NORMALIZE[8]*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*NORMALIZE[7]*NORMALIZE[8]*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*NORMALIZE[8]*NORMALIZE[8]*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*NORMALIZE[9]*NORMALIZE[8]*qK29;
	// num += (fabs(I[29]) >= tol);
    }
    else {
	// I[20] += scale*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*qK29;
	// num += (fabs(I[29]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK30 = ZERO;
     D128 qK32 = ZERO;
     D128 qK34 = ZERO;
#else
     //double qK30 = 0.0;
     //double qK31 = 0.0;
     //double qK32 = 0.0;
     //double qK33 = 0.0;
     //double qK34 = 0.0;
     //double qK35 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q30 = ZERO;
	D128 q31 = ZERO;
	D128 q32 = ZERO;
	D128 q33 = ZERO;
	D128 q34 = ZERO;
	D128 q35 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 y21 = LOAD(&Iy(a,2,1));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    D128 z21 = LOAD(&Iz(a,2,1));
	    q30 = ADD(q30, MUL(MUL(x20, y01), z01));
	    q31 = ADD(q31, MUL(MUL(x00, y21), z01));
	    q32 = ADD(q32, MUL(MUL(x00, y01), z21));
	    q33 = ADD(q33, MUL(MUL(x10, y11), z01));
	    q34 = ADD(q34, MUL(MUL(x10, y01), z11));
	    q35 = ADD(q35, MUL(MUL(x00, y11), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 y21 = LOAD1(&Iy(a,2,1));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    D128 z21 = LOAD1(&Iz(a,2,1));
	    q30 = ADD1(q30, MUL1(MUL1(x20, y01), z01));
	    q31 = ADD1(q31, MUL1(MUL1(x00, y21), z01));
	    q32 = ADD1(q32, MUL1(MUL1(x00, y01), z21));
	    q33 = ADD1(q33, MUL1(MUL1(x10, y11), z01));
	    q34 = ADD1(q34, MUL1(MUL1(x10, y01), z11));
	    q35 = ADD1(q35, MUL1(MUL1(x00, y11), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK30 = ADD(qK30, MUL(C00, HADD(q30, q31)));
	qK32 = ADD(qK32, MUL(C00, HADD(q32, q33)));
	qK34 = ADD(qK34, MUL(C00, HADD(q34, q35)));

#else // SSE
	    
	// function registers
	T q30 = 0.0;
	T q31 = 0.0;
	T q32 = 0.0;
	T q33 = 0.0;
	T q34 = 0.0;
	T q35 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q30 += Ix(a,2,0)*Iy(a,0,1)*Iz(a,0,1);
	    q31 += Ix(a,0,0)*Iy(a,2,1)*Iz(a,0,1);
	    q32 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,2,1);
	    q33 += Ix(a,1,0)*Iy(a,1,1)*Iz(a,0,1);
	    q34 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,1,1);
	    q35 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK30 += q30*C[k+0];
	//I[30] += q30*C[k+0];
	I[30] += q30*C_[0];
	//qK31 += q31*C[k+0];
	//I[31] += q31*C[k+0];
	I[31] += q31*C_[0];
	//qK32 += q32*C[k+0];
	//I[32] += q32*C[k+0];
	I[32] += q32*C_[0];
	//qK33 += q33*C[k+0];
	//I[33] += q33*C[k+0];
	I[33] += q33*C_[0];
	//qK34 += q34*C[k+0];
	//I[34] += q34*C[k+0];
	I[34] += q34*C_[0];
	//qK35 += q35*C[k+0];
	//I[35] += q35*C[k+0];
	I[35] += q35*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK30 = MUL(q, qK30);
	qK32 = MUL(q, qK32);
	qK34 = MUL(q, qK34);
	num += 6; //num += (fabs(I[34]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[30]) & 0xF) {
	    // 30
	    STOREU(&I[30], ADD(qK30, LOADU(&I[30])));
	    STOREU(&I[32], ADD(qK32, LOADU(&I[32])));
	    STOREU(&I[34], ADD(qK34, LOADU(&I[34])));
	}
	else {
	    STORE(&I[30], ADD(qK30, LOADU(&I[30])));
	    STORE(&I[32], ADD(qK32, LOADU(&I[32])));
	    STORE(&I[34], ADD(qK34, LOADU(&I[34])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[30] += scale*NORMALIZE[4]*NORMALIZE[9]*qK30;
	// num += (fabs(I[30]) >= tol);
	// I[31] += scale*NORMALIZE[5]*NORMALIZE[9]*qK31;
	// num += (fabs(I[31]) >= tol);
	// I[32] += scale*NORMALIZE[6]*NORMALIZE[9]*qK32;
	// num += (fabs(I[32]) >= tol);
	// I[33] += scale*NORMALIZE[7]*NORMALIZE[9]*qK33;
	// num += (fabs(I[33]) >= tol);
	// I[34] += scale*NORMALIZE[8]*NORMALIZE[9]*qK34;
	// num += (fabs(I[34]) >= tol);
	// I[35] += scale*NORMALIZE[9]*NORMALIZE[9]*qK35;
	// num += (fabs(I[35]) >= tol);
    }
    else {
	// I[30] += scale*qK30;
	// num += (fabs(I[30]) >= tol);
	// I[31] += scale*qK31;
	// num += (fabs(I[31]) >= tol);
	// I[32] += scale*qK32;
	// num += (fabs(I[32]) >= tol);
	// I[33] += scale*qK33;
	// num += (fabs(I[33]) >= tol);
	// I[34] += scale*qK34;
	// num += (fabs(I[34]) >= tol);
	// I[35] += scale*qK35;
	// num += (fabs(I[35]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <fd| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::F,rysq::D> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::F,rysq::D> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 4;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x12 = LOAD(&Ix(a,1,2));
	    D128 x32 = LOAD(&Ix(a,3,2));
	    D128 x02 = LOAD(&Ix(a,0,2));
	    D128 x22 = LOAD(&Ix(a,2,2));
	    D128 y30 = LOAD(&Iy(a,3,0));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z30 = LOAD(&Iz(a,3,0));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x32, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x02, y30), z00));
	    q2 = ADD(q2, MUL(MUL(x02, y00), z30));
	    q3 = ADD(q3, MUL(MUL(x22, y10), z00));
	    q4 = ADD(q4, MUL(MUL(x22, y00), z10));
	    q5 = ADD(q5, MUL(MUL(x12, y20), z00));
	    q6 = ADD(q6, MUL(MUL(x02, y20), z10));
	    q7 = ADD(q7, MUL(MUL(x12, y00), z20));
	    q8 = ADD(q8, MUL(MUL(x02, y10), z20));
	    q9 = ADD(q9, MUL(MUL(x12, y10), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x12 = LOAD1(&Ix(a,1,2));
	    D128 x32 = LOAD1(&Ix(a,3,2));
	    D128 x02 = LOAD1(&Ix(a,0,2));
	    D128 x22 = LOAD1(&Ix(a,2,2));
	    D128 y30 = LOAD1(&Iy(a,3,0));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z30 = LOAD1(&Iz(a,3,0));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x32, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x02, y30), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x02, y00), z30));
	    q3 = ADD1(q3, MUL1(MUL1(x22, y10), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x22, y00), z10));
	    q5 = ADD1(q5, MUL1(MUL1(x12, y20), z00));
	    q6 = ADD1(q6, MUL1(MUL1(x02, y20), z10));
	    q7 = ADD1(q7, MUL1(MUL1(x12, y00), z20));
	    q8 = ADD1(q8, MUL1(MUL1(x02, y10), z20));
	    q9 = ADD1(q9, MUL1(MUL1(x12, y10), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD(qK2, MUL(C00, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C00, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C00, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C00, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,3,2)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,2)*Iy(a,3,0)*Iz(a,0,0);
	    q2 += Ix(a,0,2)*Iy(a,0,0)*Iz(a,3,0);
	    q3 += Ix(a,2,2)*Iy(a,1,0)*Iz(a,0,0);
	    q4 += Ix(a,2,2)*Iy(a,0,0)*Iz(a,1,0);
	    q5 += Ix(a,1,2)*Iy(a,2,0)*Iz(a,0,0);
	    q6 += Ix(a,0,2)*Iy(a,2,0)*Iz(a,1,0);
	    q7 += Ix(a,1,2)*Iy(a,0,0)*Iz(a,2,0);
	    q8 += Ix(a,0,2)*Iy(a,1,0)*Iz(a,2,0);
	    q9 += Ix(a,1,2)*Iy(a,1,0)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+0];
	//I[3] += q3*C[k+0];
	I[3] += q3*C_[0];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+0];
	//I[5] += q5*C[k+0];
	I[5] += q5*C_[0];
	//qK6 += q6*C[k+0];
	//I[6] += q6*C[k+0];
	I[6] += q6*C_[0];
	//qK7 += q7*C[k+0];
	//I[7] += q7*C[k+0];
	I[7] += q7*C_[0];
	//qK8 += q8*C[k+0];
	//I[8] += q8*C[k+0];
	I[8] += q8*C_[0];
	//qK9 += q9*C[k+0];
	//I[9] += q9*C[k+0];
	I[9] += q9*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[10]*NORMALIZE[4]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[11]*NORMALIZE[4]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[12]*NORMALIZE[4]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[13]*NORMALIZE[4]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[14]*NORMALIZE[4]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[15]*NORMALIZE[4]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[16]*NORMALIZE[4]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[17]*NORMALIZE[4]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[18]*NORMALIZE[4]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[19]*NORMALIZE[4]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK10 = ZERO;
     D128 qK12 = ZERO;
     D128 qK14 = ZERO;
     D128 qK16 = ZERO;
     D128 qK18 = ZERO;
#else
     //double qK10 = 0.0;
     //double qK11 = 0.0;
     //double qK12 = 0.0;
     //double qK13 = 0.0;
     //double qK14 = 0.0;
     //double qK15 = 0.0;
     //double qK16 = 0.0;
     //double qK17 = 0.0;
     //double qK18 = 0.0;
     //double qK19 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q10 = ZERO;
	D128 q11 = ZERO;
	D128 q12 = ZERO;
	D128 q13 = ZERO;
	D128 q14 = ZERO;
	D128 q15 = ZERO;
	D128 q16 = ZERO;
	D128 q17 = ZERO;
	D128 q18 = ZERO;
	D128 q19 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x30 = LOAD(&Ix(a,3,0));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y12 = LOAD(&Iy(a,1,2));
	    D128 y32 = LOAD(&Iy(a,3,2));
	    D128 y02 = LOAD(&Iy(a,0,2));
	    D128 y22 = LOAD(&Iy(a,2,2));
	    D128 z30 = LOAD(&Iz(a,3,0));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q10 = ADD(q10, MUL(MUL(x30, y02), z00));
	    q11 = ADD(q11, MUL(MUL(x00, y32), z00));
	    q12 = ADD(q12, MUL(MUL(x00, y02), z30));
	    q13 = ADD(q13, MUL(MUL(x20, y12), z00));
	    q14 = ADD(q14, MUL(MUL(x20, y02), z10));
	    q15 = ADD(q15, MUL(MUL(x10, y22), z00));
	    q16 = ADD(q16, MUL(MUL(x00, y22), z10));
	    q17 = ADD(q17, MUL(MUL(x10, y02), z20));
	    q18 = ADD(q18, MUL(MUL(x00, y12), z20));
	    q19 = ADD(q19, MUL(MUL(x10, y12), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x30 = LOAD1(&Ix(a,3,0));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y12 = LOAD1(&Iy(a,1,2));
	    D128 y32 = LOAD1(&Iy(a,3,2));
	    D128 y02 = LOAD1(&Iy(a,0,2));
	    D128 y22 = LOAD1(&Iy(a,2,2));
	    D128 z30 = LOAD1(&Iz(a,3,0));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q10 = ADD1(q10, MUL1(MUL1(x30, y02), z00));
	    q11 = ADD1(q11, MUL1(MUL1(x00, y32), z00));
	    q12 = ADD1(q12, MUL1(MUL1(x00, y02), z30));
	    q13 = ADD1(q13, MUL1(MUL1(x20, y12), z00));
	    q14 = ADD1(q14, MUL1(MUL1(x20, y02), z10));
	    q15 = ADD1(q15, MUL1(MUL1(x10, y22), z00));
	    q16 = ADD1(q16, MUL1(MUL1(x00, y22), z10));
	    q17 = ADD1(q17, MUL1(MUL1(x10, y02), z20));
	    q18 = ADD1(q18, MUL1(MUL1(x00, y12), z20));
	    q19 = ADD1(q19, MUL1(MUL1(x10, y12), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK10 = ADD(qK10, MUL(C00, HADD(q10, q11)));
	qK12 = ADD(qK12, MUL(C00, HADD(q12, q13)));
	qK14 = ADD(qK14, MUL(C00, HADD(q14, q15)));
	qK16 = ADD(qK16, MUL(C00, HADD(q16, q17)));
	qK18 = ADD(qK18, MUL(C00, HADD(q18, q19)));

#else // SSE
	    
	// function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
	T q16 = 0.0;
	T q17 = 0.0;
	T q18 = 0.0;
	T q19 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += Ix(a,3,0)*Iy(a,0,2)*Iz(a,0,0);
	    q11 += Ix(a,0,0)*Iy(a,3,2)*Iz(a,0,0);
	    q12 += Ix(a,0,0)*Iy(a,0,2)*Iz(a,3,0);
	    q13 += Ix(a,2,0)*Iy(a,1,2)*Iz(a,0,0);
	    q14 += Ix(a,2,0)*Iy(a,0,2)*Iz(a,1,0);
	    q15 += Ix(a,1,0)*Iy(a,2,2)*Iz(a,0,0);
	    q16 += Ix(a,0,0)*Iy(a,2,2)*Iz(a,1,0);
	    q17 += Ix(a,1,0)*Iy(a,0,2)*Iz(a,2,0);
	    q18 += Ix(a,0,0)*Iy(a,1,2)*Iz(a,2,0);
	    q19 += Ix(a,1,0)*Iy(a,1,2)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK10 += q10*C[k+0];
	//I[10] += q10*C[k+0];
	I[10] += q10*C_[0];
	//qK11 += q11*C[k+0];
	//I[11] += q11*C[k+0];
	I[11] += q11*C_[0];
	//qK12 += q12*C[k+0];
	//I[12] += q12*C[k+0];
	I[12] += q12*C_[0];
	//qK13 += q13*C[k+0];
	//I[13] += q13*C[k+0];
	I[13] += q13*C_[0];
	//qK14 += q14*C[k+0];
	//I[14] += q14*C[k+0];
	I[14] += q14*C_[0];
	//qK15 += q15*C[k+0];
	//I[15] += q15*C[k+0];
	I[15] += q15*C_[0];
	//qK16 += q16*C[k+0];
	//I[16] += q16*C[k+0];
	I[16] += q16*C_[0];
	//qK17 += q17*C[k+0];
	//I[17] += q17*C[k+0];
	I[17] += q17*C_[0];
	//qK18 += q18*C[k+0];
	//I[18] += q18*C[k+0];
	I[18] += q18*C_[0];
	//qK19 += q19*C[k+0];
	//I[19] += q19*C[k+0];
	I[19] += q19*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK10 = MUL(q, qK10);
	qK12 = MUL(q, qK12);
	qK14 = MUL(q, qK14);
	qK16 = MUL(q, qK16);
	qK18 = MUL(q, qK18);
	num += 10; //num += (fabs(I[18]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[10]) & 0xF) {
	    // 10
	    STOREU(&I[10], ADD(qK10, LOADU(&I[10])));
	    STOREU(&I[12], ADD(qK12, LOADU(&I[12])));
	    STOREU(&I[14], ADD(qK14, LOADU(&I[14])));
	    STOREU(&I[16], ADD(qK16, LOADU(&I[16])));
	    STOREU(&I[18], ADD(qK18, LOADU(&I[18])));
	}
	else {
	    STORE(&I[10], ADD(qK10, LOADU(&I[10])));
	    STORE(&I[12], ADD(qK12, LOADU(&I[12])));
	    STORE(&I[14], ADD(qK14, LOADU(&I[14])));
	    STORE(&I[16], ADD(qK16, LOADU(&I[16])));
	    STORE(&I[18], ADD(qK18, LOADU(&I[18])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[10] += scale*NORMALIZE[10]*NORMALIZE[5]*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*NORMALIZE[11]*NORMALIZE[5]*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*NORMALIZE[12]*NORMALIZE[5]*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*NORMALIZE[13]*NORMALIZE[5]*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*NORMALIZE[14]*NORMALIZE[5]*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*NORMALIZE[15]*NORMALIZE[5]*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*NORMALIZE[16]*NORMALIZE[5]*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*NORMALIZE[17]*NORMALIZE[5]*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*NORMALIZE[18]*NORMALIZE[5]*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*NORMALIZE[19]*NORMALIZE[5]*qK19;
	// num += (fabs(I[19]) >= tol);
    }
    else {
	// I[10] += scale*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*qK19;
	// num += (fabs(I[19]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK20 = ZERO;
     D128 qK22 = ZERO;
     D128 qK24 = ZERO;
     D128 qK26 = ZERO;
     D128 qK28 = ZERO;
#else
     //double qK20 = 0.0;
     //double qK21 = 0.0;
     //double qK22 = 0.0;
     //double qK23 = 0.0;
     //double qK24 = 0.0;
     //double qK25 = 0.0;
     //double qK26 = 0.0;
     //double qK27 = 0.0;
     //double qK28 = 0.0;
     //double qK29 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q20 = ZERO;
	D128 q21 = ZERO;
	D128 q22 = ZERO;
	D128 q23 = ZERO;
	D128 q24 = ZERO;
	D128 q25 = ZERO;
	D128 q26 = ZERO;
	D128 q27 = ZERO;
	D128 q28 = ZERO;
	D128 q29 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x30 = LOAD(&Ix(a,3,0));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y30 = LOAD(&Iy(a,3,0));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z12 = LOAD(&Iz(a,1,2));
	    D128 z32 = LOAD(&Iz(a,3,2));
	    D128 z02 = LOAD(&Iz(a,0,2));
	    D128 z22 = LOAD(&Iz(a,2,2));
	    q20 = ADD(q20, MUL(MUL(x30, y00), z02));
	    q21 = ADD(q21, MUL(MUL(x00, y30), z02));
	    q22 = ADD(q22, MUL(MUL(x00, y00), z32));
	    q23 = ADD(q23, MUL(MUL(x20, y10), z02));
	    q24 = ADD(q24, MUL(MUL(x20, y00), z12));
	    q25 = ADD(q25, MUL(MUL(x10, y20), z02));
	    q26 = ADD(q26, MUL(MUL(x00, y20), z12));
	    q27 = ADD(q27, MUL(MUL(x10, y00), z22));
	    q28 = ADD(q28, MUL(MUL(x00, y10), z22));
	    q29 = ADD(q29, MUL(MUL(x10, y10), z12));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x30 = LOAD1(&Ix(a,3,0));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y30 = LOAD1(&Iy(a,3,0));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z12 = LOAD1(&Iz(a,1,2));
	    D128 z32 = LOAD1(&Iz(a,3,2));
	    D128 z02 = LOAD1(&Iz(a,0,2));
	    D128 z22 = LOAD1(&Iz(a,2,2));
	    q20 = ADD1(q20, MUL1(MUL1(x30, y00), z02));
	    q21 = ADD1(q21, MUL1(MUL1(x00, y30), z02));
	    q22 = ADD1(q22, MUL1(MUL1(x00, y00), z32));
	    q23 = ADD1(q23, MUL1(MUL1(x20, y10), z02));
	    q24 = ADD1(q24, MUL1(MUL1(x20, y00), z12));
	    q25 = ADD1(q25, MUL1(MUL1(x10, y20), z02));
	    q26 = ADD1(q26, MUL1(MUL1(x00, y20), z12));
	    q27 = ADD1(q27, MUL1(MUL1(x10, y00), z22));
	    q28 = ADD1(q28, MUL1(MUL1(x00, y10), z22));
	    q29 = ADD1(q29, MUL1(MUL1(x10, y10), z12));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK20 = ADD(qK20, MUL(C00, HADD(q20, q21)));
	qK22 = ADD(qK22, MUL(C00, HADD(q22, q23)));
	qK24 = ADD(qK24, MUL(C00, HADD(q24, q25)));
	qK26 = ADD(qK26, MUL(C00, HADD(q26, q27)));
	qK28 = ADD(qK28, MUL(C00, HADD(q28, q29)));

#else // SSE
	    
	// function registers
	T q20 = 0.0;
	T q21 = 0.0;
	T q22 = 0.0;
	T q23 = 0.0;
	T q24 = 0.0;
	T q25 = 0.0;
	T q26 = 0.0;
	T q27 = 0.0;
	T q28 = 0.0;
	T q29 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q20 += Ix(a,3,0)*Iy(a,0,0)*Iz(a,0,2);
	    q21 += Ix(a,0,0)*Iy(a,3,0)*Iz(a,0,2);
	    q22 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,3,2);
	    q23 += Ix(a,2,0)*Iy(a,1,0)*Iz(a,0,2);
	    q24 += Ix(a,2,0)*Iy(a,0,0)*Iz(a,1,2);
	    q25 += Ix(a,1,0)*Iy(a,2,0)*Iz(a,0,2);
	    q26 += Ix(a,0,0)*Iy(a,2,0)*Iz(a,1,2);
	    q27 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,2,2);
	    q28 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,2,2);
	    q29 += Ix(a,1,0)*Iy(a,1,0)*Iz(a,1,2);
	}
	    
	//contraction coefficients
	//qK20 += q20*C[k+0];
	//I[20] += q20*C[k+0];
	I[20] += q20*C_[0];
	//qK21 += q21*C[k+0];
	//I[21] += q21*C[k+0];
	I[21] += q21*C_[0];
	//qK22 += q22*C[k+0];
	//I[22] += q22*C[k+0];
	I[22] += q22*C_[0];
	//qK23 += q23*C[k+0];
	//I[23] += q23*C[k+0];
	I[23] += q23*C_[0];
	//qK24 += q24*C[k+0];
	//I[24] += q24*C[k+0];
	I[24] += q24*C_[0];
	//qK25 += q25*C[k+0];
	//I[25] += q25*C[k+0];
	I[25] += q25*C_[0];
	//qK26 += q26*C[k+0];
	//I[26] += q26*C[k+0];
	I[26] += q26*C_[0];
	//qK27 += q27*C[k+0];
	//I[27] += q27*C[k+0];
	I[27] += q27*C_[0];
	//qK28 += q28*C[k+0];
	//I[28] += q28*C[k+0];
	I[28] += q28*C_[0];
	//qK29 += q29*C[k+0];
	//I[29] += q29*C[k+0];
	I[29] += q29*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK20 = MUL(q, qK20);
	qK22 = MUL(q, qK22);
	qK24 = MUL(q, qK24);
	qK26 = MUL(q, qK26);
	qK28 = MUL(q, qK28);
	num += 10; //num += (fabs(I[28]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[20]) & 0xF) {
	    // 20
	    STOREU(&I[20], ADD(qK20, LOADU(&I[20])));
	    STOREU(&I[22], ADD(qK22, LOADU(&I[22])));
	    STOREU(&I[24], ADD(qK24, LOADU(&I[24])));
	    STOREU(&I[26], ADD(qK26, LOADU(&I[26])));
	    STOREU(&I[28], ADD(qK28, LOADU(&I[28])));
	}
	else {
	    STORE(&I[20], ADD(qK20, LOADU(&I[20])));
	    STORE(&I[22], ADD(qK22, LOADU(&I[22])));
	    STORE(&I[24], ADD(qK24, LOADU(&I[24])));
	    STORE(&I[26], ADD(qK26, LOADU(&I[26])));
	    STORE(&I[28], ADD(qK28, LOADU(&I[28])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[20] += scale*NORMALIZE[10]*NORMALIZE[6]*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*NORMALIZE[11]*NORMALIZE[6]*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*NORMALIZE[12]*NORMALIZE[6]*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*NORMALIZE[13]*NORMALIZE[6]*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*NORMALIZE[14]*NORMALIZE[6]*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*NORMALIZE[15]*NORMALIZE[6]*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*NORMALIZE[16]*NORMALIZE[6]*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*NORMALIZE[17]*NORMALIZE[6]*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*NORMALIZE[18]*NORMALIZE[6]*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*NORMALIZE[19]*NORMALIZE[6]*qK29;
	// num += (fabs(I[29]) >= tol);
    }
    else {
	// I[20] += scale*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*qK29;
	// num += (fabs(I[29]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK30 = ZERO;
     D128 qK32 = ZERO;
     D128 qK34 = ZERO;
     D128 qK36 = ZERO;
     D128 qK38 = ZERO;
#else
     //double qK30 = 0.0;
     //double qK31 = 0.0;
     //double qK32 = 0.0;
     //double qK33 = 0.0;
     //double qK34 = 0.0;
     //double qK35 = 0.0;
     //double qK36 = 0.0;
     //double qK37 = 0.0;
     //double qK38 = 0.0;
     //double qK39 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q30 = ZERO;
	D128 q31 = ZERO;
	D128 q32 = ZERO;
	D128 q33 = ZERO;
	D128 q34 = ZERO;
	D128 q35 = ZERO;
	D128 q36 = ZERO;
	D128 q37 = ZERO;
	D128 q38 = ZERO;
	D128 q39 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x31 = LOAD(&Ix(a,3,1));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 x21 = LOAD(&Ix(a,2,1));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y31 = LOAD(&Iy(a,3,1));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 y21 = LOAD(&Iy(a,2,1));
	    D128 z30 = LOAD(&Iz(a,3,0));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q30 = ADD(q30, MUL(MUL(x31, y01), z00));
	    q31 = ADD(q31, MUL(MUL(x01, y31), z00));
	    q32 = ADD(q32, MUL(MUL(x01, y01), z30));
	    q33 = ADD(q33, MUL(MUL(x21, y11), z00));
	    q34 = ADD(q34, MUL(MUL(x21, y01), z10));
	    q35 = ADD(q35, MUL(MUL(x11, y21), z00));
	    q36 = ADD(q36, MUL(MUL(x01, y21), z10));
	    q37 = ADD(q37, MUL(MUL(x11, y01), z20));
	    q38 = ADD(q38, MUL(MUL(x01, y11), z20));
	    q39 = ADD(q39, MUL(MUL(x11, y11), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x31 = LOAD1(&Ix(a,3,1));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 x21 = LOAD1(&Ix(a,2,1));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y31 = LOAD1(&Iy(a,3,1));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 y21 = LOAD1(&Iy(a,2,1));
	    D128 z30 = LOAD1(&Iz(a,3,0));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q30 = ADD1(q30, MUL1(MUL1(x31, y01), z00));
	    q31 = ADD1(q31, MUL1(MUL1(x01, y31), z00));
	    q32 = ADD1(q32, MUL1(MUL1(x01, y01), z30));
	    q33 = ADD1(q33, MUL1(MUL1(x21, y11), z00));
	    q34 = ADD1(q34, MUL1(MUL1(x21, y01), z10));
	    q35 = ADD1(q35, MUL1(MUL1(x11, y21), z00));
	    q36 = ADD1(q36, MUL1(MUL1(x01, y21), z10));
	    q37 = ADD1(q37, MUL1(MUL1(x11, y01), z20));
	    q38 = ADD1(q38, MUL1(MUL1(x01, y11), z20));
	    q39 = ADD1(q39, MUL1(MUL1(x11, y11), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK30 = ADD(qK30, MUL(C00, HADD(q30, q31)));
	qK32 = ADD(qK32, MUL(C00, HADD(q32, q33)));
	qK34 = ADD(qK34, MUL(C00, HADD(q34, q35)));
	qK36 = ADD(qK36, MUL(C00, HADD(q36, q37)));
	qK38 = ADD(qK38, MUL(C00, HADD(q38, q39)));

#else // SSE
	    
	// function registers
	T q30 = 0.0;
	T q31 = 0.0;
	T q32 = 0.0;
	T q33 = 0.0;
	T q34 = 0.0;
	T q35 = 0.0;
	T q36 = 0.0;
	T q37 = 0.0;
	T q38 = 0.0;
	T q39 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q30 += Ix(a,3,1)*Iy(a,0,1)*Iz(a,0,0);
	    q31 += Ix(a,0,1)*Iy(a,3,1)*Iz(a,0,0);
	    q32 += Ix(a,0,1)*Iy(a,0,1)*Iz(a,3,0);
	    q33 += Ix(a,2,1)*Iy(a,1,1)*Iz(a,0,0);
	    q34 += Ix(a,2,1)*Iy(a,0,1)*Iz(a,1,0);
	    q35 += Ix(a,1,1)*Iy(a,2,1)*Iz(a,0,0);
	    q36 += Ix(a,0,1)*Iy(a,2,1)*Iz(a,1,0);
	    q37 += Ix(a,1,1)*Iy(a,0,1)*Iz(a,2,0);
	    q38 += Ix(a,0,1)*Iy(a,1,1)*Iz(a,2,0);
	    q39 += Ix(a,1,1)*Iy(a,1,1)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK30 += q30*C[k+0];
	//I[30] += q30*C[k+0];
	I[30] += q30*C_[0];
	//qK31 += q31*C[k+0];
	//I[31] += q31*C[k+0];
	I[31] += q31*C_[0];
	//qK32 += q32*C[k+0];
	//I[32] += q32*C[k+0];
	I[32] += q32*C_[0];
	//qK33 += q33*C[k+0];
	//I[33] += q33*C[k+0];
	I[33] += q33*C_[0];
	//qK34 += q34*C[k+0];
	//I[34] += q34*C[k+0];
	I[34] += q34*C_[0];
	//qK35 += q35*C[k+0];
	//I[35] += q35*C[k+0];
	I[35] += q35*C_[0];
	//qK36 += q36*C[k+0];
	//I[36] += q36*C[k+0];
	I[36] += q36*C_[0];
	//qK37 += q37*C[k+0];
	//I[37] += q37*C[k+0];
	I[37] += q37*C_[0];
	//qK38 += q38*C[k+0];
	//I[38] += q38*C[k+0];
	I[38] += q38*C_[0];
	//qK39 += q39*C[k+0];
	//I[39] += q39*C[k+0];
	I[39] += q39*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK30 = MUL(q, qK30);
	qK32 = MUL(q, qK32);
	qK34 = MUL(q, qK34);
	qK36 = MUL(q, qK36);
	qK38 = MUL(q, qK38);
	num += 10; //num += (fabs(I[38]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[30]) & 0xF) {
	    // 30
	    STOREU(&I[30], ADD(qK30, LOADU(&I[30])));
	    STOREU(&I[32], ADD(qK32, LOADU(&I[32])));
	    STOREU(&I[34], ADD(qK34, LOADU(&I[34])));
	    STOREU(&I[36], ADD(qK36, LOADU(&I[36])));
	    STOREU(&I[38], ADD(qK38, LOADU(&I[38])));
	}
	else {
	    STORE(&I[30], ADD(qK30, LOADU(&I[30])));
	    STORE(&I[32], ADD(qK32, LOADU(&I[32])));
	    STORE(&I[34], ADD(qK34, LOADU(&I[34])));
	    STORE(&I[36], ADD(qK36, LOADU(&I[36])));
	    STORE(&I[38], ADD(qK38, LOADU(&I[38])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[30] += scale*NORMALIZE[10]*NORMALIZE[7]*qK30;
	// num += (fabs(I[30]) >= tol);
	// I[31] += scale*NORMALIZE[11]*NORMALIZE[7]*qK31;
	// num += (fabs(I[31]) >= tol);
	// I[32] += scale*NORMALIZE[12]*NORMALIZE[7]*qK32;
	// num += (fabs(I[32]) >= tol);
	// I[33] += scale*NORMALIZE[13]*NORMALIZE[7]*qK33;
	// num += (fabs(I[33]) >= tol);
	// I[34] += scale*NORMALIZE[14]*NORMALIZE[7]*qK34;
	// num += (fabs(I[34]) >= tol);
	// I[35] += scale*NORMALIZE[15]*NORMALIZE[7]*qK35;
	// num += (fabs(I[35]) >= tol);
	// I[36] += scale*NORMALIZE[16]*NORMALIZE[7]*qK36;
	// num += (fabs(I[36]) >= tol);
	// I[37] += scale*NORMALIZE[17]*NORMALIZE[7]*qK37;
	// num += (fabs(I[37]) >= tol);
	// I[38] += scale*NORMALIZE[18]*NORMALIZE[7]*qK38;
	// num += (fabs(I[38]) >= tol);
	// I[39] += scale*NORMALIZE[19]*NORMALIZE[7]*qK39;
	// num += (fabs(I[39]) >= tol);
    }
    else {
	// I[30] += scale*qK30;
	// num += (fabs(I[30]) >= tol);
	// I[31] += scale*qK31;
	// num += (fabs(I[31]) >= tol);
	// I[32] += scale*qK32;
	// num += (fabs(I[32]) >= tol);
	// I[33] += scale*qK33;
	// num += (fabs(I[33]) >= tol);
	// I[34] += scale*qK34;
	// num += (fabs(I[34]) >= tol);
	// I[35] += scale*qK35;
	// num += (fabs(I[35]) >= tol);
	// I[36] += scale*qK36;
	// num += (fabs(I[36]) >= tol);
	// I[37] += scale*qK37;
	// num += (fabs(I[37]) >= tol);
	// I[38] += scale*qK38;
	// num += (fabs(I[38]) >= tol);
	// I[39] += scale*qK39;
	// num += (fabs(I[39]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK40 = ZERO;
     D128 qK42 = ZERO;
     D128 qK44 = ZERO;
     D128 qK46 = ZERO;
     D128 qK48 = ZERO;
#else
     //double qK40 = 0.0;
     //double qK41 = 0.0;
     //double qK42 = 0.0;
     //double qK43 = 0.0;
     //double qK44 = 0.0;
     //double qK45 = 0.0;
     //double qK46 = 0.0;
     //double qK47 = 0.0;
     //double qK48 = 0.0;
     //double qK49 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q40 = ZERO;
	D128 q41 = ZERO;
	D128 q42 = ZERO;
	D128 q43 = ZERO;
	D128 q44 = ZERO;
	D128 q45 = ZERO;
	D128 q46 = ZERO;
	D128 q47 = ZERO;
	D128 q48 = ZERO;
	D128 q49 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x31 = LOAD(&Ix(a,3,1));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 x21 = LOAD(&Ix(a,2,1));
	    D128 y30 = LOAD(&Iy(a,3,0));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z31 = LOAD(&Iz(a,3,1));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    D128 z21 = LOAD(&Iz(a,2,1));
	    q40 = ADD(q40, MUL(MUL(x31, y00), z01));
	    q41 = ADD(q41, MUL(MUL(x01, y30), z01));
	    q42 = ADD(q42, MUL(MUL(x01, y00), z31));
	    q43 = ADD(q43, MUL(MUL(x21, y10), z01));
	    q44 = ADD(q44, MUL(MUL(x21, y00), z11));
	    q45 = ADD(q45, MUL(MUL(x11, y20), z01));
	    q46 = ADD(q46, MUL(MUL(x01, y20), z11));
	    q47 = ADD(q47, MUL(MUL(x11, y00), z21));
	    q48 = ADD(q48, MUL(MUL(x01, y10), z21));
	    q49 = ADD(q49, MUL(MUL(x11, y10), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x31 = LOAD1(&Ix(a,3,1));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 x21 = LOAD1(&Ix(a,2,1));
	    D128 y30 = LOAD1(&Iy(a,3,0));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z31 = LOAD1(&Iz(a,3,1));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    D128 z21 = LOAD1(&Iz(a,2,1));
	    q40 = ADD1(q40, MUL1(MUL1(x31, y00), z01));
	    q41 = ADD1(q41, MUL1(MUL1(x01, y30), z01));
	    q42 = ADD1(q42, MUL1(MUL1(x01, y00), z31));
	    q43 = ADD1(q43, MUL1(MUL1(x21, y10), z01));
	    q44 = ADD1(q44, MUL1(MUL1(x21, y00), z11));
	    q45 = ADD1(q45, MUL1(MUL1(x11, y20), z01));
	    q46 = ADD1(q46, MUL1(MUL1(x01, y20), z11));
	    q47 = ADD1(q47, MUL1(MUL1(x11, y00), z21));
	    q48 = ADD1(q48, MUL1(MUL1(x01, y10), z21));
	    q49 = ADD1(q49, MUL1(MUL1(x11, y10), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK40 = ADD(qK40, MUL(C00, HADD(q40, q41)));
	qK42 = ADD(qK42, MUL(C00, HADD(q42, q43)));
	qK44 = ADD(qK44, MUL(C00, HADD(q44, q45)));
	qK46 = ADD(qK46, MUL(C00, HADD(q46, q47)));
	qK48 = ADD(qK48, MUL(C00, HADD(q48, q49)));

#else // SSE
	    
	// function registers
	T q40 = 0.0;
	T q41 = 0.0;
	T q42 = 0.0;
	T q43 = 0.0;
	T q44 = 0.0;
	T q45 = 0.0;
	T q46 = 0.0;
	T q47 = 0.0;
	T q48 = 0.0;
	T q49 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q40 += Ix(a,3,1)*Iy(a,0,0)*Iz(a,0,1);
	    q41 += Ix(a,0,1)*Iy(a,3,0)*Iz(a,0,1);
	    q42 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,3,1);
	    q43 += Ix(a,2,1)*Iy(a,1,0)*Iz(a,0,1);
	    q44 += Ix(a,2,1)*Iy(a,0,0)*Iz(a,1,1);
	    q45 += Ix(a,1,1)*Iy(a,2,0)*Iz(a,0,1);
	    q46 += Ix(a,0,1)*Iy(a,2,0)*Iz(a,1,1);
	    q47 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,2,1);
	    q48 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,2,1);
	    q49 += Ix(a,1,1)*Iy(a,1,0)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK40 += q40*C[k+0];
	//I[40] += q40*C[k+0];
	I[40] += q40*C_[0];
	//qK41 += q41*C[k+0];
	//I[41] += q41*C[k+0];
	I[41] += q41*C_[0];
	//qK42 += q42*C[k+0];
	//I[42] += q42*C[k+0];
	I[42] += q42*C_[0];
	//qK43 += q43*C[k+0];
	//I[43] += q43*C[k+0];
	I[43] += q43*C_[0];
	//qK44 += q44*C[k+0];
	//I[44] += q44*C[k+0];
	I[44] += q44*C_[0];
	//qK45 += q45*C[k+0];
	//I[45] += q45*C[k+0];
	I[45] += q45*C_[0];
	//qK46 += q46*C[k+0];
	//I[46] += q46*C[k+0];
	I[46] += q46*C_[0];
	//qK47 += q47*C[k+0];
	//I[47] += q47*C[k+0];
	I[47] += q47*C_[0];
	//qK48 += q48*C[k+0];
	//I[48] += q48*C[k+0];
	I[48] += q48*C_[0];
	//qK49 += q49*C[k+0];
	//I[49] += q49*C[k+0];
	I[49] += q49*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK40 = MUL(q, qK40);
	qK42 = MUL(q, qK42);
	qK44 = MUL(q, qK44);
	qK46 = MUL(q, qK46);
	qK48 = MUL(q, qK48);
	num += 10; //num += (fabs(I[48]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[40]) & 0xF) {
	    // 40
	    STOREU(&I[40], ADD(qK40, LOADU(&I[40])));
	    STOREU(&I[42], ADD(qK42, LOADU(&I[42])));
	    STOREU(&I[44], ADD(qK44, LOADU(&I[44])));
	    STOREU(&I[46], ADD(qK46, LOADU(&I[46])));
	    STOREU(&I[48], ADD(qK48, LOADU(&I[48])));
	}
	else {
	    STORE(&I[40], ADD(qK40, LOADU(&I[40])));
	    STORE(&I[42], ADD(qK42, LOADU(&I[42])));
	    STORE(&I[44], ADD(qK44, LOADU(&I[44])));
	    STORE(&I[46], ADD(qK46, LOADU(&I[46])));
	    STORE(&I[48], ADD(qK48, LOADU(&I[48])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[40] += scale*NORMALIZE[10]*NORMALIZE[8]*qK40;
	// num += (fabs(I[40]) >= tol);
	// I[41] += scale*NORMALIZE[11]*NORMALIZE[8]*qK41;
	// num += (fabs(I[41]) >= tol);
	// I[42] += scale*NORMALIZE[12]*NORMALIZE[8]*qK42;
	// num += (fabs(I[42]) >= tol);
	// I[43] += scale*NORMALIZE[13]*NORMALIZE[8]*qK43;
	// num += (fabs(I[43]) >= tol);
	// I[44] += scale*NORMALIZE[14]*NORMALIZE[8]*qK44;
	// num += (fabs(I[44]) >= tol);
	// I[45] += scale*NORMALIZE[15]*NORMALIZE[8]*qK45;
	// num += (fabs(I[45]) >= tol);
	// I[46] += scale*NORMALIZE[16]*NORMALIZE[8]*qK46;
	// num += (fabs(I[46]) >= tol);
	// I[47] += scale*NORMALIZE[17]*NORMALIZE[8]*qK47;
	// num += (fabs(I[47]) >= tol);
	// I[48] += scale*NORMALIZE[18]*NORMALIZE[8]*qK48;
	// num += (fabs(I[48]) >= tol);
	// I[49] += scale*NORMALIZE[19]*NORMALIZE[8]*qK49;
	// num += (fabs(I[49]) >= tol);
    }
    else {
	// I[40] += scale*qK40;
	// num += (fabs(I[40]) >= tol);
	// I[41] += scale*qK41;
	// num += (fabs(I[41]) >= tol);
	// I[42] += scale*qK42;
	// num += (fabs(I[42]) >= tol);
	// I[43] += scale*qK43;
	// num += (fabs(I[43]) >= tol);
	// I[44] += scale*qK44;
	// num += (fabs(I[44]) >= tol);
	// I[45] += scale*qK45;
	// num += (fabs(I[45]) >= tol);
	// I[46] += scale*qK46;
	// num += (fabs(I[46]) >= tol);
	// I[47] += scale*qK47;
	// num += (fabs(I[47]) >= tol);
	// I[48] += scale*qK48;
	// num += (fabs(I[48]) >= tol);
	// I[49] += scale*qK49;
	// num += (fabs(I[49]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK50 = ZERO;
     D128 qK52 = ZERO;
     D128 qK54 = ZERO;
     D128 qK56 = ZERO;
     D128 qK58 = ZERO;
#else
     //double qK50 = 0.0;
     //double qK51 = 0.0;
     //double qK52 = 0.0;
     //double qK53 = 0.0;
     //double qK54 = 0.0;
     //double qK55 = 0.0;
     //double qK56 = 0.0;
     //double qK57 = 0.0;
     //double qK58 = 0.0;
     //double qK59 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q50 = ZERO;
	D128 q51 = ZERO;
	D128 q52 = ZERO;
	D128 q53 = ZERO;
	D128 q54 = ZERO;
	D128 q55 = ZERO;
	D128 q56 = ZERO;
	D128 q57 = ZERO;
	D128 q58 = ZERO;
	D128 q59 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x30 = LOAD(&Ix(a,3,0));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y31 = LOAD(&Iy(a,3,1));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 y21 = LOAD(&Iy(a,2,1));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z31 = LOAD(&Iz(a,3,1));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    D128 z21 = LOAD(&Iz(a,2,1));
	    q50 = ADD(q50, MUL(MUL(x30, y01), z01));
	    q51 = ADD(q51, MUL(MUL(x00, y31), z01));
	    q52 = ADD(q52, MUL(MUL(x00, y01), z31));
	    q53 = ADD(q53, MUL(MUL(x20, y11), z01));
	    q54 = ADD(q54, MUL(MUL(x20, y01), z11));
	    q55 = ADD(q55, MUL(MUL(x10, y21), z01));
	    q56 = ADD(q56, MUL(MUL(x00, y21), z11));
	    q57 = ADD(q57, MUL(MUL(x10, y01), z21));
	    q58 = ADD(q58, MUL(MUL(x00, y11), z21));
	    q59 = ADD(q59, MUL(MUL(x10, y11), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x30 = LOAD1(&Ix(a,3,0));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y31 = LOAD1(&Iy(a,3,1));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 y21 = LOAD1(&Iy(a,2,1));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z31 = LOAD1(&Iz(a,3,1));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    D128 z21 = LOAD1(&Iz(a,2,1));
	    q50 = ADD1(q50, MUL1(MUL1(x30, y01), z01));
	    q51 = ADD1(q51, MUL1(MUL1(x00, y31), z01));
	    q52 = ADD1(q52, MUL1(MUL1(x00, y01), z31));
	    q53 = ADD1(q53, MUL1(MUL1(x20, y11), z01));
	    q54 = ADD1(q54, MUL1(MUL1(x20, y01), z11));
	    q55 = ADD1(q55, MUL1(MUL1(x10, y21), z01));
	    q56 = ADD1(q56, MUL1(MUL1(x00, y21), z11));
	    q57 = ADD1(q57, MUL1(MUL1(x10, y01), z21));
	    q58 = ADD1(q58, MUL1(MUL1(x00, y11), z21));
	    q59 = ADD1(q59, MUL1(MUL1(x10, y11), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK50 = ADD(qK50, MUL(C00, HADD(q50, q51)));
	qK52 = ADD(qK52, MUL(C00, HADD(q52, q53)));
	qK54 = ADD(qK54, MUL(C00, HADD(q54, q55)));
	qK56 = ADD(qK56, MUL(C00, HADD(q56, q57)));
	qK58 = ADD(qK58, MUL(C00, HADD(q58, q59)));

#else // SSE
	    
	// function registers
	T q50 = 0.0;
	T q51 = 0.0;
	T q52 = 0.0;
	T q53 = 0.0;
	T q54 = 0.0;
	T q55 = 0.0;
	T q56 = 0.0;
	T q57 = 0.0;
	T q58 = 0.0;
	T q59 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q50 += Ix(a,3,0)*Iy(a,0,1)*Iz(a,0,1);
	    q51 += Ix(a,0,0)*Iy(a,3,1)*Iz(a,0,1);
	    q52 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,3,1);
	    q53 += Ix(a,2,0)*Iy(a,1,1)*Iz(a,0,1);
	    q54 += Ix(a,2,0)*Iy(a,0,1)*Iz(a,1,1);
	    q55 += Ix(a,1,0)*Iy(a,2,1)*Iz(a,0,1);
	    q56 += Ix(a,0,0)*Iy(a,2,1)*Iz(a,1,1);
	    q57 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,2,1);
	    q58 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,2,1);
	    q59 += Ix(a,1,0)*Iy(a,1,1)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK50 += q50*C[k+0];
	//I[50] += q50*C[k+0];
	I[50] += q50*C_[0];
	//qK51 += q51*C[k+0];
	//I[51] += q51*C[k+0];
	I[51] += q51*C_[0];
	//qK52 += q52*C[k+0];
	//I[52] += q52*C[k+0];
	I[52] += q52*C_[0];
	//qK53 += q53*C[k+0];
	//I[53] += q53*C[k+0];
	I[53] += q53*C_[0];
	//qK54 += q54*C[k+0];
	//I[54] += q54*C[k+0];
	I[54] += q54*C_[0];
	//qK55 += q55*C[k+0];
	//I[55] += q55*C[k+0];
	I[55] += q55*C_[0];
	//qK56 += q56*C[k+0];
	//I[56] += q56*C[k+0];
	I[56] += q56*C_[0];
	//qK57 += q57*C[k+0];
	//I[57] += q57*C[k+0];
	I[57] += q57*C_[0];
	//qK58 += q58*C[k+0];
	//I[58] += q58*C[k+0];
	I[58] += q58*C_[0];
	//qK59 += q59*C[k+0];
	//I[59] += q59*C[k+0];
	I[59] += q59*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK50 = MUL(q, qK50);
	qK52 = MUL(q, qK52);
	qK54 = MUL(q, qK54);
	qK56 = MUL(q, qK56);
	qK58 = MUL(q, qK58);
	num += 10; //num += (fabs(I[58]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[50]) & 0xF) {
	    // 50
	    STOREU(&I[50], ADD(qK50, LOADU(&I[50])));
	    STOREU(&I[52], ADD(qK52, LOADU(&I[52])));
	    STOREU(&I[54], ADD(qK54, LOADU(&I[54])));
	    STOREU(&I[56], ADD(qK56, LOADU(&I[56])));
	    STOREU(&I[58], ADD(qK58, LOADU(&I[58])));
	}
	else {
	    STORE(&I[50], ADD(qK50, LOADU(&I[50])));
	    STORE(&I[52], ADD(qK52, LOADU(&I[52])));
	    STORE(&I[54], ADD(qK54, LOADU(&I[54])));
	    STORE(&I[56], ADD(qK56, LOADU(&I[56])));
	    STORE(&I[58], ADD(qK58, LOADU(&I[58])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[50] += scale*NORMALIZE[10]*NORMALIZE[9]*qK50;
	// num += (fabs(I[50]) >= tol);
	// I[51] += scale*NORMALIZE[11]*NORMALIZE[9]*qK51;
	// num += (fabs(I[51]) >= tol);
	// I[52] += scale*NORMALIZE[12]*NORMALIZE[9]*qK52;
	// num += (fabs(I[52]) >= tol);
	// I[53] += scale*NORMALIZE[13]*NORMALIZE[9]*qK53;
	// num += (fabs(I[53]) >= tol);
	// I[54] += scale*NORMALIZE[14]*NORMALIZE[9]*qK54;
	// num += (fabs(I[54]) >= tol);
	// I[55] += scale*NORMALIZE[15]*NORMALIZE[9]*qK55;
	// num += (fabs(I[55]) >= tol);
	// I[56] += scale*NORMALIZE[16]*NORMALIZE[9]*qK56;
	// num += (fabs(I[56]) >= tol);
	// I[57] += scale*NORMALIZE[17]*NORMALIZE[9]*qK57;
	// num += (fabs(I[57]) >= tol);
	// I[58] += scale*NORMALIZE[18]*NORMALIZE[9]*qK58;
	// num += (fabs(I[58]) >= tol);
	// I[59] += scale*NORMALIZE[19]*NORMALIZE[9]*qK59;
	// num += (fabs(I[59]) >= tol);
    }
    else {
	// I[50] += scale*qK50;
	// num += (fabs(I[50]) >= tol);
	// I[51] += scale*qK51;
	// num += (fabs(I[51]) >= tol);
	// I[52] += scale*qK52;
	// num += (fabs(I[52]) >= tol);
	// I[53] += scale*qK53;
	// num += (fabs(I[53]) >= tol);
	// I[54] += scale*qK54;
	// num += (fabs(I[54]) >= tol);
	// I[55] += scale*qK55;
	// num += (fabs(I[55]) >= tol);
	// I[56] += scale*qK56;
	// num += (fabs(I[56]) >= tol);
	// I[57] += scale*qK57;
	// num += (fabs(I[57]) >= tol);
	// I[58] += scale*qK58;
	// num += (fabs(I[58]) >= tol);
	// I[59] += scale*qK59;
	// num += (fabs(I[59]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <spd| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::SP,rysq::D> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::SP,rysq::D> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 2;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x12 = LOAD(&Ix(a,1,2));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x02 = LOAD(&Ix(a,0,2));
	    D128 y12 = LOAD(&Iy(a,1,2));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y02 = LOAD(&Iy(a,0,2));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z02 = LOAD(&Iz(a,0,2));
	    q0 = ADD(q0, MUL(MUL(x02, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x12, y00), z00));
	    q2 = ADD(q2, MUL(MUL(x02, y10), z00));
	    q3 = ADD(q3, MUL(MUL(x02, y00), z10));
	    q4 = ADD(q4, MUL(MUL(x00, y02), z00));
	    q5 = ADD(q5, MUL(MUL(x10, y02), z00));
	    q6 = ADD(q6, MUL(MUL(x00, y12), z00));
	    q7 = ADD(q7, MUL(MUL(x00, y02), z10));
	    q8 = ADD(q8, MUL(MUL(x00, y00), z02));
	    q9 = ADD(q9, MUL(MUL(x10, y00), z02));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x12 = LOAD1(&Ix(a,1,2));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x02 = LOAD1(&Ix(a,0,2));
	    D128 y12 = LOAD1(&Iy(a,1,2));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y02 = LOAD1(&Iy(a,0,2));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z02 = LOAD1(&Iz(a,0,2));
	    q0 = ADD1(q0, MUL1(MUL1(x02, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x12, y00), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x02, y10), z00));
	    q3 = ADD1(q3, MUL1(MUL1(x02, y00), z10));
	    q4 = ADD1(q4, MUL1(MUL1(x00, y02), z00));
	    q5 = ADD1(q5, MUL1(MUL1(x10, y02), z00));
	    q6 = ADD1(q6, MUL1(MUL1(x00, y12), z00));
	    q7 = ADD1(q7, MUL1(MUL1(x00, y02), z10));
	    q8 = ADD1(q8, MUL1(MUL1(x00, y00), z02));
	    q9 = ADD1(q9, MUL1(MUL1(x10, y00), z02));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C01 = LOADU(&C[k+0]);
	qK0 = ADD(qK0, MUL(C01, HADD(q0, q1)));
	D128 C11 = LOADDUP(&C[k+1]);
	qK2 = ADD(qK2, MUL(C11, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C01, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C11, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C01, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,0,2)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,1,2)*Iy(a,0,0)*Iz(a,0,0);
	    q2 += Ix(a,0,2)*Iy(a,1,0)*Iz(a,0,0);
	    q3 += Ix(a,0,2)*Iy(a,0,0)*Iz(a,1,0);
	    q4 += Ix(a,0,0)*Iy(a,0,2)*Iz(a,0,0);
	    q5 += Ix(a,1,0)*Iy(a,0,2)*Iz(a,0,0);
	    q6 += Ix(a,0,0)*Iy(a,1,2)*Iz(a,0,0);
	    q7 += Ix(a,0,0)*Iy(a,0,2)*Iz(a,1,0);
	    q8 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,0,2);
	    q9 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,0,2);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+1];
	//I[1] += q1*C[k+1];
	I[1] += q1*C_[1];
	//qK2 += q2*C[k+1];
	//I[2] += q2*C[k+1];
	I[2] += q2*C_[1];
	//qK3 += q3*C[k+1];
	//I[3] += q3*C[k+1];
	I[3] += q3*C_[1];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+1];
	//I[5] += q5*C[k+1];
	I[5] += q5*C_[1];
	//qK6 += q6*C[k+1];
	//I[6] += q6*C[k+1];
	I[6] += q6*C_[1];
	//qK7 += q7*C[k+1];
	//I[7] += q7*C[k+1];
	I[7] += q7*C_[1];
	//qK8 += q8*C[k+0];
	//I[8] += q8*C[k+0];
	I[8] += q8*C_[0];
	//qK9 += q9*C[k+1];
	//I[9] += q9*C[k+1];
	I[9] += q9*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[0]*NORMALIZE[4]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[1]*NORMALIZE[4]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[2]*NORMALIZE[4]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[3]*NORMALIZE[4]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[0]*NORMALIZE[5]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[1]*NORMALIZE[5]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[2]*NORMALIZE[5]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[3]*NORMALIZE[5]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[0]*NORMALIZE[6]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[1]*NORMALIZE[6]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK10 = ZERO;
     D128 qK12 = ZERO;
     D128 qK14 = ZERO;
     D128 qK16 = ZERO;
     D128 qK18 = ZERO;
#else
     //double qK10 = 0.0;
     //double qK11 = 0.0;
     //double qK12 = 0.0;
     //double qK13 = 0.0;
     //double qK14 = 0.0;
     //double qK15 = 0.0;
     //double qK16 = 0.0;
     //double qK17 = 0.0;
     //double qK18 = 0.0;
     //double qK19 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q10 = ZERO;
	D128 q11 = ZERO;
	D128 q12 = ZERO;
	D128 q13 = ZERO;
	D128 q14 = ZERO;
	D128 q15 = ZERO;
	D128 q16 = ZERO;
	D128 q17 = ZERO;
	D128 q18 = ZERO;
	D128 q19 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z12 = LOAD(&Iz(a,1,2));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z02 = LOAD(&Iz(a,0,2));
	    q10 = ADD(q10, MUL(MUL(x00, y10), z02));
	    q11 = ADD(q11, MUL(MUL(x00, y00), z12));
	    q12 = ADD(q12, MUL(MUL(x01, y01), z00));
	    q13 = ADD(q13, MUL(MUL(x11, y01), z00));
	    q14 = ADD(q14, MUL(MUL(x01, y11), z00));
	    q15 = ADD(q15, MUL(MUL(x01, y01), z10));
	    q16 = ADD(q16, MUL(MUL(x01, y00), z01));
	    q17 = ADD(q17, MUL(MUL(x11, y00), z01));
	    q18 = ADD(q18, MUL(MUL(x01, y10), z01));
	    q19 = ADD(q19, MUL(MUL(x01, y00), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z12 = LOAD1(&Iz(a,1,2));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z02 = LOAD1(&Iz(a,0,2));
	    q10 = ADD1(q10, MUL1(MUL1(x00, y10), z02));
	    q11 = ADD1(q11, MUL1(MUL1(x00, y00), z12));
	    q12 = ADD1(q12, MUL1(MUL1(x01, y01), z00));
	    q13 = ADD1(q13, MUL1(MUL1(x11, y01), z00));
	    q14 = ADD1(q14, MUL1(MUL1(x01, y11), z00));
	    q15 = ADD1(q15, MUL1(MUL1(x01, y01), z10));
	    q16 = ADD1(q16, MUL1(MUL1(x01, y00), z01));
	    q17 = ADD1(q17, MUL1(MUL1(x11, y00), z01));
	    q18 = ADD1(q18, MUL1(MUL1(x01, y10), z01));
	    q19 = ADD1(q19, MUL1(MUL1(x01, y00), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C11 = LOADDUP(&C[k+1]);
	qK10 = ADD(qK10, MUL(C11, HADD(q10, q11)));
	D128 C01 = LOADU(&C[k+0]);
	qK12 = ADD(qK12, MUL(C01, HADD(q12, q13)));
	qK14 = ADD(qK14, MUL(C11, HADD(q14, q15)));
	qK16 = ADD(qK16, MUL(C01, HADD(q16, q17)));
	qK18 = ADD(qK18, MUL(C11, HADD(q18, q19)));

#else // SSE
	    
	// function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
	T q16 = 0.0;
	T q17 = 0.0;
	T q18 = 0.0;
	T q19 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,0,2);
	    q11 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,1,2);
	    q12 += Ix(a,0,1)*Iy(a,0,1)*Iz(a,0,0);
	    q13 += Ix(a,1,1)*Iy(a,0,1)*Iz(a,0,0);
	    q14 += Ix(a,0,1)*Iy(a,1,1)*Iz(a,0,0);
	    q15 += Ix(a,0,1)*Iy(a,0,1)*Iz(a,1,0);
	    q16 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,0,1);
	    q17 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,0,1);
	    q18 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,0,1);
	    q19 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK10 += q10*C[k+1];
	//I[10] += q10*C[k+1];
	I[10] += q10*C_[1];
	//qK11 += q11*C[k+1];
	//I[11] += q11*C[k+1];
	I[11] += q11*C_[1];
	//qK12 += q12*C[k+0];
	//I[12] += q12*C[k+0];
	I[12] += q12*C_[0];
	//qK13 += q13*C[k+1];
	//I[13] += q13*C[k+1];
	I[13] += q13*C_[1];
	//qK14 += q14*C[k+1];
	//I[14] += q14*C[k+1];
	I[14] += q14*C_[1];
	//qK15 += q15*C[k+1];
	//I[15] += q15*C[k+1];
	I[15] += q15*C_[1];
	//qK16 += q16*C[k+0];
	//I[16] += q16*C[k+0];
	I[16] += q16*C_[0];
	//qK17 += q17*C[k+1];
	//I[17] += q17*C[k+1];
	I[17] += q17*C_[1];
	//qK18 += q18*C[k+1];
	//I[18] += q18*C[k+1];
	I[18] += q18*C_[1];
	//qK19 += q19*C[k+1];
	//I[19] += q19*C[k+1];
	I[19] += q19*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK10 = MUL(q, qK10);
	qK12 = MUL(q, qK12);
	qK14 = MUL(q, qK14);
	qK16 = MUL(q, qK16);
	qK18 = MUL(q, qK18);
	num += 10; //num += (fabs(I[18]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[10]) & 0xF) {
	    // 10
	    STOREU(&I[10], ADD(qK10, LOADU(&I[10])));
	    STOREU(&I[12], ADD(qK12, LOADU(&I[12])));
	    STOREU(&I[14], ADD(qK14, LOADU(&I[14])));
	    STOREU(&I[16], ADD(qK16, LOADU(&I[16])));
	    STOREU(&I[18], ADD(qK18, LOADU(&I[18])));
	}
	else {
	    STORE(&I[10], ADD(qK10, LOADU(&I[10])));
	    STORE(&I[12], ADD(qK12, LOADU(&I[12])));
	    STORE(&I[14], ADD(qK14, LOADU(&I[14])));
	    STORE(&I[16], ADD(qK16, LOADU(&I[16])));
	    STORE(&I[18], ADD(qK18, LOADU(&I[18])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[10] += scale*NORMALIZE[2]*NORMALIZE[6]*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*NORMALIZE[3]*NORMALIZE[6]*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*NORMALIZE[0]*NORMALIZE[7]*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*NORMALIZE[1]*NORMALIZE[7]*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*NORMALIZE[2]*NORMALIZE[7]*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*NORMALIZE[3]*NORMALIZE[7]*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*NORMALIZE[0]*NORMALIZE[8]*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*NORMALIZE[1]*NORMALIZE[8]*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*NORMALIZE[2]*NORMALIZE[8]*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*NORMALIZE[3]*NORMALIZE[8]*qK19;
	// num += (fabs(I[19]) >= tol);
    }
    else {
	// I[10] += scale*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*qK19;
	// num += (fabs(I[19]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK20 = ZERO;
     D128 qK22 = ZERO;
#else
     //double qK20 = 0.0;
     //double qK21 = 0.0;
     //double qK22 = 0.0;
     //double qK23 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q20 = ZERO;
	D128 q21 = ZERO;
	D128 q22 = ZERO;
	D128 q23 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    q20 = ADD(q20, MUL(MUL(x00, y01), z01));
	    q21 = ADD(q21, MUL(MUL(x10, y01), z01));
	    q22 = ADD(q22, MUL(MUL(x00, y11), z01));
	    q23 = ADD(q23, MUL(MUL(x00, y01), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    q20 = ADD1(q20, MUL1(MUL1(x00, y01), z01));
	    q21 = ADD1(q21, MUL1(MUL1(x10, y01), z01));
	    q22 = ADD1(q22, MUL1(MUL1(x00, y11), z01));
	    q23 = ADD1(q23, MUL1(MUL1(x00, y01), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C01 = LOADU(&C[k+0]);
	qK20 = ADD(qK20, MUL(C01, HADD(q20, q21)));
	D128 C11 = LOADDUP(&C[k+1]);
	qK22 = ADD(qK22, MUL(C11, HADD(q22, q23)));

#else // SSE
	    
	// function registers
	T q20 = 0.0;
	T q21 = 0.0;
	T q22 = 0.0;
	T q23 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q20 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,0,1);
	    q21 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,0,1);
	    q22 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,0,1);
	    q23 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK20 += q20*C[k+0];
	//I[20] += q20*C[k+0];
	I[20] += q20*C_[0];
	//qK21 += q21*C[k+1];
	//I[21] += q21*C[k+1];
	I[21] += q21*C_[1];
	//qK22 += q22*C[k+1];
	//I[22] += q22*C[k+1];
	I[22] += q22*C_[1];
	//qK23 += q23*C[k+1];
	//I[23] += q23*C[k+1];
	I[23] += q23*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK20 = MUL(q, qK20);
	qK22 = MUL(q, qK22);
	num += 4; //num += (fabs(I[22]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[20]) & 0xF) {
	    // 20
	    STOREU(&I[20], ADD(qK20, LOADU(&I[20])));
	    STOREU(&I[22], ADD(qK22, LOADU(&I[22])));
	}
	else {
	    STORE(&I[20], ADD(qK20, LOADU(&I[20])));
	    STORE(&I[22], ADD(qK22, LOADU(&I[22])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[20] += scale*NORMALIZE[0]*NORMALIZE[9]*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*NORMALIZE[1]*NORMALIZE[9]*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*NORMALIZE[2]*NORMALIZE[9]*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*NORMALIZE[3]*NORMALIZE[9]*qK23;
	// num += (fabs(I[23]) >= tol);
    }
    else {
	// I[20] += scale*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*qK23;
	// num += (fabs(I[23]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <sf| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::S,rysq::F> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::S,rysq::F> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 1;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x03 = LOAD(&Ix(a,0,3));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x02 = LOAD(&Ix(a,0,2));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y03 = LOAD(&Iy(a,0,3));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y02 = LOAD(&Iy(a,0,2));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z03 = LOAD(&Iz(a,0,3));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z02 = LOAD(&Iz(a,0,2));
	    q0 = ADD(q0, MUL(MUL(x03, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x00, y03), z00));
	    q2 = ADD(q2, MUL(MUL(x00, y00), z03));
	    q3 = ADD(q3, MUL(MUL(x02, y01), z00));
	    q4 = ADD(q4, MUL(MUL(x02, y00), z01));
	    q5 = ADD(q5, MUL(MUL(x01, y02), z00));
	    q6 = ADD(q6, MUL(MUL(x00, y02), z01));
	    q7 = ADD(q7, MUL(MUL(x01, y00), z02));
	    q8 = ADD(q8, MUL(MUL(x00, y01), z02));
	    q9 = ADD(q9, MUL(MUL(x01, y01), z01));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x03 = LOAD1(&Ix(a,0,3));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x02 = LOAD1(&Ix(a,0,2));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y03 = LOAD1(&Iy(a,0,3));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y02 = LOAD1(&Iy(a,0,2));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z03 = LOAD1(&Iz(a,0,3));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z02 = LOAD1(&Iz(a,0,2));
	    q0 = ADD1(q0, MUL1(MUL1(x03, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x00, y03), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x00, y00), z03));
	    q3 = ADD1(q3, MUL1(MUL1(x02, y01), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x02, y00), z01));
	    q5 = ADD1(q5, MUL1(MUL1(x01, y02), z00));
	    q6 = ADD1(q6, MUL1(MUL1(x00, y02), z01));
	    q7 = ADD1(q7, MUL1(MUL1(x01, y00), z02));
	    q8 = ADD1(q8, MUL1(MUL1(x00, y01), z02));
	    q9 = ADD1(q9, MUL1(MUL1(x01, y01), z01));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD(qK2, MUL(C00, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C00, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C00, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C00, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,0,3)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,0)*Iy(a,0,3)*Iz(a,0,0);
	    q2 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,0,3);
	    q3 += Ix(a,0,2)*Iy(a,0,1)*Iz(a,0,0);
	    q4 += Ix(a,0,2)*Iy(a,0,0)*Iz(a,0,1);
	    q5 += Ix(a,0,1)*Iy(a,0,2)*Iz(a,0,0);
	    q6 += Ix(a,0,0)*Iy(a,0,2)*Iz(a,0,1);
	    q7 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,0,2);
	    q8 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,0,2);
	    q9 += Ix(a,0,1)*Iy(a,0,1)*Iz(a,0,1);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+0];
	//I[3] += q3*C[k+0];
	I[3] += q3*C_[0];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+0];
	//I[5] += q5*C[k+0];
	I[5] += q5*C_[0];
	//qK6 += q6*C[k+0];
	//I[6] += q6*C[k+0];
	I[6] += q6*C_[0];
	//qK7 += q7*C[k+0];
	//I[7] += q7*C[k+0];
	I[7] += q7*C_[0];
	//qK8 += q8*C[k+0];
	//I[8] += q8*C[k+0];
	I[8] += q8*C_[0];
	//qK9 += q9*C[k+0];
	//I[9] += q9*C[k+0];
	I[9] += q9*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[0]*NORMALIZE[10]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[0]*NORMALIZE[11]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[0]*NORMALIZE[12]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[0]*NORMALIZE[13]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[0]*NORMALIZE[14]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[0]*NORMALIZE[15]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[0]*NORMALIZE[16]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[0]*NORMALIZE[17]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[0]*NORMALIZE[18]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[0]*NORMALIZE[19]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <pf| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::P,rysq::F> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::P,rysq::F> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 2;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x12 = LOAD(&Ix(a,1,2));
	    D128 x03 = LOAD(&Ix(a,0,3));
	    D128 x13 = LOAD(&Ix(a,1,3));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y13 = LOAD(&Iy(a,1,3));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y03 = LOAD(&Iy(a,0,3));
	    D128 z13 = LOAD(&Iz(a,1,3));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z03 = LOAD(&Iz(a,0,3));
	    q0 = ADD(q0, MUL(MUL(x13, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x03, y10), z00));
	    q2 = ADD(q2, MUL(MUL(x03, y00), z10));
	    q3 = ADD(q3, MUL(MUL(x10, y03), z00));
	    q4 = ADD(q4, MUL(MUL(x00, y13), z00));
	    q5 = ADD(q5, MUL(MUL(x00, y03), z10));
	    q6 = ADD(q6, MUL(MUL(x10, y00), z03));
	    q7 = ADD(q7, MUL(MUL(x00, y10), z03));
	    q8 = ADD(q8, MUL(MUL(x00, y00), z13));
	    q9 = ADD(q9, MUL(MUL(x12, y01), z00));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x12 = LOAD1(&Ix(a,1,2));
	    D128 x03 = LOAD1(&Ix(a,0,3));
	    D128 x13 = LOAD1(&Ix(a,1,3));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y13 = LOAD1(&Iy(a,1,3));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y03 = LOAD1(&Iy(a,0,3));
	    D128 z13 = LOAD1(&Iz(a,1,3));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z03 = LOAD1(&Iz(a,0,3));
	    q0 = ADD1(q0, MUL1(MUL1(x13, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x03, y10), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x03, y00), z10));
	    q3 = ADD1(q3, MUL1(MUL1(x10, y03), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x00, y13), z00));
	    q5 = ADD1(q5, MUL1(MUL1(x00, y03), z10));
	    q6 = ADD1(q6, MUL1(MUL1(x10, y00), z03));
	    q7 = ADD1(q7, MUL1(MUL1(x00, y10), z03));
	    q8 = ADD1(q8, MUL1(MUL1(x00, y00), z13));
	    q9 = ADD1(q9, MUL1(MUL1(x12, y01), z00));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD(qK2, MUL(C00, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C00, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C00, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C00, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,1,3)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,3)*Iy(a,1,0)*Iz(a,0,0);
	    q2 += Ix(a,0,3)*Iy(a,0,0)*Iz(a,1,0);
	    q3 += Ix(a,1,0)*Iy(a,0,3)*Iz(a,0,0);
	    q4 += Ix(a,0,0)*Iy(a,1,3)*Iz(a,0,0);
	    q5 += Ix(a,0,0)*Iy(a,0,3)*Iz(a,1,0);
	    q6 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,0,3);
	    q7 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,0,3);
	    q8 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,1,3);
	    q9 += Ix(a,1,2)*Iy(a,0,1)*Iz(a,0,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+0];
	//I[3] += q3*C[k+0];
	I[3] += q3*C_[0];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+0];
	//I[5] += q5*C[k+0];
	I[5] += q5*C_[0];
	//qK6 += q6*C[k+0];
	//I[6] += q6*C[k+0];
	I[6] += q6*C_[0];
	//qK7 += q7*C[k+0];
	//I[7] += q7*C[k+0];
	I[7] += q7*C_[0];
	//qK8 += q8*C[k+0];
	//I[8] += q8*C[k+0];
	I[8] += q8*C_[0];
	//qK9 += q9*C[k+0];
	//I[9] += q9*C[k+0];
	I[9] += q9*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[1]*NORMALIZE[10]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[2]*NORMALIZE[10]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[3]*NORMALIZE[10]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[1]*NORMALIZE[11]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[2]*NORMALIZE[11]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[3]*NORMALIZE[11]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[1]*NORMALIZE[12]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[2]*NORMALIZE[12]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[3]*NORMALIZE[12]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[1]*NORMALIZE[13]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK10 = ZERO;
     D128 qK12 = ZERO;
     D128 qK14 = ZERO;
     D128 qK16 = ZERO;
     D128 qK18 = ZERO;
#else
     //double qK10 = 0.0;
     //double qK11 = 0.0;
     //double qK12 = 0.0;
     //double qK13 = 0.0;
     //double qK14 = 0.0;
     //double qK15 = 0.0;
     //double qK16 = 0.0;
     //double qK17 = 0.0;
     //double qK18 = 0.0;
     //double qK19 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q10 = ZERO;
	D128 q11 = ZERO;
	D128 q12 = ZERO;
	D128 q13 = ZERO;
	D128 q14 = ZERO;
	D128 q15 = ZERO;
	D128 q16 = ZERO;
	D128 q17 = ZERO;
	D128 q18 = ZERO;
	D128 q19 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x12 = LOAD(&Ix(a,1,2));
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x02 = LOAD(&Ix(a,0,2));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y12 = LOAD(&Iy(a,1,2));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y02 = LOAD(&Iy(a,0,2));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    q10 = ADD(q10, MUL(MUL(x02, y11), z00));
	    q11 = ADD(q11, MUL(MUL(x02, y01), z10));
	    q12 = ADD(q12, MUL(MUL(x12, y00), z01));
	    q13 = ADD(q13, MUL(MUL(x02, y10), z01));
	    q14 = ADD(q14, MUL(MUL(x02, y00), z11));
	    q15 = ADD(q15, MUL(MUL(x11, y02), z00));
	    q16 = ADD(q16, MUL(MUL(x01, y12), z00));
	    q17 = ADD(q17, MUL(MUL(x01, y02), z10));
	    q18 = ADD(q18, MUL(MUL(x10, y02), z01));
	    q19 = ADD(q19, MUL(MUL(x00, y12), z01));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x12 = LOAD1(&Ix(a,1,2));
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x02 = LOAD1(&Ix(a,0,2));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y12 = LOAD1(&Iy(a,1,2));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y02 = LOAD1(&Iy(a,0,2));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    q10 = ADD1(q10, MUL1(MUL1(x02, y11), z00));
	    q11 = ADD1(q11, MUL1(MUL1(x02, y01), z10));
	    q12 = ADD1(q12, MUL1(MUL1(x12, y00), z01));
	    q13 = ADD1(q13, MUL1(MUL1(x02, y10), z01));
	    q14 = ADD1(q14, MUL1(MUL1(x02, y00), z11));
	    q15 = ADD1(q15, MUL1(MUL1(x11, y02), z00));
	    q16 = ADD1(q16, MUL1(MUL1(x01, y12), z00));
	    q17 = ADD1(q17, MUL1(MUL1(x01, y02), z10));
	    q18 = ADD1(q18, MUL1(MUL1(x10, y02), z01));
	    q19 = ADD1(q19, MUL1(MUL1(x00, y12), z01));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK10 = ADD(qK10, MUL(C00, HADD(q10, q11)));
	qK12 = ADD(qK12, MUL(C00, HADD(q12, q13)));
	qK14 = ADD(qK14, MUL(C00, HADD(q14, q15)));
	qK16 = ADD(qK16, MUL(C00, HADD(q16, q17)));
	qK18 = ADD(qK18, MUL(C00, HADD(q18, q19)));

#else // SSE
	    
	// function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
	T q16 = 0.0;
	T q17 = 0.0;
	T q18 = 0.0;
	T q19 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += Ix(a,0,2)*Iy(a,1,1)*Iz(a,0,0);
	    q11 += Ix(a,0,2)*Iy(a,0,1)*Iz(a,1,0);
	    q12 += Ix(a,1,2)*Iy(a,0,0)*Iz(a,0,1);
	    q13 += Ix(a,0,2)*Iy(a,1,0)*Iz(a,0,1);
	    q14 += Ix(a,0,2)*Iy(a,0,0)*Iz(a,1,1);
	    q15 += Ix(a,1,1)*Iy(a,0,2)*Iz(a,0,0);
	    q16 += Ix(a,0,1)*Iy(a,1,2)*Iz(a,0,0);
	    q17 += Ix(a,0,1)*Iy(a,0,2)*Iz(a,1,0);
	    q18 += Ix(a,1,0)*Iy(a,0,2)*Iz(a,0,1);
	    q19 += Ix(a,0,0)*Iy(a,1,2)*Iz(a,0,1);
	}
	    
	//contraction coefficients
	//qK10 += q10*C[k+0];
	//I[10] += q10*C[k+0];
	I[10] += q10*C_[0];
	//qK11 += q11*C[k+0];
	//I[11] += q11*C[k+0];
	I[11] += q11*C_[0];
	//qK12 += q12*C[k+0];
	//I[12] += q12*C[k+0];
	I[12] += q12*C_[0];
	//qK13 += q13*C[k+0];
	//I[13] += q13*C[k+0];
	I[13] += q13*C_[0];
	//qK14 += q14*C[k+0];
	//I[14] += q14*C[k+0];
	I[14] += q14*C_[0];
	//qK15 += q15*C[k+0];
	//I[15] += q15*C[k+0];
	I[15] += q15*C_[0];
	//qK16 += q16*C[k+0];
	//I[16] += q16*C[k+0];
	I[16] += q16*C_[0];
	//qK17 += q17*C[k+0];
	//I[17] += q17*C[k+0];
	I[17] += q17*C_[0];
	//qK18 += q18*C[k+0];
	//I[18] += q18*C[k+0];
	I[18] += q18*C_[0];
	//qK19 += q19*C[k+0];
	//I[19] += q19*C[k+0];
	I[19] += q19*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK10 = MUL(q, qK10);
	qK12 = MUL(q, qK12);
	qK14 = MUL(q, qK14);
	qK16 = MUL(q, qK16);
	qK18 = MUL(q, qK18);
	num += 10; //num += (fabs(I[18]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[10]) & 0xF) {
	    // 10
	    STOREU(&I[10], ADD(qK10, LOADU(&I[10])));
	    STOREU(&I[12], ADD(qK12, LOADU(&I[12])));
	    STOREU(&I[14], ADD(qK14, LOADU(&I[14])));
	    STOREU(&I[16], ADD(qK16, LOADU(&I[16])));
	    STOREU(&I[18], ADD(qK18, LOADU(&I[18])));
	}
	else {
	    STORE(&I[10], ADD(qK10, LOADU(&I[10])));
	    STORE(&I[12], ADD(qK12, LOADU(&I[12])));
	    STORE(&I[14], ADD(qK14, LOADU(&I[14])));
	    STORE(&I[16], ADD(qK16, LOADU(&I[16])));
	    STORE(&I[18], ADD(qK18, LOADU(&I[18])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[10] += scale*NORMALIZE[2]*NORMALIZE[13]*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*NORMALIZE[3]*NORMALIZE[13]*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*NORMALIZE[1]*NORMALIZE[14]*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*NORMALIZE[2]*NORMALIZE[14]*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*NORMALIZE[3]*NORMALIZE[14]*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*NORMALIZE[1]*NORMALIZE[15]*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*NORMALIZE[2]*NORMALIZE[15]*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*NORMALIZE[3]*NORMALIZE[15]*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*NORMALIZE[1]*NORMALIZE[16]*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*NORMALIZE[2]*NORMALIZE[16]*qK19;
	// num += (fabs(I[19]) >= tol);
    }
    else {
	// I[10] += scale*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*qK19;
	// num += (fabs(I[19]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK20 = ZERO;
     D128 qK22 = ZERO;
     D128 qK24 = ZERO;
     D128 qK26 = ZERO;
     D128 qK28 = ZERO;
#else
     //double qK20 = 0.0;
     //double qK21 = 0.0;
     //double qK22 = 0.0;
     //double qK23 = 0.0;
     //double qK24 = 0.0;
     //double qK25 = 0.0;
     //double qK26 = 0.0;
     //double qK27 = 0.0;
     //double qK28 = 0.0;
     //double qK29 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q20 = ZERO;
	D128 q21 = ZERO;
	D128 q22 = ZERO;
	D128 q23 = ZERO;
	D128 q24 = ZERO;
	D128 q25 = ZERO;
	D128 q26 = ZERO;
	D128 q27 = ZERO;
	D128 q28 = ZERO;
	D128 q29 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y02 = LOAD(&Iy(a,0,2));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z12 = LOAD(&Iz(a,1,2));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    D128 z02 = LOAD(&Iz(a,0,2));
	    q20 = ADD(q20, MUL(MUL(x00, y02), z11));
	    q21 = ADD(q21, MUL(MUL(x11, y00), z02));
	    q22 = ADD(q22, MUL(MUL(x01, y10), z02));
	    q23 = ADD(q23, MUL(MUL(x01, y00), z12));
	    q24 = ADD(q24, MUL(MUL(x10, y01), z02));
	    q25 = ADD(q25, MUL(MUL(x00, y11), z02));
	    q26 = ADD(q26, MUL(MUL(x00, y01), z12));
	    q27 = ADD(q27, MUL(MUL(x11, y01), z01));
	    q28 = ADD(q28, MUL(MUL(x01, y11), z01));
	    q29 = ADD(q29, MUL(MUL(x01, y01), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y02 = LOAD1(&Iy(a,0,2));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z12 = LOAD1(&Iz(a,1,2));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    D128 z02 = LOAD1(&Iz(a,0,2));
	    q20 = ADD1(q20, MUL1(MUL1(x00, y02), z11));
	    q21 = ADD1(q21, MUL1(MUL1(x11, y00), z02));
	    q22 = ADD1(q22, MUL1(MUL1(x01, y10), z02));
	    q23 = ADD1(q23, MUL1(MUL1(x01, y00), z12));
	    q24 = ADD1(q24, MUL1(MUL1(x10, y01), z02));
	    q25 = ADD1(q25, MUL1(MUL1(x00, y11), z02));
	    q26 = ADD1(q26, MUL1(MUL1(x00, y01), z12));
	    q27 = ADD1(q27, MUL1(MUL1(x11, y01), z01));
	    q28 = ADD1(q28, MUL1(MUL1(x01, y11), z01));
	    q29 = ADD1(q29, MUL1(MUL1(x01, y01), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK20 = ADD(qK20, MUL(C00, HADD(q20, q21)));
	qK22 = ADD(qK22, MUL(C00, HADD(q22, q23)));
	qK24 = ADD(qK24, MUL(C00, HADD(q24, q25)));
	qK26 = ADD(qK26, MUL(C00, HADD(q26, q27)));
	qK28 = ADD(qK28, MUL(C00, HADD(q28, q29)));

#else // SSE
	    
	// function registers
	T q20 = 0.0;
	T q21 = 0.0;
	T q22 = 0.0;
	T q23 = 0.0;
	T q24 = 0.0;
	T q25 = 0.0;
	T q26 = 0.0;
	T q27 = 0.0;
	T q28 = 0.0;
	T q29 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q20 += Ix(a,0,0)*Iy(a,0,2)*Iz(a,1,1);
	    q21 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,0,2);
	    q22 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,0,2);
	    q23 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,1,2);
	    q24 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,0,2);
	    q25 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,0,2);
	    q26 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,1,2);
	    q27 += Ix(a,1,1)*Iy(a,0,1)*Iz(a,0,1);
	    q28 += Ix(a,0,1)*Iy(a,1,1)*Iz(a,0,1);
	    q29 += Ix(a,0,1)*Iy(a,0,1)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK20 += q20*C[k+0];
	//I[20] += q20*C[k+0];
	I[20] += q20*C_[0];
	//qK21 += q21*C[k+0];
	//I[21] += q21*C[k+0];
	I[21] += q21*C_[0];
	//qK22 += q22*C[k+0];
	//I[22] += q22*C[k+0];
	I[22] += q22*C_[0];
	//qK23 += q23*C[k+0];
	//I[23] += q23*C[k+0];
	I[23] += q23*C_[0];
	//qK24 += q24*C[k+0];
	//I[24] += q24*C[k+0];
	I[24] += q24*C_[0];
	//qK25 += q25*C[k+0];
	//I[25] += q25*C[k+0];
	I[25] += q25*C_[0];
	//qK26 += q26*C[k+0];
	//I[26] += q26*C[k+0];
	I[26] += q26*C_[0];
	//qK27 += q27*C[k+0];
	//I[27] += q27*C[k+0];
	I[27] += q27*C_[0];
	//qK28 += q28*C[k+0];
	//I[28] += q28*C[k+0];
	I[28] += q28*C_[0];
	//qK29 += q29*C[k+0];
	//I[29] += q29*C[k+0];
	I[29] += q29*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK20 = MUL(q, qK20);
	qK22 = MUL(q, qK22);
	qK24 = MUL(q, qK24);
	qK26 = MUL(q, qK26);
	qK28 = MUL(q, qK28);
	num += 10; //num += (fabs(I[28]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[20]) & 0xF) {
	    // 20
	    STOREU(&I[20], ADD(qK20, LOADU(&I[20])));
	    STOREU(&I[22], ADD(qK22, LOADU(&I[22])));
	    STOREU(&I[24], ADD(qK24, LOADU(&I[24])));
	    STOREU(&I[26], ADD(qK26, LOADU(&I[26])));
	    STOREU(&I[28], ADD(qK28, LOADU(&I[28])));
	}
	else {
	    STORE(&I[20], ADD(qK20, LOADU(&I[20])));
	    STORE(&I[22], ADD(qK22, LOADU(&I[22])));
	    STORE(&I[24], ADD(qK24, LOADU(&I[24])));
	    STORE(&I[26], ADD(qK26, LOADU(&I[26])));
	    STORE(&I[28], ADD(qK28, LOADU(&I[28])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[20] += scale*NORMALIZE[3]*NORMALIZE[16]*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*NORMALIZE[1]*NORMALIZE[17]*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*NORMALIZE[2]*NORMALIZE[17]*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*NORMALIZE[3]*NORMALIZE[17]*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*NORMALIZE[1]*NORMALIZE[18]*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*NORMALIZE[2]*NORMALIZE[18]*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*NORMALIZE[3]*NORMALIZE[18]*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*NORMALIZE[1]*NORMALIZE[19]*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*NORMALIZE[2]*NORMALIZE[19]*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*NORMALIZE[3]*NORMALIZE[19]*qK29;
	// num += (fabs(I[29]) >= tol);
    }
    else {
	// I[20] += scale*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*qK29;
	// num += (fabs(I[29]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <df| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::D,rysq::F> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::D,rysq::F> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 3;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x13 = LOAD(&Ix(a,1,3));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x23 = LOAD(&Ix(a,2,3));
	    D128 x03 = LOAD(&Ix(a,0,3));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 y13 = LOAD(&Iy(a,1,3));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y23 = LOAD(&Iy(a,2,3));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y03 = LOAD(&Iy(a,0,3));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x23, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x03, y20), z00));
	    q2 = ADD(q2, MUL(MUL(x03, y00), z20));
	    q3 = ADD(q3, MUL(MUL(x13, y10), z00));
	    q4 = ADD(q4, MUL(MUL(x13, y00), z10));
	    q5 = ADD(q5, MUL(MUL(x03, y10), z10));
	    q6 = ADD(q6, MUL(MUL(x20, y03), z00));
	    q7 = ADD(q7, MUL(MUL(x00, y23), z00));
	    q8 = ADD(q8, MUL(MUL(x00, y03), z20));
	    q9 = ADD(q9, MUL(MUL(x10, y13), z00));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x13 = LOAD1(&Ix(a,1,3));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x23 = LOAD1(&Ix(a,2,3));
	    D128 x03 = LOAD1(&Ix(a,0,3));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 y13 = LOAD1(&Iy(a,1,3));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y23 = LOAD1(&Iy(a,2,3));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y03 = LOAD1(&Iy(a,0,3));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x23, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x03, y20), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x03, y00), z20));
	    q3 = ADD1(q3, MUL1(MUL1(x13, y10), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x13, y00), z10));
	    q5 = ADD1(q5, MUL1(MUL1(x03, y10), z10));
	    q6 = ADD1(q6, MUL1(MUL1(x20, y03), z00));
	    q7 = ADD1(q7, MUL1(MUL1(x00, y23), z00));
	    q8 = ADD1(q8, MUL1(MUL1(x00, y03), z20));
	    q9 = ADD1(q9, MUL1(MUL1(x10, y13), z00));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD(qK2, MUL(C00, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C00, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C00, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C00, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,2,3)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,3)*Iy(a,2,0)*Iz(a,0,0);
	    q2 += Ix(a,0,3)*Iy(a,0,0)*Iz(a,2,0);
	    q3 += Ix(a,1,3)*Iy(a,1,0)*Iz(a,0,0);
	    q4 += Ix(a,1,3)*Iy(a,0,0)*Iz(a,1,0);
	    q5 += Ix(a,0,3)*Iy(a,1,0)*Iz(a,1,0);
	    q6 += Ix(a,2,0)*Iy(a,0,3)*Iz(a,0,0);
	    q7 += Ix(a,0,0)*Iy(a,2,3)*Iz(a,0,0);
	    q8 += Ix(a,0,0)*Iy(a,0,3)*Iz(a,2,0);
	    q9 += Ix(a,1,0)*Iy(a,1,3)*Iz(a,0,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+0];
	//I[3] += q3*C[k+0];
	I[3] += q3*C_[0];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+0];
	//I[5] += q5*C[k+0];
	I[5] += q5*C_[0];
	//qK6 += q6*C[k+0];
	//I[6] += q6*C[k+0];
	I[6] += q6*C_[0];
	//qK7 += q7*C[k+0];
	//I[7] += q7*C[k+0];
	I[7] += q7*C_[0];
	//qK8 += q8*C[k+0];
	//I[8] += q8*C[k+0];
	I[8] += q8*C_[0];
	//qK9 += q9*C[k+0];
	//I[9] += q9*C[k+0];
	I[9] += q9*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[4]*NORMALIZE[10]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[5]*NORMALIZE[10]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[6]*NORMALIZE[10]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[7]*NORMALIZE[10]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[8]*NORMALIZE[10]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[9]*NORMALIZE[10]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[4]*NORMALIZE[11]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[5]*NORMALIZE[11]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[6]*NORMALIZE[11]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[7]*NORMALIZE[11]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK10 = ZERO;
     D128 qK12 = ZERO;
     D128 qK14 = ZERO;
     D128 qK16 = ZERO;
     D128 qK18 = ZERO;
#else
     //double qK10 = 0.0;
     //double qK11 = 0.0;
     //double qK12 = 0.0;
     //double qK13 = 0.0;
     //double qK14 = 0.0;
     //double qK15 = 0.0;
     //double qK16 = 0.0;
     //double qK17 = 0.0;
     //double qK18 = 0.0;
     //double qK19 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q10 = ZERO;
	D128 q11 = ZERO;
	D128 q12 = ZERO;
	D128 q13 = ZERO;
	D128 q14 = ZERO;
	D128 q15 = ZERO;
	D128 q16 = ZERO;
	D128 q17 = ZERO;
	D128 q18 = ZERO;
	D128 q19 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x02 = LOAD(&Ix(a,0,2));
	    D128 x22 = LOAD(&Ix(a,2,2));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y13 = LOAD(&Iy(a,1,3));
	    D128 y21 = LOAD(&Iy(a,2,1));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y03 = LOAD(&Iy(a,0,3));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z13 = LOAD(&Iz(a,1,3));
	    D128 z23 = LOAD(&Iz(a,2,3));
	    D128 z03 = LOAD(&Iz(a,0,3));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q10 = ADD(q10, MUL(MUL(x10, y03), z10));
	    q11 = ADD(q11, MUL(MUL(x00, y13), z10));
	    q12 = ADD(q12, MUL(MUL(x20, y00), z03));
	    q13 = ADD(q13, MUL(MUL(x00, y20), z03));
	    q14 = ADD(q14, MUL(MUL(x00, y00), z23));
	    q15 = ADD(q15, MUL(MUL(x10, y10), z03));
	    q16 = ADD(q16, MUL(MUL(x10, y00), z13));
	    q17 = ADD(q17, MUL(MUL(x00, y10), z13));
	    q18 = ADD(q18, MUL(MUL(x22, y01), z00));
	    q19 = ADD(q19, MUL(MUL(x02, y21), z00));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x02 = LOAD1(&Ix(a,0,2));
	    D128 x22 = LOAD1(&Ix(a,2,2));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y13 = LOAD1(&Iy(a,1,3));
	    D128 y21 = LOAD1(&Iy(a,2,1));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y03 = LOAD1(&Iy(a,0,3));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z13 = LOAD1(&Iz(a,1,3));
	    D128 z23 = LOAD1(&Iz(a,2,3));
	    D128 z03 = LOAD1(&Iz(a,0,3));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q10 = ADD1(q10, MUL1(MUL1(x10, y03), z10));
	    q11 = ADD1(q11, MUL1(MUL1(x00, y13), z10));
	    q12 = ADD1(q12, MUL1(MUL1(x20, y00), z03));
	    q13 = ADD1(q13, MUL1(MUL1(x00, y20), z03));
	    q14 = ADD1(q14, MUL1(MUL1(x00, y00), z23));
	    q15 = ADD1(q15, MUL1(MUL1(x10, y10), z03));
	    q16 = ADD1(q16, MUL1(MUL1(x10, y00), z13));
	    q17 = ADD1(q17, MUL1(MUL1(x00, y10), z13));
	    q18 = ADD1(q18, MUL1(MUL1(x22, y01), z00));
	    q19 = ADD1(q19, MUL1(MUL1(x02, y21), z00));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK10 = ADD(qK10, MUL(C00, HADD(q10, q11)));
	qK12 = ADD(qK12, MUL(C00, HADD(q12, q13)));
	qK14 = ADD(qK14, MUL(C00, HADD(q14, q15)));
	qK16 = ADD(qK16, MUL(C00, HADD(q16, q17)));
	qK18 = ADD(qK18, MUL(C00, HADD(q18, q19)));

#else // SSE
	    
	// function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
	T q16 = 0.0;
	T q17 = 0.0;
	T q18 = 0.0;
	T q19 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += Ix(a,1,0)*Iy(a,0,3)*Iz(a,1,0);
	    q11 += Ix(a,0,0)*Iy(a,1,3)*Iz(a,1,0);
	    q12 += Ix(a,2,0)*Iy(a,0,0)*Iz(a,0,3);
	    q13 += Ix(a,0,0)*Iy(a,2,0)*Iz(a,0,3);
	    q14 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,2,3);
	    q15 += Ix(a,1,0)*Iy(a,1,0)*Iz(a,0,3);
	    q16 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,1,3);
	    q17 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,1,3);
	    q18 += Ix(a,2,2)*Iy(a,0,1)*Iz(a,0,0);
	    q19 += Ix(a,0,2)*Iy(a,2,1)*Iz(a,0,0);
	}
	    
	//contraction coefficients
	//qK10 += q10*C[k+0];
	//I[10] += q10*C[k+0];
	I[10] += q10*C_[0];
	//qK11 += q11*C[k+0];
	//I[11] += q11*C[k+0];
	I[11] += q11*C_[0];
	//qK12 += q12*C[k+0];
	//I[12] += q12*C[k+0];
	I[12] += q12*C_[0];
	//qK13 += q13*C[k+0];
	//I[13] += q13*C[k+0];
	I[13] += q13*C_[0];
	//qK14 += q14*C[k+0];
	//I[14] += q14*C[k+0];
	I[14] += q14*C_[0];
	//qK15 += q15*C[k+0];
	//I[15] += q15*C[k+0];
	I[15] += q15*C_[0];
	//qK16 += q16*C[k+0];
	//I[16] += q16*C[k+0];
	I[16] += q16*C_[0];
	//qK17 += q17*C[k+0];
	//I[17] += q17*C[k+0];
	I[17] += q17*C_[0];
	//qK18 += q18*C[k+0];
	//I[18] += q18*C[k+0];
	I[18] += q18*C_[0];
	//qK19 += q19*C[k+0];
	//I[19] += q19*C[k+0];
	I[19] += q19*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK10 = MUL(q, qK10);
	qK12 = MUL(q, qK12);
	qK14 = MUL(q, qK14);
	qK16 = MUL(q, qK16);
	qK18 = MUL(q, qK18);
	num += 10; //num += (fabs(I[18]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[10]) & 0xF) {
	    // 10
	    STOREU(&I[10], ADD(qK10, LOADU(&I[10])));
	    STOREU(&I[12], ADD(qK12, LOADU(&I[12])));
	    STOREU(&I[14], ADD(qK14, LOADU(&I[14])));
	    STOREU(&I[16], ADD(qK16, LOADU(&I[16])));
	    STOREU(&I[18], ADD(qK18, LOADU(&I[18])));
	}
	else {
	    STORE(&I[10], ADD(qK10, LOADU(&I[10])));
	    STORE(&I[12], ADD(qK12, LOADU(&I[12])));
	    STORE(&I[14], ADD(qK14, LOADU(&I[14])));
	    STORE(&I[16], ADD(qK16, LOADU(&I[16])));
	    STORE(&I[18], ADD(qK18, LOADU(&I[18])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[10] += scale*NORMALIZE[8]*NORMALIZE[11]*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*NORMALIZE[9]*NORMALIZE[11]*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*NORMALIZE[4]*NORMALIZE[12]*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*NORMALIZE[5]*NORMALIZE[12]*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*NORMALIZE[6]*NORMALIZE[12]*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*NORMALIZE[7]*NORMALIZE[12]*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*NORMALIZE[8]*NORMALIZE[12]*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*NORMALIZE[9]*NORMALIZE[12]*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*NORMALIZE[4]*NORMALIZE[13]*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*NORMALIZE[5]*NORMALIZE[13]*qK19;
	// num += (fabs(I[19]) >= tol);
    }
    else {
	// I[10] += scale*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*qK19;
	// num += (fabs(I[19]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK20 = ZERO;
     D128 qK22 = ZERO;
     D128 qK24 = ZERO;
     D128 qK26 = ZERO;
     D128 qK28 = ZERO;
#else
     //double qK20 = 0.0;
     //double qK21 = 0.0;
     //double qK22 = 0.0;
     //double qK23 = 0.0;
     //double qK24 = 0.0;
     //double qK25 = 0.0;
     //double qK26 = 0.0;
     //double qK27 = 0.0;
     //double qK28 = 0.0;
     //double qK29 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q20 = ZERO;
	D128 q21 = ZERO;
	D128 q22 = ZERO;
	D128 q23 = ZERO;
	D128 q24 = ZERO;
	D128 q25 = ZERO;
	D128 q26 = ZERO;
	D128 q27 = ZERO;
	D128 q28 = ZERO;
	D128 q29 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x12 = LOAD(&Ix(a,1,2));
	    D128 x02 = LOAD(&Ix(a,0,2));
	    D128 x22 = LOAD(&Ix(a,2,2));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z21 = LOAD(&Iz(a,2,1));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    q20 = ADD(q20, MUL(MUL(x02, y01), z20));
	    q21 = ADD(q21, MUL(MUL(x12, y11), z00));
	    q22 = ADD(q22, MUL(MUL(x12, y01), z10));
	    q23 = ADD(q23, MUL(MUL(x02, y11), z10));
	    q24 = ADD(q24, MUL(MUL(x22, y00), z01));
	    q25 = ADD(q25, MUL(MUL(x02, y20), z01));
	    q26 = ADD(q26, MUL(MUL(x02, y00), z21));
	    q27 = ADD(q27, MUL(MUL(x12, y10), z01));
	    q28 = ADD(q28, MUL(MUL(x12, y00), z11));
	    q29 = ADD(q29, MUL(MUL(x02, y10), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x12 = LOAD1(&Ix(a,1,2));
	    D128 x02 = LOAD1(&Ix(a,0,2));
	    D128 x22 = LOAD1(&Ix(a,2,2));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z21 = LOAD1(&Iz(a,2,1));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    q20 = ADD1(q20, MUL1(MUL1(x02, y01), z20));
	    q21 = ADD1(q21, MUL1(MUL1(x12, y11), z00));
	    q22 = ADD1(q22, MUL1(MUL1(x12, y01), z10));
	    q23 = ADD1(q23, MUL1(MUL1(x02, y11), z10));
	    q24 = ADD1(q24, MUL1(MUL1(x22, y00), z01));
	    q25 = ADD1(q25, MUL1(MUL1(x02, y20), z01));
	    q26 = ADD1(q26, MUL1(MUL1(x02, y00), z21));
	    q27 = ADD1(q27, MUL1(MUL1(x12, y10), z01));
	    q28 = ADD1(q28, MUL1(MUL1(x12, y00), z11));
	    q29 = ADD1(q29, MUL1(MUL1(x02, y10), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK20 = ADD(qK20, MUL(C00, HADD(q20, q21)));
	qK22 = ADD(qK22, MUL(C00, HADD(q22, q23)));
	qK24 = ADD(qK24, MUL(C00, HADD(q24, q25)));
	qK26 = ADD(qK26, MUL(C00, HADD(q26, q27)));
	qK28 = ADD(qK28, MUL(C00, HADD(q28, q29)));

#else // SSE
	    
	// function registers
	T q20 = 0.0;
	T q21 = 0.0;
	T q22 = 0.0;
	T q23 = 0.0;
	T q24 = 0.0;
	T q25 = 0.0;
	T q26 = 0.0;
	T q27 = 0.0;
	T q28 = 0.0;
	T q29 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q20 += Ix(a,0,2)*Iy(a,0,1)*Iz(a,2,0);
	    q21 += Ix(a,1,2)*Iy(a,1,1)*Iz(a,0,0);
	    q22 += Ix(a,1,2)*Iy(a,0,1)*Iz(a,1,0);
	    q23 += Ix(a,0,2)*Iy(a,1,1)*Iz(a,1,0);
	    q24 += Ix(a,2,2)*Iy(a,0,0)*Iz(a,0,1);
	    q25 += Ix(a,0,2)*Iy(a,2,0)*Iz(a,0,1);
	    q26 += Ix(a,0,2)*Iy(a,0,0)*Iz(a,2,1);
	    q27 += Ix(a,1,2)*Iy(a,1,0)*Iz(a,0,1);
	    q28 += Ix(a,1,2)*Iy(a,0,0)*Iz(a,1,1);
	    q29 += Ix(a,0,2)*Iy(a,1,0)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK20 += q20*C[k+0];
	//I[20] += q20*C[k+0];
	I[20] += q20*C_[0];
	//qK21 += q21*C[k+0];
	//I[21] += q21*C[k+0];
	I[21] += q21*C_[0];
	//qK22 += q22*C[k+0];
	//I[22] += q22*C[k+0];
	I[22] += q22*C_[0];
	//qK23 += q23*C[k+0];
	//I[23] += q23*C[k+0];
	I[23] += q23*C_[0];
	//qK24 += q24*C[k+0];
	//I[24] += q24*C[k+0];
	I[24] += q24*C_[0];
	//qK25 += q25*C[k+0];
	//I[25] += q25*C[k+0];
	I[25] += q25*C_[0];
	//qK26 += q26*C[k+0];
	//I[26] += q26*C[k+0];
	I[26] += q26*C_[0];
	//qK27 += q27*C[k+0];
	//I[27] += q27*C[k+0];
	I[27] += q27*C_[0];
	//qK28 += q28*C[k+0];
	//I[28] += q28*C[k+0];
	I[28] += q28*C_[0];
	//qK29 += q29*C[k+0];
	//I[29] += q29*C[k+0];
	I[29] += q29*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK20 = MUL(q, qK20);
	qK22 = MUL(q, qK22);
	qK24 = MUL(q, qK24);
	qK26 = MUL(q, qK26);
	qK28 = MUL(q, qK28);
	num += 10; //num += (fabs(I[28]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[20]) & 0xF) {
	    // 20
	    STOREU(&I[20], ADD(qK20, LOADU(&I[20])));
	    STOREU(&I[22], ADD(qK22, LOADU(&I[22])));
	    STOREU(&I[24], ADD(qK24, LOADU(&I[24])));
	    STOREU(&I[26], ADD(qK26, LOADU(&I[26])));
	    STOREU(&I[28], ADD(qK28, LOADU(&I[28])));
	}
	else {
	    STORE(&I[20], ADD(qK20, LOADU(&I[20])));
	    STORE(&I[22], ADD(qK22, LOADU(&I[22])));
	    STORE(&I[24], ADD(qK24, LOADU(&I[24])));
	    STORE(&I[26], ADD(qK26, LOADU(&I[26])));
	    STORE(&I[28], ADD(qK28, LOADU(&I[28])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[20] += scale*NORMALIZE[6]*NORMALIZE[13]*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*NORMALIZE[7]*NORMALIZE[13]*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*NORMALIZE[8]*NORMALIZE[13]*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*NORMALIZE[9]*NORMALIZE[13]*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*NORMALIZE[4]*NORMALIZE[14]*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*NORMALIZE[5]*NORMALIZE[14]*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*NORMALIZE[6]*NORMALIZE[14]*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*NORMALIZE[7]*NORMALIZE[14]*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*NORMALIZE[8]*NORMALIZE[14]*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*NORMALIZE[9]*NORMALIZE[14]*qK29;
	// num += (fabs(I[29]) >= tol);
    }
    else {
	// I[20] += scale*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*qK29;
	// num += (fabs(I[29]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK30 = ZERO;
     D128 qK32 = ZERO;
     D128 qK34 = ZERO;
     D128 qK36 = ZERO;
     D128 qK38 = ZERO;
#else
     //double qK30 = 0.0;
     //double qK31 = 0.0;
     //double qK32 = 0.0;
     //double qK33 = 0.0;
     //double qK34 = 0.0;
     //double qK35 = 0.0;
     //double qK36 = 0.0;
     //double qK37 = 0.0;
     //double qK38 = 0.0;
     //double qK39 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q30 = ZERO;
	D128 q31 = ZERO;
	D128 q32 = ZERO;
	D128 q33 = ZERO;
	D128 q34 = ZERO;
	D128 q35 = ZERO;
	D128 q36 = ZERO;
	D128 q37 = ZERO;
	D128 q38 = ZERO;
	D128 q39 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x21 = LOAD(&Ix(a,2,1));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 y12 = LOAD(&Iy(a,1,2));
	    D128 y02 = LOAD(&Iy(a,0,2));
	    D128 y22 = LOAD(&Iy(a,2,2));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z21 = LOAD(&Iz(a,2,1));
	    q30 = ADD(q30, MUL(MUL(x21, y02), z00));
	    q31 = ADD(q31, MUL(MUL(x01, y22), z00));
	    q32 = ADD(q32, MUL(MUL(x01, y02), z20));
	    q33 = ADD(q33, MUL(MUL(x11, y12), z00));
	    q34 = ADD(q34, MUL(MUL(x11, y02), z10));
	    q35 = ADD(q35, MUL(MUL(x01, y12), z10));
	    q36 = ADD(q36, MUL(MUL(x20, y02), z01));
	    q37 = ADD(q37, MUL(MUL(x00, y22), z01));
	    q38 = ADD(q38, MUL(MUL(x00, y02), z21));
	    q39 = ADD(q39, MUL(MUL(x10, y12), z01));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x21 = LOAD1(&Ix(a,2,1));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 y12 = LOAD1(&Iy(a,1,2));
	    D128 y02 = LOAD1(&Iy(a,0,2));
	    D128 y22 = LOAD1(&Iy(a,2,2));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z21 = LOAD1(&Iz(a,2,1));
	    q30 = ADD1(q30, MUL1(MUL1(x21, y02), z00));
	    q31 = ADD1(q31, MUL1(MUL1(x01, y22), z00));
	    q32 = ADD1(q32, MUL1(MUL1(x01, y02), z20));
	    q33 = ADD1(q33, MUL1(MUL1(x11, y12), z00));
	    q34 = ADD1(q34, MUL1(MUL1(x11, y02), z10));
	    q35 = ADD1(q35, MUL1(MUL1(x01, y12), z10));
	    q36 = ADD1(q36, MUL1(MUL1(x20, y02), z01));
	    q37 = ADD1(q37, MUL1(MUL1(x00, y22), z01));
	    q38 = ADD1(q38, MUL1(MUL1(x00, y02), z21));
	    q39 = ADD1(q39, MUL1(MUL1(x10, y12), z01));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK30 = ADD(qK30, MUL(C00, HADD(q30, q31)));
	qK32 = ADD(qK32, MUL(C00, HADD(q32, q33)));
	qK34 = ADD(qK34, MUL(C00, HADD(q34, q35)));
	qK36 = ADD(qK36, MUL(C00, HADD(q36, q37)));
	qK38 = ADD(qK38, MUL(C00, HADD(q38, q39)));

#else // SSE
	    
	// function registers
	T q30 = 0.0;
	T q31 = 0.0;
	T q32 = 0.0;
	T q33 = 0.0;
	T q34 = 0.0;
	T q35 = 0.0;
	T q36 = 0.0;
	T q37 = 0.0;
	T q38 = 0.0;
	T q39 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q30 += Ix(a,2,1)*Iy(a,0,2)*Iz(a,0,0);
	    q31 += Ix(a,0,1)*Iy(a,2,2)*Iz(a,0,0);
	    q32 += Ix(a,0,1)*Iy(a,0,2)*Iz(a,2,0);
	    q33 += Ix(a,1,1)*Iy(a,1,2)*Iz(a,0,0);
	    q34 += Ix(a,1,1)*Iy(a,0,2)*Iz(a,1,0);
	    q35 += Ix(a,0,1)*Iy(a,1,2)*Iz(a,1,0);
	    q36 += Ix(a,2,0)*Iy(a,0,2)*Iz(a,0,1);
	    q37 += Ix(a,0,0)*Iy(a,2,2)*Iz(a,0,1);
	    q38 += Ix(a,0,0)*Iy(a,0,2)*Iz(a,2,1);
	    q39 += Ix(a,1,0)*Iy(a,1,2)*Iz(a,0,1);
	}
	    
	//contraction coefficients
	//qK30 += q30*C[k+0];
	//I[30] += q30*C[k+0];
	I[30] += q30*C_[0];
	//qK31 += q31*C[k+0];
	//I[31] += q31*C[k+0];
	I[31] += q31*C_[0];
	//qK32 += q32*C[k+0];
	//I[32] += q32*C[k+0];
	I[32] += q32*C_[0];
	//qK33 += q33*C[k+0];
	//I[33] += q33*C[k+0];
	I[33] += q33*C_[0];
	//qK34 += q34*C[k+0];
	//I[34] += q34*C[k+0];
	I[34] += q34*C_[0];
	//qK35 += q35*C[k+0];
	//I[35] += q35*C[k+0];
	I[35] += q35*C_[0];
	//qK36 += q36*C[k+0];
	//I[36] += q36*C[k+0];
	I[36] += q36*C_[0];
	//qK37 += q37*C[k+0];
	//I[37] += q37*C[k+0];
	I[37] += q37*C_[0];
	//qK38 += q38*C[k+0];
	//I[38] += q38*C[k+0];
	I[38] += q38*C_[0];
	//qK39 += q39*C[k+0];
	//I[39] += q39*C[k+0];
	I[39] += q39*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK30 = MUL(q, qK30);
	qK32 = MUL(q, qK32);
	qK34 = MUL(q, qK34);
	qK36 = MUL(q, qK36);
	qK38 = MUL(q, qK38);
	num += 10; //num += (fabs(I[38]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[30]) & 0xF) {
	    // 30
	    STOREU(&I[30], ADD(qK30, LOADU(&I[30])));
	    STOREU(&I[32], ADD(qK32, LOADU(&I[32])));
	    STOREU(&I[34], ADD(qK34, LOADU(&I[34])));
	    STOREU(&I[36], ADD(qK36, LOADU(&I[36])));
	    STOREU(&I[38], ADD(qK38, LOADU(&I[38])));
	}
	else {
	    STORE(&I[30], ADD(qK30, LOADU(&I[30])));
	    STORE(&I[32], ADD(qK32, LOADU(&I[32])));
	    STORE(&I[34], ADD(qK34, LOADU(&I[34])));
	    STORE(&I[36], ADD(qK36, LOADU(&I[36])));
	    STORE(&I[38], ADD(qK38, LOADU(&I[38])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[30] += scale*NORMALIZE[4]*NORMALIZE[15]*qK30;
	// num += (fabs(I[30]) >= tol);
	// I[31] += scale*NORMALIZE[5]*NORMALIZE[15]*qK31;
	// num += (fabs(I[31]) >= tol);
	// I[32] += scale*NORMALIZE[6]*NORMALIZE[15]*qK32;
	// num += (fabs(I[32]) >= tol);
	// I[33] += scale*NORMALIZE[7]*NORMALIZE[15]*qK33;
	// num += (fabs(I[33]) >= tol);
	// I[34] += scale*NORMALIZE[8]*NORMALIZE[15]*qK34;
	// num += (fabs(I[34]) >= tol);
	// I[35] += scale*NORMALIZE[9]*NORMALIZE[15]*qK35;
	// num += (fabs(I[35]) >= tol);
	// I[36] += scale*NORMALIZE[4]*NORMALIZE[16]*qK36;
	// num += (fabs(I[36]) >= tol);
	// I[37] += scale*NORMALIZE[5]*NORMALIZE[16]*qK37;
	// num += (fabs(I[37]) >= tol);
	// I[38] += scale*NORMALIZE[6]*NORMALIZE[16]*qK38;
	// num += (fabs(I[38]) >= tol);
	// I[39] += scale*NORMALIZE[7]*NORMALIZE[16]*qK39;
	// num += (fabs(I[39]) >= tol);
    }
    else {
	// I[30] += scale*qK30;
	// num += (fabs(I[30]) >= tol);
	// I[31] += scale*qK31;
	// num += (fabs(I[31]) >= tol);
	// I[32] += scale*qK32;
	// num += (fabs(I[32]) >= tol);
	// I[33] += scale*qK33;
	// num += (fabs(I[33]) >= tol);
	// I[34] += scale*qK34;
	// num += (fabs(I[34]) >= tol);
	// I[35] += scale*qK35;
	// num += (fabs(I[35]) >= tol);
	// I[36] += scale*qK36;
	// num += (fabs(I[36]) >= tol);
	// I[37] += scale*qK37;
	// num += (fabs(I[37]) >= tol);
	// I[38] += scale*qK38;
	// num += (fabs(I[38]) >= tol);
	// I[39] += scale*qK39;
	// num += (fabs(I[39]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK40 = ZERO;
     D128 qK42 = ZERO;
     D128 qK44 = ZERO;
     D128 qK46 = ZERO;
     D128 qK48 = ZERO;
#else
     //double qK40 = 0.0;
     //double qK41 = 0.0;
     //double qK42 = 0.0;
     //double qK43 = 0.0;
     //double qK44 = 0.0;
     //double qK45 = 0.0;
     //double qK46 = 0.0;
     //double qK47 = 0.0;
     //double qK48 = 0.0;
     //double qK49 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q40 = ZERO;
	D128 q41 = ZERO;
	D128 q42 = ZERO;
	D128 q43 = ZERO;
	D128 q44 = ZERO;
	D128 q45 = ZERO;
	D128 q46 = ZERO;
	D128 q47 = ZERO;
	D128 q48 = ZERO;
	D128 q49 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x21 = LOAD(&Ix(a,2,1));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 y12 = LOAD(&Iy(a,1,2));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y21 = LOAD(&Iy(a,2,1));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y02 = LOAD(&Iy(a,0,2));
	    D128 z12 = LOAD(&Iz(a,1,2));
	    D128 z22 = LOAD(&Iz(a,2,2));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    D128 z02 = LOAD(&Iz(a,0,2));
	    q40 = ADD(q40, MUL(MUL(x10, y02), z11));
	    q41 = ADD(q41, MUL(MUL(x00, y12), z11));
	    q42 = ADD(q42, MUL(MUL(x21, y00), z02));
	    q43 = ADD(q43, MUL(MUL(x01, y20), z02));
	    q44 = ADD(q44, MUL(MUL(x01, y00), z22));
	    q45 = ADD(q45, MUL(MUL(x11, y10), z02));
	    q46 = ADD(q46, MUL(MUL(x11, y00), z12));
	    q47 = ADD(q47, MUL(MUL(x01, y10), z12));
	    q48 = ADD(q48, MUL(MUL(x20, y01), z02));
	    q49 = ADD(q49, MUL(MUL(x00, y21), z02));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x21 = LOAD1(&Ix(a,2,1));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 y12 = LOAD1(&Iy(a,1,2));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y21 = LOAD1(&Iy(a,2,1));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y02 = LOAD1(&Iy(a,0,2));
	    D128 z12 = LOAD1(&Iz(a,1,2));
	    D128 z22 = LOAD1(&Iz(a,2,2));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    D128 z02 = LOAD1(&Iz(a,0,2));
	    q40 = ADD1(q40, MUL1(MUL1(x10, y02), z11));
	    q41 = ADD1(q41, MUL1(MUL1(x00, y12), z11));
	    q42 = ADD1(q42, MUL1(MUL1(x21, y00), z02));
	    q43 = ADD1(q43, MUL1(MUL1(x01, y20), z02));
	    q44 = ADD1(q44, MUL1(MUL1(x01, y00), z22));
	    q45 = ADD1(q45, MUL1(MUL1(x11, y10), z02));
	    q46 = ADD1(q46, MUL1(MUL1(x11, y00), z12));
	    q47 = ADD1(q47, MUL1(MUL1(x01, y10), z12));
	    q48 = ADD1(q48, MUL1(MUL1(x20, y01), z02));
	    q49 = ADD1(q49, MUL1(MUL1(x00, y21), z02));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK40 = ADD(qK40, MUL(C00, HADD(q40, q41)));
	qK42 = ADD(qK42, MUL(C00, HADD(q42, q43)));
	qK44 = ADD(qK44, MUL(C00, HADD(q44, q45)));
	qK46 = ADD(qK46, MUL(C00, HADD(q46, q47)));
	qK48 = ADD(qK48, MUL(C00, HADD(q48, q49)));

#else // SSE
	    
	// function registers
	T q40 = 0.0;
	T q41 = 0.0;
	T q42 = 0.0;
	T q43 = 0.0;
	T q44 = 0.0;
	T q45 = 0.0;
	T q46 = 0.0;
	T q47 = 0.0;
	T q48 = 0.0;
	T q49 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q40 += Ix(a,1,0)*Iy(a,0,2)*Iz(a,1,1);
	    q41 += Ix(a,0,0)*Iy(a,1,2)*Iz(a,1,1);
	    q42 += Ix(a,2,1)*Iy(a,0,0)*Iz(a,0,2);
	    q43 += Ix(a,0,1)*Iy(a,2,0)*Iz(a,0,2);
	    q44 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,2,2);
	    q45 += Ix(a,1,1)*Iy(a,1,0)*Iz(a,0,2);
	    q46 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,1,2);
	    q47 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,1,2);
	    q48 += Ix(a,2,0)*Iy(a,0,1)*Iz(a,0,2);
	    q49 += Ix(a,0,0)*Iy(a,2,1)*Iz(a,0,2);
	}
	    
	//contraction coefficients
	//qK40 += q40*C[k+0];
	//I[40] += q40*C[k+0];
	I[40] += q40*C_[0];
	//qK41 += q41*C[k+0];
	//I[41] += q41*C[k+0];
	I[41] += q41*C_[0];
	//qK42 += q42*C[k+0];
	//I[42] += q42*C[k+0];
	I[42] += q42*C_[0];
	//qK43 += q43*C[k+0];
	//I[43] += q43*C[k+0];
	I[43] += q43*C_[0];
	//qK44 += q44*C[k+0];
	//I[44] += q44*C[k+0];
	I[44] += q44*C_[0];
	//qK45 += q45*C[k+0];
	//I[45] += q45*C[k+0];
	I[45] += q45*C_[0];
	//qK46 += q46*C[k+0];
	//I[46] += q46*C[k+0];
	I[46] += q46*C_[0];
	//qK47 += q47*C[k+0];
	//I[47] += q47*C[k+0];
	I[47] += q47*C_[0];
	//qK48 += q48*C[k+0];
	//I[48] += q48*C[k+0];
	I[48] += q48*C_[0];
	//qK49 += q49*C[k+0];
	//I[49] += q49*C[k+0];
	I[49] += q49*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK40 = MUL(q, qK40);
	qK42 = MUL(q, qK42);
	qK44 = MUL(q, qK44);
	qK46 = MUL(q, qK46);
	qK48 = MUL(q, qK48);
	num += 10; //num += (fabs(I[48]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[40]) & 0xF) {
	    // 40
	    STOREU(&I[40], ADD(qK40, LOADU(&I[40])));
	    STOREU(&I[42], ADD(qK42, LOADU(&I[42])));
	    STOREU(&I[44], ADD(qK44, LOADU(&I[44])));
	    STOREU(&I[46], ADD(qK46, LOADU(&I[46])));
	    STOREU(&I[48], ADD(qK48, LOADU(&I[48])));
	}
	else {
	    STORE(&I[40], ADD(qK40, LOADU(&I[40])));
	    STORE(&I[42], ADD(qK42, LOADU(&I[42])));
	    STORE(&I[44], ADD(qK44, LOADU(&I[44])));
	    STORE(&I[46], ADD(qK46, LOADU(&I[46])));
	    STORE(&I[48], ADD(qK48, LOADU(&I[48])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[40] += scale*NORMALIZE[8]*NORMALIZE[16]*qK40;
	// num += (fabs(I[40]) >= tol);
	// I[41] += scale*NORMALIZE[9]*NORMALIZE[16]*qK41;
	// num += (fabs(I[41]) >= tol);
	// I[42] += scale*NORMALIZE[4]*NORMALIZE[17]*qK42;
	// num += (fabs(I[42]) >= tol);
	// I[43] += scale*NORMALIZE[5]*NORMALIZE[17]*qK43;
	// num += (fabs(I[43]) >= tol);
	// I[44] += scale*NORMALIZE[6]*NORMALIZE[17]*qK44;
	// num += (fabs(I[44]) >= tol);
	// I[45] += scale*NORMALIZE[7]*NORMALIZE[17]*qK45;
	// num += (fabs(I[45]) >= tol);
	// I[46] += scale*NORMALIZE[8]*NORMALIZE[17]*qK46;
	// num += (fabs(I[46]) >= tol);
	// I[47] += scale*NORMALIZE[9]*NORMALIZE[17]*qK47;
	// num += (fabs(I[47]) >= tol);
	// I[48] += scale*NORMALIZE[4]*NORMALIZE[18]*qK48;
	// num += (fabs(I[48]) >= tol);
	// I[49] += scale*NORMALIZE[5]*NORMALIZE[18]*qK49;
	// num += (fabs(I[49]) >= tol);
    }
    else {
	// I[40] += scale*qK40;
	// num += (fabs(I[40]) >= tol);
	// I[41] += scale*qK41;
	// num += (fabs(I[41]) >= tol);
	// I[42] += scale*qK42;
	// num += (fabs(I[42]) >= tol);
	// I[43] += scale*qK43;
	// num += (fabs(I[43]) >= tol);
	// I[44] += scale*qK44;
	// num += (fabs(I[44]) >= tol);
	// I[45] += scale*qK45;
	// num += (fabs(I[45]) >= tol);
	// I[46] += scale*qK46;
	// num += (fabs(I[46]) >= tol);
	// I[47] += scale*qK47;
	// num += (fabs(I[47]) >= tol);
	// I[48] += scale*qK48;
	// num += (fabs(I[48]) >= tol);
	// I[49] += scale*qK49;
	// num += (fabs(I[49]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK50 = ZERO;
     D128 qK52 = ZERO;
     D128 qK54 = ZERO;
     D128 qK56 = ZERO;
     D128 qK58 = ZERO;
#else
     //double qK50 = 0.0;
     //double qK51 = 0.0;
     //double qK52 = 0.0;
     //double qK53 = 0.0;
     //double qK54 = 0.0;
     //double qK55 = 0.0;
     //double qK56 = 0.0;
     //double qK57 = 0.0;
     //double qK58 = 0.0;
     //double qK59 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q50 = ZERO;
	D128 q51 = ZERO;
	D128 q52 = ZERO;
	D128 q53 = ZERO;
	D128 q54 = ZERO;
	D128 q55 = ZERO;
	D128 q56 = ZERO;
	D128 q57 = ZERO;
	D128 q58 = ZERO;
	D128 q59 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 x21 = LOAD(&Ix(a,2,1));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 y21 = LOAD(&Iy(a,2,1));
	    D128 z12 = LOAD(&Iz(a,1,2));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z21 = LOAD(&Iz(a,2,1));
	    D128 z02 = LOAD(&Iz(a,0,2));
	    D128 z22 = LOAD(&Iz(a,2,2));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    q50 = ADD(q50, MUL(MUL(x00, y01), z22));
	    q51 = ADD(q51, MUL(MUL(x10, y11), z02));
	    q52 = ADD(q52, MUL(MUL(x10, y01), z12));
	    q53 = ADD(q53, MUL(MUL(x00, y11), z12));
	    q54 = ADD(q54, MUL(MUL(x21, y01), z01));
	    q55 = ADD(q55, MUL(MUL(x01, y21), z01));
	    q56 = ADD(q56, MUL(MUL(x01, y01), z21));
	    q57 = ADD(q57, MUL(MUL(x11, y11), z01));
	    q58 = ADD(q58, MUL(MUL(x11, y01), z11));
	    q59 = ADD(q59, MUL(MUL(x01, y11), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 x21 = LOAD1(&Ix(a,2,1));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 y21 = LOAD1(&Iy(a,2,1));
	    D128 z12 = LOAD1(&Iz(a,1,2));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z21 = LOAD1(&Iz(a,2,1));
	    D128 z02 = LOAD1(&Iz(a,0,2));
	    D128 z22 = LOAD1(&Iz(a,2,2));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    q50 = ADD1(q50, MUL1(MUL1(x00, y01), z22));
	    q51 = ADD1(q51, MUL1(MUL1(x10, y11), z02));
	    q52 = ADD1(q52, MUL1(MUL1(x10, y01), z12));
	    q53 = ADD1(q53, MUL1(MUL1(x00, y11), z12));
	    q54 = ADD1(q54, MUL1(MUL1(x21, y01), z01));
	    q55 = ADD1(q55, MUL1(MUL1(x01, y21), z01));
	    q56 = ADD1(q56, MUL1(MUL1(x01, y01), z21));
	    q57 = ADD1(q57, MUL1(MUL1(x11, y11), z01));
	    q58 = ADD1(q58, MUL1(MUL1(x11, y01), z11));
	    q59 = ADD1(q59, MUL1(MUL1(x01, y11), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK50 = ADD(qK50, MUL(C00, HADD(q50, q51)));
	qK52 = ADD(qK52, MUL(C00, HADD(q52, q53)));
	qK54 = ADD(qK54, MUL(C00, HADD(q54, q55)));
	qK56 = ADD(qK56, MUL(C00, HADD(q56, q57)));
	qK58 = ADD(qK58, MUL(C00, HADD(q58, q59)));

#else // SSE
	    
	// function registers
	T q50 = 0.0;
	T q51 = 0.0;
	T q52 = 0.0;
	T q53 = 0.0;
	T q54 = 0.0;
	T q55 = 0.0;
	T q56 = 0.0;
	T q57 = 0.0;
	T q58 = 0.0;
	T q59 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q50 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,2,2);
	    q51 += Ix(a,1,0)*Iy(a,1,1)*Iz(a,0,2);
	    q52 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,1,2);
	    q53 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,1,2);
	    q54 += Ix(a,2,1)*Iy(a,0,1)*Iz(a,0,1);
	    q55 += Ix(a,0,1)*Iy(a,2,1)*Iz(a,0,1);
	    q56 += Ix(a,0,1)*Iy(a,0,1)*Iz(a,2,1);
	    q57 += Ix(a,1,1)*Iy(a,1,1)*Iz(a,0,1);
	    q58 += Ix(a,1,1)*Iy(a,0,1)*Iz(a,1,1);
	    q59 += Ix(a,0,1)*Iy(a,1,1)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK50 += q50*C[k+0];
	//I[50] += q50*C[k+0];
	I[50] += q50*C_[0];
	//qK51 += q51*C[k+0];
	//I[51] += q51*C[k+0];
	I[51] += q51*C_[0];
	//qK52 += q52*C[k+0];
	//I[52] += q52*C[k+0];
	I[52] += q52*C_[0];
	//qK53 += q53*C[k+0];
	//I[53] += q53*C[k+0];
	I[53] += q53*C_[0];
	//qK54 += q54*C[k+0];
	//I[54] += q54*C[k+0];
	I[54] += q54*C_[0];
	//qK55 += q55*C[k+0];
	//I[55] += q55*C[k+0];
	I[55] += q55*C_[0];
	//qK56 += q56*C[k+0];
	//I[56] += q56*C[k+0];
	I[56] += q56*C_[0];
	//qK57 += q57*C[k+0];
	//I[57] += q57*C[k+0];
	I[57] += q57*C_[0];
	//qK58 += q58*C[k+0];
	//I[58] += q58*C[k+0];
	I[58] += q58*C_[0];
	//qK59 += q59*C[k+0];
	//I[59] += q59*C[k+0];
	I[59] += q59*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK50 = MUL(q, qK50);
	qK52 = MUL(q, qK52);
	qK54 = MUL(q, qK54);
	qK56 = MUL(q, qK56);
	qK58 = MUL(q, qK58);
	num += 10; //num += (fabs(I[58]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[50]) & 0xF) {
	    // 50
	    STOREU(&I[50], ADD(qK50, LOADU(&I[50])));
	    STOREU(&I[52], ADD(qK52, LOADU(&I[52])));
	    STOREU(&I[54], ADD(qK54, LOADU(&I[54])));
	    STOREU(&I[56], ADD(qK56, LOADU(&I[56])));
	    STOREU(&I[58], ADD(qK58, LOADU(&I[58])));
	}
	else {
	    STORE(&I[50], ADD(qK50, LOADU(&I[50])));
	    STORE(&I[52], ADD(qK52, LOADU(&I[52])));
	    STORE(&I[54], ADD(qK54, LOADU(&I[54])));
	    STORE(&I[56], ADD(qK56, LOADU(&I[56])));
	    STORE(&I[58], ADD(qK58, LOADU(&I[58])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[50] += scale*NORMALIZE[6]*NORMALIZE[18]*qK50;
	// num += (fabs(I[50]) >= tol);
	// I[51] += scale*NORMALIZE[7]*NORMALIZE[18]*qK51;
	// num += (fabs(I[51]) >= tol);
	// I[52] += scale*NORMALIZE[8]*NORMALIZE[18]*qK52;
	// num += (fabs(I[52]) >= tol);
	// I[53] += scale*NORMALIZE[9]*NORMALIZE[18]*qK53;
	// num += (fabs(I[53]) >= tol);
	// I[54] += scale*NORMALIZE[4]*NORMALIZE[19]*qK54;
	// num += (fabs(I[54]) >= tol);
	// I[55] += scale*NORMALIZE[5]*NORMALIZE[19]*qK55;
	// num += (fabs(I[55]) >= tol);
	// I[56] += scale*NORMALIZE[6]*NORMALIZE[19]*qK56;
	// num += (fabs(I[56]) >= tol);
	// I[57] += scale*NORMALIZE[7]*NORMALIZE[19]*qK57;
	// num += (fabs(I[57]) >= tol);
	// I[58] += scale*NORMALIZE[8]*NORMALIZE[19]*qK58;
	// num += (fabs(I[58]) >= tol);
	// I[59] += scale*NORMALIZE[9]*NORMALIZE[19]*qK59;
	// num += (fabs(I[59]) >= tol);
    }
    else {
	// I[50] += scale*qK50;
	// num += (fabs(I[50]) >= tol);
	// I[51] += scale*qK51;
	// num += (fabs(I[51]) >= tol);
	// I[52] += scale*qK52;
	// num += (fabs(I[52]) >= tol);
	// I[53] += scale*qK53;
	// num += (fabs(I[53]) >= tol);
	// I[54] += scale*qK54;
	// num += (fabs(I[54]) >= tol);
	// I[55] += scale*qK55;
	// num += (fabs(I[55]) >= tol);
	// I[56] += scale*qK56;
	// num += (fabs(I[56]) >= tol);
	// I[57] += scale*qK57;
	// num += (fabs(I[57]) >= tol);
	// I[58] += scale*qK58;
	// num += (fabs(I[58]) >= tol);
	// I[59] += scale*qK59;
	// num += (fabs(I[59]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <ff| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::F,rysq::F> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::F,rysq::F> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 4;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x03 = LOAD(&Ix(a,0,3));
	    D128 x13 = LOAD(&Ix(a,1,3));
	    D128 x23 = LOAD(&Ix(a,2,3));
	    D128 x33 = LOAD(&Ix(a,3,3));
	    D128 y30 = LOAD(&Iy(a,3,0));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z30 = LOAD(&Iz(a,3,0));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x33, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x03, y30), z00));
	    q2 = ADD(q2, MUL(MUL(x03, y00), z30));
	    q3 = ADD(q3, MUL(MUL(x23, y10), z00));
	    q4 = ADD(q4, MUL(MUL(x23, y00), z10));
	    q5 = ADD(q5, MUL(MUL(x13, y20), z00));
	    q6 = ADD(q6, MUL(MUL(x03, y20), z10));
	    q7 = ADD(q7, MUL(MUL(x13, y00), z20));
	    q8 = ADD(q8, MUL(MUL(x03, y10), z20));
	    q9 = ADD(q9, MUL(MUL(x13, y10), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x03 = LOAD1(&Ix(a,0,3));
	    D128 x13 = LOAD1(&Ix(a,1,3));
	    D128 x23 = LOAD1(&Ix(a,2,3));
	    D128 x33 = LOAD1(&Ix(a,3,3));
	    D128 y30 = LOAD1(&Iy(a,3,0));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z30 = LOAD1(&Iz(a,3,0));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x33, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x03, y30), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x03, y00), z30));
	    q3 = ADD1(q3, MUL1(MUL1(x23, y10), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x23, y00), z10));
	    q5 = ADD1(q5, MUL1(MUL1(x13, y20), z00));
	    q6 = ADD1(q6, MUL1(MUL1(x03, y20), z10));
	    q7 = ADD1(q7, MUL1(MUL1(x13, y00), z20));
	    q8 = ADD1(q8, MUL1(MUL1(x03, y10), z20));
	    q9 = ADD1(q9, MUL1(MUL1(x13, y10), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD(qK2, MUL(C00, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C00, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C00, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C00, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,3,3)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,3)*Iy(a,3,0)*Iz(a,0,0);
	    q2 += Ix(a,0,3)*Iy(a,0,0)*Iz(a,3,0);
	    q3 += Ix(a,2,3)*Iy(a,1,0)*Iz(a,0,0);
	    q4 += Ix(a,2,3)*Iy(a,0,0)*Iz(a,1,0);
	    q5 += Ix(a,1,3)*Iy(a,2,0)*Iz(a,0,0);
	    q6 += Ix(a,0,3)*Iy(a,2,0)*Iz(a,1,0);
	    q7 += Ix(a,1,3)*Iy(a,0,0)*Iz(a,2,0);
	    q8 += Ix(a,0,3)*Iy(a,1,0)*Iz(a,2,0);
	    q9 += Ix(a,1,3)*Iy(a,1,0)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+0];
	//I[3] += q3*C[k+0];
	I[3] += q3*C_[0];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+0];
	//I[5] += q5*C[k+0];
	I[5] += q5*C_[0];
	//qK6 += q6*C[k+0];
	//I[6] += q6*C[k+0];
	I[6] += q6*C_[0];
	//qK7 += q7*C[k+0];
	//I[7] += q7*C[k+0];
	I[7] += q7*C_[0];
	//qK8 += q8*C[k+0];
	//I[8] += q8*C[k+0];
	I[8] += q8*C_[0];
	//qK9 += q9*C[k+0];
	//I[9] += q9*C[k+0];
	I[9] += q9*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[10]*NORMALIZE[10]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[11]*NORMALIZE[10]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[12]*NORMALIZE[10]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[13]*NORMALIZE[10]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[14]*NORMALIZE[10]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[15]*NORMALIZE[10]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[16]*NORMALIZE[10]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[17]*NORMALIZE[10]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[18]*NORMALIZE[10]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[19]*NORMALIZE[10]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK10 = ZERO;
     D128 qK12 = ZERO;
     D128 qK14 = ZERO;
     D128 qK16 = ZERO;
     D128 qK18 = ZERO;
#else
     //double qK10 = 0.0;
     //double qK11 = 0.0;
     //double qK12 = 0.0;
     //double qK13 = 0.0;
     //double qK14 = 0.0;
     //double qK15 = 0.0;
     //double qK16 = 0.0;
     //double qK17 = 0.0;
     //double qK18 = 0.0;
     //double qK19 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q10 = ZERO;
	D128 q11 = ZERO;
	D128 q12 = ZERO;
	D128 q13 = ZERO;
	D128 q14 = ZERO;
	D128 q15 = ZERO;
	D128 q16 = ZERO;
	D128 q17 = ZERO;
	D128 q18 = ZERO;
	D128 q19 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x30 = LOAD(&Ix(a,3,0));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y03 = LOAD(&Iy(a,0,3));
	    D128 y13 = LOAD(&Iy(a,1,3));
	    D128 y23 = LOAD(&Iy(a,2,3));
	    D128 y33 = LOAD(&Iy(a,3,3));
	    D128 z30 = LOAD(&Iz(a,3,0));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q10 = ADD(q10, MUL(MUL(x30, y03), z00));
	    q11 = ADD(q11, MUL(MUL(x00, y33), z00));
	    q12 = ADD(q12, MUL(MUL(x00, y03), z30));
	    q13 = ADD(q13, MUL(MUL(x20, y13), z00));
	    q14 = ADD(q14, MUL(MUL(x20, y03), z10));
	    q15 = ADD(q15, MUL(MUL(x10, y23), z00));
	    q16 = ADD(q16, MUL(MUL(x00, y23), z10));
	    q17 = ADD(q17, MUL(MUL(x10, y03), z20));
	    q18 = ADD(q18, MUL(MUL(x00, y13), z20));
	    q19 = ADD(q19, MUL(MUL(x10, y13), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x30 = LOAD1(&Ix(a,3,0));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y03 = LOAD1(&Iy(a,0,3));
	    D128 y13 = LOAD1(&Iy(a,1,3));
	    D128 y23 = LOAD1(&Iy(a,2,3));
	    D128 y33 = LOAD1(&Iy(a,3,3));
	    D128 z30 = LOAD1(&Iz(a,3,0));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q10 = ADD1(q10, MUL1(MUL1(x30, y03), z00));
	    q11 = ADD1(q11, MUL1(MUL1(x00, y33), z00));
	    q12 = ADD1(q12, MUL1(MUL1(x00, y03), z30));
	    q13 = ADD1(q13, MUL1(MUL1(x20, y13), z00));
	    q14 = ADD1(q14, MUL1(MUL1(x20, y03), z10));
	    q15 = ADD1(q15, MUL1(MUL1(x10, y23), z00));
	    q16 = ADD1(q16, MUL1(MUL1(x00, y23), z10));
	    q17 = ADD1(q17, MUL1(MUL1(x10, y03), z20));
	    q18 = ADD1(q18, MUL1(MUL1(x00, y13), z20));
	    q19 = ADD1(q19, MUL1(MUL1(x10, y13), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK10 = ADD(qK10, MUL(C00, HADD(q10, q11)));
	qK12 = ADD(qK12, MUL(C00, HADD(q12, q13)));
	qK14 = ADD(qK14, MUL(C00, HADD(q14, q15)));
	qK16 = ADD(qK16, MUL(C00, HADD(q16, q17)));
	qK18 = ADD(qK18, MUL(C00, HADD(q18, q19)));

#else // SSE
	    
	// function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
	T q16 = 0.0;
	T q17 = 0.0;
	T q18 = 0.0;
	T q19 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += Ix(a,3,0)*Iy(a,0,3)*Iz(a,0,0);
	    q11 += Ix(a,0,0)*Iy(a,3,3)*Iz(a,0,0);
	    q12 += Ix(a,0,0)*Iy(a,0,3)*Iz(a,3,0);
	    q13 += Ix(a,2,0)*Iy(a,1,3)*Iz(a,0,0);
	    q14 += Ix(a,2,0)*Iy(a,0,3)*Iz(a,1,0);
	    q15 += Ix(a,1,0)*Iy(a,2,3)*Iz(a,0,0);
	    q16 += Ix(a,0,0)*Iy(a,2,3)*Iz(a,1,0);
	    q17 += Ix(a,1,0)*Iy(a,0,3)*Iz(a,2,0);
	    q18 += Ix(a,0,0)*Iy(a,1,3)*Iz(a,2,0);
	    q19 += Ix(a,1,0)*Iy(a,1,3)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK10 += q10*C[k+0];
	//I[10] += q10*C[k+0];
	I[10] += q10*C_[0];
	//qK11 += q11*C[k+0];
	//I[11] += q11*C[k+0];
	I[11] += q11*C_[0];
	//qK12 += q12*C[k+0];
	//I[12] += q12*C[k+0];
	I[12] += q12*C_[0];
	//qK13 += q13*C[k+0];
	//I[13] += q13*C[k+0];
	I[13] += q13*C_[0];
	//qK14 += q14*C[k+0];
	//I[14] += q14*C[k+0];
	I[14] += q14*C_[0];
	//qK15 += q15*C[k+0];
	//I[15] += q15*C[k+0];
	I[15] += q15*C_[0];
	//qK16 += q16*C[k+0];
	//I[16] += q16*C[k+0];
	I[16] += q16*C_[0];
	//qK17 += q17*C[k+0];
	//I[17] += q17*C[k+0];
	I[17] += q17*C_[0];
	//qK18 += q18*C[k+0];
	//I[18] += q18*C[k+0];
	I[18] += q18*C_[0];
	//qK19 += q19*C[k+0];
	//I[19] += q19*C[k+0];
	I[19] += q19*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK10 = MUL(q, qK10);
	qK12 = MUL(q, qK12);
	qK14 = MUL(q, qK14);
	qK16 = MUL(q, qK16);
	qK18 = MUL(q, qK18);
	num += 10; //num += (fabs(I[18]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[10]) & 0xF) {
	    // 10
	    STOREU(&I[10], ADD(qK10, LOADU(&I[10])));
	    STOREU(&I[12], ADD(qK12, LOADU(&I[12])));
	    STOREU(&I[14], ADD(qK14, LOADU(&I[14])));
	    STOREU(&I[16], ADD(qK16, LOADU(&I[16])));
	    STOREU(&I[18], ADD(qK18, LOADU(&I[18])));
	}
	else {
	    STORE(&I[10], ADD(qK10, LOADU(&I[10])));
	    STORE(&I[12], ADD(qK12, LOADU(&I[12])));
	    STORE(&I[14], ADD(qK14, LOADU(&I[14])));
	    STORE(&I[16], ADD(qK16, LOADU(&I[16])));
	    STORE(&I[18], ADD(qK18, LOADU(&I[18])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[10] += scale*NORMALIZE[10]*NORMALIZE[11]*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*NORMALIZE[11]*NORMALIZE[11]*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*NORMALIZE[12]*NORMALIZE[11]*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*NORMALIZE[13]*NORMALIZE[11]*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*NORMALIZE[14]*NORMALIZE[11]*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*NORMALIZE[15]*NORMALIZE[11]*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*NORMALIZE[16]*NORMALIZE[11]*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*NORMALIZE[17]*NORMALIZE[11]*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*NORMALIZE[18]*NORMALIZE[11]*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*NORMALIZE[19]*NORMALIZE[11]*qK19;
	// num += (fabs(I[19]) >= tol);
    }
    else {
	// I[10] += scale*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*qK19;
	// num += (fabs(I[19]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK20 = ZERO;
     D128 qK22 = ZERO;
     D128 qK24 = ZERO;
     D128 qK26 = ZERO;
     D128 qK28 = ZERO;
#else
     //double qK20 = 0.0;
     //double qK21 = 0.0;
     //double qK22 = 0.0;
     //double qK23 = 0.0;
     //double qK24 = 0.0;
     //double qK25 = 0.0;
     //double qK26 = 0.0;
     //double qK27 = 0.0;
     //double qK28 = 0.0;
     //double qK29 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q20 = ZERO;
	D128 q21 = ZERO;
	D128 q22 = ZERO;
	D128 q23 = ZERO;
	D128 q24 = ZERO;
	D128 q25 = ZERO;
	D128 q26 = ZERO;
	D128 q27 = ZERO;
	D128 q28 = ZERO;
	D128 q29 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x30 = LOAD(&Ix(a,3,0));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y30 = LOAD(&Iy(a,3,0));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z03 = LOAD(&Iz(a,0,3));
	    D128 z13 = LOAD(&Iz(a,1,3));
	    D128 z23 = LOAD(&Iz(a,2,3));
	    D128 z33 = LOAD(&Iz(a,3,3));
	    q20 = ADD(q20, MUL(MUL(x30, y00), z03));
	    q21 = ADD(q21, MUL(MUL(x00, y30), z03));
	    q22 = ADD(q22, MUL(MUL(x00, y00), z33));
	    q23 = ADD(q23, MUL(MUL(x20, y10), z03));
	    q24 = ADD(q24, MUL(MUL(x20, y00), z13));
	    q25 = ADD(q25, MUL(MUL(x10, y20), z03));
	    q26 = ADD(q26, MUL(MUL(x00, y20), z13));
	    q27 = ADD(q27, MUL(MUL(x10, y00), z23));
	    q28 = ADD(q28, MUL(MUL(x00, y10), z23));
	    q29 = ADD(q29, MUL(MUL(x10, y10), z13));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x30 = LOAD1(&Ix(a,3,0));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y30 = LOAD1(&Iy(a,3,0));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z03 = LOAD1(&Iz(a,0,3));
	    D128 z13 = LOAD1(&Iz(a,1,3));
	    D128 z23 = LOAD1(&Iz(a,2,3));
	    D128 z33 = LOAD1(&Iz(a,3,3));
	    q20 = ADD1(q20, MUL1(MUL1(x30, y00), z03));
	    q21 = ADD1(q21, MUL1(MUL1(x00, y30), z03));
	    q22 = ADD1(q22, MUL1(MUL1(x00, y00), z33));
	    q23 = ADD1(q23, MUL1(MUL1(x20, y10), z03));
	    q24 = ADD1(q24, MUL1(MUL1(x20, y00), z13));
	    q25 = ADD1(q25, MUL1(MUL1(x10, y20), z03));
	    q26 = ADD1(q26, MUL1(MUL1(x00, y20), z13));
	    q27 = ADD1(q27, MUL1(MUL1(x10, y00), z23));
	    q28 = ADD1(q28, MUL1(MUL1(x00, y10), z23));
	    q29 = ADD1(q29, MUL1(MUL1(x10, y10), z13));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK20 = ADD(qK20, MUL(C00, HADD(q20, q21)));
	qK22 = ADD(qK22, MUL(C00, HADD(q22, q23)));
	qK24 = ADD(qK24, MUL(C00, HADD(q24, q25)));
	qK26 = ADD(qK26, MUL(C00, HADD(q26, q27)));
	qK28 = ADD(qK28, MUL(C00, HADD(q28, q29)));

#else // SSE
	    
	// function registers
	T q20 = 0.0;
	T q21 = 0.0;
	T q22 = 0.0;
	T q23 = 0.0;
	T q24 = 0.0;
	T q25 = 0.0;
	T q26 = 0.0;
	T q27 = 0.0;
	T q28 = 0.0;
	T q29 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q20 += Ix(a,3,0)*Iy(a,0,0)*Iz(a,0,3);
	    q21 += Ix(a,0,0)*Iy(a,3,0)*Iz(a,0,3);
	    q22 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,3,3);
	    q23 += Ix(a,2,0)*Iy(a,1,0)*Iz(a,0,3);
	    q24 += Ix(a,2,0)*Iy(a,0,0)*Iz(a,1,3);
	    q25 += Ix(a,1,0)*Iy(a,2,0)*Iz(a,0,3);
	    q26 += Ix(a,0,0)*Iy(a,2,0)*Iz(a,1,3);
	    q27 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,2,3);
	    q28 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,2,3);
	    q29 += Ix(a,1,0)*Iy(a,1,0)*Iz(a,1,3);
	}
	    
	//contraction coefficients
	//qK20 += q20*C[k+0];
	//I[20] += q20*C[k+0];
	I[20] += q20*C_[0];
	//qK21 += q21*C[k+0];
	//I[21] += q21*C[k+0];
	I[21] += q21*C_[0];
	//qK22 += q22*C[k+0];
	//I[22] += q22*C[k+0];
	I[22] += q22*C_[0];
	//qK23 += q23*C[k+0];
	//I[23] += q23*C[k+0];
	I[23] += q23*C_[0];
	//qK24 += q24*C[k+0];
	//I[24] += q24*C[k+0];
	I[24] += q24*C_[0];
	//qK25 += q25*C[k+0];
	//I[25] += q25*C[k+0];
	I[25] += q25*C_[0];
	//qK26 += q26*C[k+0];
	//I[26] += q26*C[k+0];
	I[26] += q26*C_[0];
	//qK27 += q27*C[k+0];
	//I[27] += q27*C[k+0];
	I[27] += q27*C_[0];
	//qK28 += q28*C[k+0];
	//I[28] += q28*C[k+0];
	I[28] += q28*C_[0];
	//qK29 += q29*C[k+0];
	//I[29] += q29*C[k+0];
	I[29] += q29*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK20 = MUL(q, qK20);
	qK22 = MUL(q, qK22);
	qK24 = MUL(q, qK24);
	qK26 = MUL(q, qK26);
	qK28 = MUL(q, qK28);
	num += 10; //num += (fabs(I[28]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[20]) & 0xF) {
	    // 20
	    STOREU(&I[20], ADD(qK20, LOADU(&I[20])));
	    STOREU(&I[22], ADD(qK22, LOADU(&I[22])));
	    STOREU(&I[24], ADD(qK24, LOADU(&I[24])));
	    STOREU(&I[26], ADD(qK26, LOADU(&I[26])));
	    STOREU(&I[28], ADD(qK28, LOADU(&I[28])));
	}
	else {
	    STORE(&I[20], ADD(qK20, LOADU(&I[20])));
	    STORE(&I[22], ADD(qK22, LOADU(&I[22])));
	    STORE(&I[24], ADD(qK24, LOADU(&I[24])));
	    STORE(&I[26], ADD(qK26, LOADU(&I[26])));
	    STORE(&I[28], ADD(qK28, LOADU(&I[28])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[20] += scale*NORMALIZE[10]*NORMALIZE[12]*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*NORMALIZE[11]*NORMALIZE[12]*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*NORMALIZE[12]*NORMALIZE[12]*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*NORMALIZE[13]*NORMALIZE[12]*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*NORMALIZE[14]*NORMALIZE[12]*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*NORMALIZE[15]*NORMALIZE[12]*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*NORMALIZE[16]*NORMALIZE[12]*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*NORMALIZE[17]*NORMALIZE[12]*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*NORMALIZE[18]*NORMALIZE[12]*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*NORMALIZE[19]*NORMALIZE[12]*qK29;
	// num += (fabs(I[29]) >= tol);
    }
    else {
	// I[20] += scale*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*qK29;
	// num += (fabs(I[29]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK30 = ZERO;
     D128 qK32 = ZERO;
     D128 qK34 = ZERO;
     D128 qK36 = ZERO;
     D128 qK38 = ZERO;
#else
     //double qK30 = 0.0;
     //double qK31 = 0.0;
     //double qK32 = 0.0;
     //double qK33 = 0.0;
     //double qK34 = 0.0;
     //double qK35 = 0.0;
     //double qK36 = 0.0;
     //double qK37 = 0.0;
     //double qK38 = 0.0;
     //double qK39 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q30 = ZERO;
	D128 q31 = ZERO;
	D128 q32 = ZERO;
	D128 q33 = ZERO;
	D128 q34 = ZERO;
	D128 q35 = ZERO;
	D128 q36 = ZERO;
	D128 q37 = ZERO;
	D128 q38 = ZERO;
	D128 q39 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x12 = LOAD(&Ix(a,1,2));
	    D128 x32 = LOAD(&Ix(a,3,2));
	    D128 x02 = LOAD(&Ix(a,0,2));
	    D128 x22 = LOAD(&Ix(a,2,2));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y31 = LOAD(&Iy(a,3,1));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 y21 = LOAD(&Iy(a,2,1));
	    D128 z30 = LOAD(&Iz(a,3,0));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q30 = ADD(q30, MUL(MUL(x32, y01), z00));
	    q31 = ADD(q31, MUL(MUL(x02, y31), z00));
	    q32 = ADD(q32, MUL(MUL(x02, y01), z30));
	    q33 = ADD(q33, MUL(MUL(x22, y11), z00));
	    q34 = ADD(q34, MUL(MUL(x22, y01), z10));
	    q35 = ADD(q35, MUL(MUL(x12, y21), z00));
	    q36 = ADD(q36, MUL(MUL(x02, y21), z10));
	    q37 = ADD(q37, MUL(MUL(x12, y01), z20));
	    q38 = ADD(q38, MUL(MUL(x02, y11), z20));
	    q39 = ADD(q39, MUL(MUL(x12, y11), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x12 = LOAD1(&Ix(a,1,2));
	    D128 x32 = LOAD1(&Ix(a,3,2));
	    D128 x02 = LOAD1(&Ix(a,0,2));
	    D128 x22 = LOAD1(&Ix(a,2,2));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y31 = LOAD1(&Iy(a,3,1));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 y21 = LOAD1(&Iy(a,2,1));
	    D128 z30 = LOAD1(&Iz(a,3,0));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q30 = ADD1(q30, MUL1(MUL1(x32, y01), z00));
	    q31 = ADD1(q31, MUL1(MUL1(x02, y31), z00));
	    q32 = ADD1(q32, MUL1(MUL1(x02, y01), z30));
	    q33 = ADD1(q33, MUL1(MUL1(x22, y11), z00));
	    q34 = ADD1(q34, MUL1(MUL1(x22, y01), z10));
	    q35 = ADD1(q35, MUL1(MUL1(x12, y21), z00));
	    q36 = ADD1(q36, MUL1(MUL1(x02, y21), z10));
	    q37 = ADD1(q37, MUL1(MUL1(x12, y01), z20));
	    q38 = ADD1(q38, MUL1(MUL1(x02, y11), z20));
	    q39 = ADD1(q39, MUL1(MUL1(x12, y11), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK30 = ADD(qK30, MUL(C00, HADD(q30, q31)));
	qK32 = ADD(qK32, MUL(C00, HADD(q32, q33)));
	qK34 = ADD(qK34, MUL(C00, HADD(q34, q35)));
	qK36 = ADD(qK36, MUL(C00, HADD(q36, q37)));
	qK38 = ADD(qK38, MUL(C00, HADD(q38, q39)));

#else // SSE
	    
	// function registers
	T q30 = 0.0;
	T q31 = 0.0;
	T q32 = 0.0;
	T q33 = 0.0;
	T q34 = 0.0;
	T q35 = 0.0;
	T q36 = 0.0;
	T q37 = 0.0;
	T q38 = 0.0;
	T q39 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q30 += Ix(a,3,2)*Iy(a,0,1)*Iz(a,0,0);
	    q31 += Ix(a,0,2)*Iy(a,3,1)*Iz(a,0,0);
	    q32 += Ix(a,0,2)*Iy(a,0,1)*Iz(a,3,0);
	    q33 += Ix(a,2,2)*Iy(a,1,1)*Iz(a,0,0);
	    q34 += Ix(a,2,2)*Iy(a,0,1)*Iz(a,1,0);
	    q35 += Ix(a,1,2)*Iy(a,2,1)*Iz(a,0,0);
	    q36 += Ix(a,0,2)*Iy(a,2,1)*Iz(a,1,0);
	    q37 += Ix(a,1,2)*Iy(a,0,1)*Iz(a,2,0);
	    q38 += Ix(a,0,2)*Iy(a,1,1)*Iz(a,2,0);
	    q39 += Ix(a,1,2)*Iy(a,1,1)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK30 += q30*C[k+0];
	//I[30] += q30*C[k+0];
	I[30] += q30*C_[0];
	//qK31 += q31*C[k+0];
	//I[31] += q31*C[k+0];
	I[31] += q31*C_[0];
	//qK32 += q32*C[k+0];
	//I[32] += q32*C[k+0];
	I[32] += q32*C_[0];
	//qK33 += q33*C[k+0];
	//I[33] += q33*C[k+0];
	I[33] += q33*C_[0];
	//qK34 += q34*C[k+0];
	//I[34] += q34*C[k+0];
	I[34] += q34*C_[0];
	//qK35 += q35*C[k+0];
	//I[35] += q35*C[k+0];
	I[35] += q35*C_[0];
	//qK36 += q36*C[k+0];
	//I[36] += q36*C[k+0];
	I[36] += q36*C_[0];
	//qK37 += q37*C[k+0];
	//I[37] += q37*C[k+0];
	I[37] += q37*C_[0];
	//qK38 += q38*C[k+0];
	//I[38] += q38*C[k+0];
	I[38] += q38*C_[0];
	//qK39 += q39*C[k+0];
	//I[39] += q39*C[k+0];
	I[39] += q39*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK30 = MUL(q, qK30);
	qK32 = MUL(q, qK32);
	qK34 = MUL(q, qK34);
	qK36 = MUL(q, qK36);
	qK38 = MUL(q, qK38);
	num += 10; //num += (fabs(I[38]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[30]) & 0xF) {
	    // 30
	    STOREU(&I[30], ADD(qK30, LOADU(&I[30])));
	    STOREU(&I[32], ADD(qK32, LOADU(&I[32])));
	    STOREU(&I[34], ADD(qK34, LOADU(&I[34])));
	    STOREU(&I[36], ADD(qK36, LOADU(&I[36])));
	    STOREU(&I[38], ADD(qK38, LOADU(&I[38])));
	}
	else {
	    STORE(&I[30], ADD(qK30, LOADU(&I[30])));
	    STORE(&I[32], ADD(qK32, LOADU(&I[32])));
	    STORE(&I[34], ADD(qK34, LOADU(&I[34])));
	    STORE(&I[36], ADD(qK36, LOADU(&I[36])));
	    STORE(&I[38], ADD(qK38, LOADU(&I[38])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[30] += scale*NORMALIZE[10]*NORMALIZE[13]*qK30;
	// num += (fabs(I[30]) >= tol);
	// I[31] += scale*NORMALIZE[11]*NORMALIZE[13]*qK31;
	// num += (fabs(I[31]) >= tol);
	// I[32] += scale*NORMALIZE[12]*NORMALIZE[13]*qK32;
	// num += (fabs(I[32]) >= tol);
	// I[33] += scale*NORMALIZE[13]*NORMALIZE[13]*qK33;
	// num += (fabs(I[33]) >= tol);
	// I[34] += scale*NORMALIZE[14]*NORMALIZE[13]*qK34;
	// num += (fabs(I[34]) >= tol);
	// I[35] += scale*NORMALIZE[15]*NORMALIZE[13]*qK35;
	// num += (fabs(I[35]) >= tol);
	// I[36] += scale*NORMALIZE[16]*NORMALIZE[13]*qK36;
	// num += (fabs(I[36]) >= tol);
	// I[37] += scale*NORMALIZE[17]*NORMALIZE[13]*qK37;
	// num += (fabs(I[37]) >= tol);
	// I[38] += scale*NORMALIZE[18]*NORMALIZE[13]*qK38;
	// num += (fabs(I[38]) >= tol);
	// I[39] += scale*NORMALIZE[19]*NORMALIZE[13]*qK39;
	// num += (fabs(I[39]) >= tol);
    }
    else {
	// I[30] += scale*qK30;
	// num += (fabs(I[30]) >= tol);
	// I[31] += scale*qK31;
	// num += (fabs(I[31]) >= tol);
	// I[32] += scale*qK32;
	// num += (fabs(I[32]) >= tol);
	// I[33] += scale*qK33;
	// num += (fabs(I[33]) >= tol);
	// I[34] += scale*qK34;
	// num += (fabs(I[34]) >= tol);
	// I[35] += scale*qK35;
	// num += (fabs(I[35]) >= tol);
	// I[36] += scale*qK36;
	// num += (fabs(I[36]) >= tol);
	// I[37] += scale*qK37;
	// num += (fabs(I[37]) >= tol);
	// I[38] += scale*qK38;
	// num += (fabs(I[38]) >= tol);
	// I[39] += scale*qK39;
	// num += (fabs(I[39]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK40 = ZERO;
     D128 qK42 = ZERO;
     D128 qK44 = ZERO;
     D128 qK46 = ZERO;
     D128 qK48 = ZERO;
#else
     //double qK40 = 0.0;
     //double qK41 = 0.0;
     //double qK42 = 0.0;
     //double qK43 = 0.0;
     //double qK44 = 0.0;
     //double qK45 = 0.0;
     //double qK46 = 0.0;
     //double qK47 = 0.0;
     //double qK48 = 0.0;
     //double qK49 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q40 = ZERO;
	D128 q41 = ZERO;
	D128 q42 = ZERO;
	D128 q43 = ZERO;
	D128 q44 = ZERO;
	D128 q45 = ZERO;
	D128 q46 = ZERO;
	D128 q47 = ZERO;
	D128 q48 = ZERO;
	D128 q49 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x12 = LOAD(&Ix(a,1,2));
	    D128 x32 = LOAD(&Ix(a,3,2));
	    D128 x02 = LOAD(&Ix(a,0,2));
	    D128 x22 = LOAD(&Ix(a,2,2));
	    D128 y30 = LOAD(&Iy(a,3,0));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z31 = LOAD(&Iz(a,3,1));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    D128 z21 = LOAD(&Iz(a,2,1));
	    q40 = ADD(q40, MUL(MUL(x32, y00), z01));
	    q41 = ADD(q41, MUL(MUL(x02, y30), z01));
	    q42 = ADD(q42, MUL(MUL(x02, y00), z31));
	    q43 = ADD(q43, MUL(MUL(x22, y10), z01));
	    q44 = ADD(q44, MUL(MUL(x22, y00), z11));
	    q45 = ADD(q45, MUL(MUL(x12, y20), z01));
	    q46 = ADD(q46, MUL(MUL(x02, y20), z11));
	    q47 = ADD(q47, MUL(MUL(x12, y00), z21));
	    q48 = ADD(q48, MUL(MUL(x02, y10), z21));
	    q49 = ADD(q49, MUL(MUL(x12, y10), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x12 = LOAD1(&Ix(a,1,2));
	    D128 x32 = LOAD1(&Ix(a,3,2));
	    D128 x02 = LOAD1(&Ix(a,0,2));
	    D128 x22 = LOAD1(&Ix(a,2,2));
	    D128 y30 = LOAD1(&Iy(a,3,0));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z31 = LOAD1(&Iz(a,3,1));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    D128 z21 = LOAD1(&Iz(a,2,1));
	    q40 = ADD1(q40, MUL1(MUL1(x32, y00), z01));
	    q41 = ADD1(q41, MUL1(MUL1(x02, y30), z01));
	    q42 = ADD1(q42, MUL1(MUL1(x02, y00), z31));
	    q43 = ADD1(q43, MUL1(MUL1(x22, y10), z01));
	    q44 = ADD1(q44, MUL1(MUL1(x22, y00), z11));
	    q45 = ADD1(q45, MUL1(MUL1(x12, y20), z01));
	    q46 = ADD1(q46, MUL1(MUL1(x02, y20), z11));
	    q47 = ADD1(q47, MUL1(MUL1(x12, y00), z21));
	    q48 = ADD1(q48, MUL1(MUL1(x02, y10), z21));
	    q49 = ADD1(q49, MUL1(MUL1(x12, y10), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK40 = ADD(qK40, MUL(C00, HADD(q40, q41)));
	qK42 = ADD(qK42, MUL(C00, HADD(q42, q43)));
	qK44 = ADD(qK44, MUL(C00, HADD(q44, q45)));
	qK46 = ADD(qK46, MUL(C00, HADD(q46, q47)));
	qK48 = ADD(qK48, MUL(C00, HADD(q48, q49)));

#else // SSE
	    
	// function registers
	T q40 = 0.0;
	T q41 = 0.0;
	T q42 = 0.0;
	T q43 = 0.0;
	T q44 = 0.0;
	T q45 = 0.0;
	T q46 = 0.0;
	T q47 = 0.0;
	T q48 = 0.0;
	T q49 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q40 += Ix(a,3,2)*Iy(a,0,0)*Iz(a,0,1);
	    q41 += Ix(a,0,2)*Iy(a,3,0)*Iz(a,0,1);
	    q42 += Ix(a,0,2)*Iy(a,0,0)*Iz(a,3,1);
	    q43 += Ix(a,2,2)*Iy(a,1,0)*Iz(a,0,1);
	    q44 += Ix(a,2,2)*Iy(a,0,0)*Iz(a,1,1);
	    q45 += Ix(a,1,2)*Iy(a,2,0)*Iz(a,0,1);
	    q46 += Ix(a,0,2)*Iy(a,2,0)*Iz(a,1,1);
	    q47 += Ix(a,1,2)*Iy(a,0,0)*Iz(a,2,1);
	    q48 += Ix(a,0,2)*Iy(a,1,0)*Iz(a,2,1);
	    q49 += Ix(a,1,2)*Iy(a,1,0)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK40 += q40*C[k+0];
	//I[40] += q40*C[k+0];
	I[40] += q40*C_[0];
	//qK41 += q41*C[k+0];
	//I[41] += q41*C[k+0];
	I[41] += q41*C_[0];
	//qK42 += q42*C[k+0];
	//I[42] += q42*C[k+0];
	I[42] += q42*C_[0];
	//qK43 += q43*C[k+0];
	//I[43] += q43*C[k+0];
	I[43] += q43*C_[0];
	//qK44 += q44*C[k+0];
	//I[44] += q44*C[k+0];
	I[44] += q44*C_[0];
	//qK45 += q45*C[k+0];
	//I[45] += q45*C[k+0];
	I[45] += q45*C_[0];
	//qK46 += q46*C[k+0];
	//I[46] += q46*C[k+0];
	I[46] += q46*C_[0];
	//qK47 += q47*C[k+0];
	//I[47] += q47*C[k+0];
	I[47] += q47*C_[0];
	//qK48 += q48*C[k+0];
	//I[48] += q48*C[k+0];
	I[48] += q48*C_[0];
	//qK49 += q49*C[k+0];
	//I[49] += q49*C[k+0];
	I[49] += q49*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK40 = MUL(q, qK40);
	qK42 = MUL(q, qK42);
	qK44 = MUL(q, qK44);
	qK46 = MUL(q, qK46);
	qK48 = MUL(q, qK48);
	num += 10; //num += (fabs(I[48]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[40]) & 0xF) {
	    // 40
	    STOREU(&I[40], ADD(qK40, LOADU(&I[40])));
	    STOREU(&I[42], ADD(qK42, LOADU(&I[42])));
	    STOREU(&I[44], ADD(qK44, LOADU(&I[44])));
	    STOREU(&I[46], ADD(qK46, LOADU(&I[46])));
	    STOREU(&I[48], ADD(qK48, LOADU(&I[48])));
	}
	else {
	    STORE(&I[40], ADD(qK40, LOADU(&I[40])));
	    STORE(&I[42], ADD(qK42, LOADU(&I[42])));
	    STORE(&I[44], ADD(qK44, LOADU(&I[44])));
	    STORE(&I[46], ADD(qK46, LOADU(&I[46])));
	    STORE(&I[48], ADD(qK48, LOADU(&I[48])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[40] += scale*NORMALIZE[10]*NORMALIZE[14]*qK40;
	// num += (fabs(I[40]) >= tol);
	// I[41] += scale*NORMALIZE[11]*NORMALIZE[14]*qK41;
	// num += (fabs(I[41]) >= tol);
	// I[42] += scale*NORMALIZE[12]*NORMALIZE[14]*qK42;
	// num += (fabs(I[42]) >= tol);
	// I[43] += scale*NORMALIZE[13]*NORMALIZE[14]*qK43;
	// num += (fabs(I[43]) >= tol);
	// I[44] += scale*NORMALIZE[14]*NORMALIZE[14]*qK44;
	// num += (fabs(I[44]) >= tol);
	// I[45] += scale*NORMALIZE[15]*NORMALIZE[14]*qK45;
	// num += (fabs(I[45]) >= tol);
	// I[46] += scale*NORMALIZE[16]*NORMALIZE[14]*qK46;
	// num += (fabs(I[46]) >= tol);
	// I[47] += scale*NORMALIZE[17]*NORMALIZE[14]*qK47;
	// num += (fabs(I[47]) >= tol);
	// I[48] += scale*NORMALIZE[18]*NORMALIZE[14]*qK48;
	// num += (fabs(I[48]) >= tol);
	// I[49] += scale*NORMALIZE[19]*NORMALIZE[14]*qK49;
	// num += (fabs(I[49]) >= tol);
    }
    else {
	// I[40] += scale*qK40;
	// num += (fabs(I[40]) >= tol);
	// I[41] += scale*qK41;
	// num += (fabs(I[41]) >= tol);
	// I[42] += scale*qK42;
	// num += (fabs(I[42]) >= tol);
	// I[43] += scale*qK43;
	// num += (fabs(I[43]) >= tol);
	// I[44] += scale*qK44;
	// num += (fabs(I[44]) >= tol);
	// I[45] += scale*qK45;
	// num += (fabs(I[45]) >= tol);
	// I[46] += scale*qK46;
	// num += (fabs(I[46]) >= tol);
	// I[47] += scale*qK47;
	// num += (fabs(I[47]) >= tol);
	// I[48] += scale*qK48;
	// num += (fabs(I[48]) >= tol);
	// I[49] += scale*qK49;
	// num += (fabs(I[49]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK50 = ZERO;
     D128 qK52 = ZERO;
     D128 qK54 = ZERO;
     D128 qK56 = ZERO;
     D128 qK58 = ZERO;
#else
     //double qK50 = 0.0;
     //double qK51 = 0.0;
     //double qK52 = 0.0;
     //double qK53 = 0.0;
     //double qK54 = 0.0;
     //double qK55 = 0.0;
     //double qK56 = 0.0;
     //double qK57 = 0.0;
     //double qK58 = 0.0;
     //double qK59 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q50 = ZERO;
	D128 q51 = ZERO;
	D128 q52 = ZERO;
	D128 q53 = ZERO;
	D128 q54 = ZERO;
	D128 q55 = ZERO;
	D128 q56 = ZERO;
	D128 q57 = ZERO;
	D128 q58 = ZERO;
	D128 q59 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x31 = LOAD(&Ix(a,3,1));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 x21 = LOAD(&Ix(a,2,1));
	    D128 y12 = LOAD(&Iy(a,1,2));
	    D128 y32 = LOAD(&Iy(a,3,2));
	    D128 y02 = LOAD(&Iy(a,0,2));
	    D128 y22 = LOAD(&Iy(a,2,2));
	    D128 z30 = LOAD(&Iz(a,3,0));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q50 = ADD(q50, MUL(MUL(x31, y02), z00));
	    q51 = ADD(q51, MUL(MUL(x01, y32), z00));
	    q52 = ADD(q52, MUL(MUL(x01, y02), z30));
	    q53 = ADD(q53, MUL(MUL(x21, y12), z00));
	    q54 = ADD(q54, MUL(MUL(x21, y02), z10));
	    q55 = ADD(q55, MUL(MUL(x11, y22), z00));
	    q56 = ADD(q56, MUL(MUL(x01, y22), z10));
	    q57 = ADD(q57, MUL(MUL(x11, y02), z20));
	    q58 = ADD(q58, MUL(MUL(x01, y12), z20));
	    q59 = ADD(q59, MUL(MUL(x11, y12), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x31 = LOAD1(&Ix(a,3,1));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 x21 = LOAD1(&Ix(a,2,1));
	    D128 y12 = LOAD1(&Iy(a,1,2));
	    D128 y32 = LOAD1(&Iy(a,3,2));
	    D128 y02 = LOAD1(&Iy(a,0,2));
	    D128 y22 = LOAD1(&Iy(a,2,2));
	    D128 z30 = LOAD1(&Iz(a,3,0));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q50 = ADD1(q50, MUL1(MUL1(x31, y02), z00));
	    q51 = ADD1(q51, MUL1(MUL1(x01, y32), z00));
	    q52 = ADD1(q52, MUL1(MUL1(x01, y02), z30));
	    q53 = ADD1(q53, MUL1(MUL1(x21, y12), z00));
	    q54 = ADD1(q54, MUL1(MUL1(x21, y02), z10));
	    q55 = ADD1(q55, MUL1(MUL1(x11, y22), z00));
	    q56 = ADD1(q56, MUL1(MUL1(x01, y22), z10));
	    q57 = ADD1(q57, MUL1(MUL1(x11, y02), z20));
	    q58 = ADD1(q58, MUL1(MUL1(x01, y12), z20));
	    q59 = ADD1(q59, MUL1(MUL1(x11, y12), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK50 = ADD(qK50, MUL(C00, HADD(q50, q51)));
	qK52 = ADD(qK52, MUL(C00, HADD(q52, q53)));
	qK54 = ADD(qK54, MUL(C00, HADD(q54, q55)));
	qK56 = ADD(qK56, MUL(C00, HADD(q56, q57)));
	qK58 = ADD(qK58, MUL(C00, HADD(q58, q59)));

#else // SSE
	    
	// function registers
	T q50 = 0.0;
	T q51 = 0.0;
	T q52 = 0.0;
	T q53 = 0.0;
	T q54 = 0.0;
	T q55 = 0.0;
	T q56 = 0.0;
	T q57 = 0.0;
	T q58 = 0.0;
	T q59 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q50 += Ix(a,3,1)*Iy(a,0,2)*Iz(a,0,0);
	    q51 += Ix(a,0,1)*Iy(a,3,2)*Iz(a,0,0);
	    q52 += Ix(a,0,1)*Iy(a,0,2)*Iz(a,3,0);
	    q53 += Ix(a,2,1)*Iy(a,1,2)*Iz(a,0,0);
	    q54 += Ix(a,2,1)*Iy(a,0,2)*Iz(a,1,0);
	    q55 += Ix(a,1,1)*Iy(a,2,2)*Iz(a,0,0);
	    q56 += Ix(a,0,1)*Iy(a,2,2)*Iz(a,1,0);
	    q57 += Ix(a,1,1)*Iy(a,0,2)*Iz(a,2,0);
	    q58 += Ix(a,0,1)*Iy(a,1,2)*Iz(a,2,0);
	    q59 += Ix(a,1,1)*Iy(a,1,2)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK50 += q50*C[k+0];
	//I[50] += q50*C[k+0];
	I[50] += q50*C_[0];
	//qK51 += q51*C[k+0];
	//I[51] += q51*C[k+0];
	I[51] += q51*C_[0];
	//qK52 += q52*C[k+0];
	//I[52] += q52*C[k+0];
	I[52] += q52*C_[0];
	//qK53 += q53*C[k+0];
	//I[53] += q53*C[k+0];
	I[53] += q53*C_[0];
	//qK54 += q54*C[k+0];
	//I[54] += q54*C[k+0];
	I[54] += q54*C_[0];
	//qK55 += q55*C[k+0];
	//I[55] += q55*C[k+0];
	I[55] += q55*C_[0];
	//qK56 += q56*C[k+0];
	//I[56] += q56*C[k+0];
	I[56] += q56*C_[0];
	//qK57 += q57*C[k+0];
	//I[57] += q57*C[k+0];
	I[57] += q57*C_[0];
	//qK58 += q58*C[k+0];
	//I[58] += q58*C[k+0];
	I[58] += q58*C_[0];
	//qK59 += q59*C[k+0];
	//I[59] += q59*C[k+0];
	I[59] += q59*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK50 = MUL(q, qK50);
	qK52 = MUL(q, qK52);
	qK54 = MUL(q, qK54);
	qK56 = MUL(q, qK56);
	qK58 = MUL(q, qK58);
	num += 10; //num += (fabs(I[58]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[50]) & 0xF) {
	    // 50
	    STOREU(&I[50], ADD(qK50, LOADU(&I[50])));
	    STOREU(&I[52], ADD(qK52, LOADU(&I[52])));
	    STOREU(&I[54], ADD(qK54, LOADU(&I[54])));
	    STOREU(&I[56], ADD(qK56, LOADU(&I[56])));
	    STOREU(&I[58], ADD(qK58, LOADU(&I[58])));
	}
	else {
	    STORE(&I[50], ADD(qK50, LOADU(&I[50])));
	    STORE(&I[52], ADD(qK52, LOADU(&I[52])));
	    STORE(&I[54], ADD(qK54, LOADU(&I[54])));
	    STORE(&I[56], ADD(qK56, LOADU(&I[56])));
	    STORE(&I[58], ADD(qK58, LOADU(&I[58])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[50] += scale*NORMALIZE[10]*NORMALIZE[15]*qK50;
	// num += (fabs(I[50]) >= tol);
	// I[51] += scale*NORMALIZE[11]*NORMALIZE[15]*qK51;
	// num += (fabs(I[51]) >= tol);
	// I[52] += scale*NORMALIZE[12]*NORMALIZE[15]*qK52;
	// num += (fabs(I[52]) >= tol);
	// I[53] += scale*NORMALIZE[13]*NORMALIZE[15]*qK53;
	// num += (fabs(I[53]) >= tol);
	// I[54] += scale*NORMALIZE[14]*NORMALIZE[15]*qK54;
	// num += (fabs(I[54]) >= tol);
	// I[55] += scale*NORMALIZE[15]*NORMALIZE[15]*qK55;
	// num += (fabs(I[55]) >= tol);
	// I[56] += scale*NORMALIZE[16]*NORMALIZE[15]*qK56;
	// num += (fabs(I[56]) >= tol);
	// I[57] += scale*NORMALIZE[17]*NORMALIZE[15]*qK57;
	// num += (fabs(I[57]) >= tol);
	// I[58] += scale*NORMALIZE[18]*NORMALIZE[15]*qK58;
	// num += (fabs(I[58]) >= tol);
	// I[59] += scale*NORMALIZE[19]*NORMALIZE[15]*qK59;
	// num += (fabs(I[59]) >= tol);
    }
    else {
	// I[50] += scale*qK50;
	// num += (fabs(I[50]) >= tol);
	// I[51] += scale*qK51;
	// num += (fabs(I[51]) >= tol);
	// I[52] += scale*qK52;
	// num += (fabs(I[52]) >= tol);
	// I[53] += scale*qK53;
	// num += (fabs(I[53]) >= tol);
	// I[54] += scale*qK54;
	// num += (fabs(I[54]) >= tol);
	// I[55] += scale*qK55;
	// num += (fabs(I[55]) >= tol);
	// I[56] += scale*qK56;
	// num += (fabs(I[56]) >= tol);
	// I[57] += scale*qK57;
	// num += (fabs(I[57]) >= tol);
	// I[58] += scale*qK58;
	// num += (fabs(I[58]) >= tol);
	// I[59] += scale*qK59;
	// num += (fabs(I[59]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK60 = ZERO;
     D128 qK62 = ZERO;
     D128 qK64 = ZERO;
     D128 qK66 = ZERO;
     D128 qK68 = ZERO;
#else
     //double qK60 = 0.0;
     //double qK61 = 0.0;
     //double qK62 = 0.0;
     //double qK63 = 0.0;
     //double qK64 = 0.0;
     //double qK65 = 0.0;
     //double qK66 = 0.0;
     //double qK67 = 0.0;
     //double qK68 = 0.0;
     //double qK69 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q60 = ZERO;
	D128 q61 = ZERO;
	D128 q62 = ZERO;
	D128 q63 = ZERO;
	D128 q64 = ZERO;
	D128 q65 = ZERO;
	D128 q66 = ZERO;
	D128 q67 = ZERO;
	D128 q68 = ZERO;
	D128 q69 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x30 = LOAD(&Ix(a,3,0));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y12 = LOAD(&Iy(a,1,2));
	    D128 y32 = LOAD(&Iy(a,3,2));
	    D128 y02 = LOAD(&Iy(a,0,2));
	    D128 y22 = LOAD(&Iy(a,2,2));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z31 = LOAD(&Iz(a,3,1));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    D128 z21 = LOAD(&Iz(a,2,1));
	    q60 = ADD(q60, MUL(MUL(x30, y02), z01));
	    q61 = ADD(q61, MUL(MUL(x00, y32), z01));
	    q62 = ADD(q62, MUL(MUL(x00, y02), z31));
	    q63 = ADD(q63, MUL(MUL(x20, y12), z01));
	    q64 = ADD(q64, MUL(MUL(x20, y02), z11));
	    q65 = ADD(q65, MUL(MUL(x10, y22), z01));
	    q66 = ADD(q66, MUL(MUL(x00, y22), z11));
	    q67 = ADD(q67, MUL(MUL(x10, y02), z21));
	    q68 = ADD(q68, MUL(MUL(x00, y12), z21));
	    q69 = ADD(q69, MUL(MUL(x10, y12), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x30 = LOAD1(&Ix(a,3,0));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y12 = LOAD1(&Iy(a,1,2));
	    D128 y32 = LOAD1(&Iy(a,3,2));
	    D128 y02 = LOAD1(&Iy(a,0,2));
	    D128 y22 = LOAD1(&Iy(a,2,2));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z31 = LOAD1(&Iz(a,3,1));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    D128 z21 = LOAD1(&Iz(a,2,1));
	    q60 = ADD1(q60, MUL1(MUL1(x30, y02), z01));
	    q61 = ADD1(q61, MUL1(MUL1(x00, y32), z01));
	    q62 = ADD1(q62, MUL1(MUL1(x00, y02), z31));
	    q63 = ADD1(q63, MUL1(MUL1(x20, y12), z01));
	    q64 = ADD1(q64, MUL1(MUL1(x20, y02), z11));
	    q65 = ADD1(q65, MUL1(MUL1(x10, y22), z01));
	    q66 = ADD1(q66, MUL1(MUL1(x00, y22), z11));
	    q67 = ADD1(q67, MUL1(MUL1(x10, y02), z21));
	    q68 = ADD1(q68, MUL1(MUL1(x00, y12), z21));
	    q69 = ADD1(q69, MUL1(MUL1(x10, y12), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK60 = ADD(qK60, MUL(C00, HADD(q60, q61)));
	qK62 = ADD(qK62, MUL(C00, HADD(q62, q63)));
	qK64 = ADD(qK64, MUL(C00, HADD(q64, q65)));
	qK66 = ADD(qK66, MUL(C00, HADD(q66, q67)));
	qK68 = ADD(qK68, MUL(C00, HADD(q68, q69)));

#else // SSE
	    
	// function registers
	T q60 = 0.0;
	T q61 = 0.0;
	T q62 = 0.0;
	T q63 = 0.0;
	T q64 = 0.0;
	T q65 = 0.0;
	T q66 = 0.0;
	T q67 = 0.0;
	T q68 = 0.0;
	T q69 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q60 += Ix(a,3,0)*Iy(a,0,2)*Iz(a,0,1);
	    q61 += Ix(a,0,0)*Iy(a,3,2)*Iz(a,0,1);
	    q62 += Ix(a,0,0)*Iy(a,0,2)*Iz(a,3,1);
	    q63 += Ix(a,2,0)*Iy(a,1,2)*Iz(a,0,1);
	    q64 += Ix(a,2,0)*Iy(a,0,2)*Iz(a,1,1);
	    q65 += Ix(a,1,0)*Iy(a,2,2)*Iz(a,0,1);
	    q66 += Ix(a,0,0)*Iy(a,2,2)*Iz(a,1,1);
	    q67 += Ix(a,1,0)*Iy(a,0,2)*Iz(a,2,1);
	    q68 += Ix(a,0,0)*Iy(a,1,2)*Iz(a,2,1);
	    q69 += Ix(a,1,0)*Iy(a,1,2)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK60 += q60*C[k+0];
	//I[60] += q60*C[k+0];
	I[60] += q60*C_[0];
	//qK61 += q61*C[k+0];
	//I[61] += q61*C[k+0];
	I[61] += q61*C_[0];
	//qK62 += q62*C[k+0];
	//I[62] += q62*C[k+0];
	I[62] += q62*C_[0];
	//qK63 += q63*C[k+0];
	//I[63] += q63*C[k+0];
	I[63] += q63*C_[0];
	//qK64 += q64*C[k+0];
	//I[64] += q64*C[k+0];
	I[64] += q64*C_[0];
	//qK65 += q65*C[k+0];
	//I[65] += q65*C[k+0];
	I[65] += q65*C_[0];
	//qK66 += q66*C[k+0];
	//I[66] += q66*C[k+0];
	I[66] += q66*C_[0];
	//qK67 += q67*C[k+0];
	//I[67] += q67*C[k+0];
	I[67] += q67*C_[0];
	//qK68 += q68*C[k+0];
	//I[68] += q68*C[k+0];
	I[68] += q68*C_[0];
	//qK69 += q69*C[k+0];
	//I[69] += q69*C[k+0];
	I[69] += q69*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK60 = MUL(q, qK60);
	qK62 = MUL(q, qK62);
	qK64 = MUL(q, qK64);
	qK66 = MUL(q, qK66);
	qK68 = MUL(q, qK68);
	num += 10; //num += (fabs(I[68]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[60]) & 0xF) {
	    // 60
	    STOREU(&I[60], ADD(qK60, LOADU(&I[60])));
	    STOREU(&I[62], ADD(qK62, LOADU(&I[62])));
	    STOREU(&I[64], ADD(qK64, LOADU(&I[64])));
	    STOREU(&I[66], ADD(qK66, LOADU(&I[66])));
	    STOREU(&I[68], ADD(qK68, LOADU(&I[68])));
	}
	else {
	    STORE(&I[60], ADD(qK60, LOADU(&I[60])));
	    STORE(&I[62], ADD(qK62, LOADU(&I[62])));
	    STORE(&I[64], ADD(qK64, LOADU(&I[64])));
	    STORE(&I[66], ADD(qK66, LOADU(&I[66])));
	    STORE(&I[68], ADD(qK68, LOADU(&I[68])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[60] += scale*NORMALIZE[10]*NORMALIZE[16]*qK60;
	// num += (fabs(I[60]) >= tol);
	// I[61] += scale*NORMALIZE[11]*NORMALIZE[16]*qK61;
	// num += (fabs(I[61]) >= tol);
	// I[62] += scale*NORMALIZE[12]*NORMALIZE[16]*qK62;
	// num += (fabs(I[62]) >= tol);
	// I[63] += scale*NORMALIZE[13]*NORMALIZE[16]*qK63;
	// num += (fabs(I[63]) >= tol);
	// I[64] += scale*NORMALIZE[14]*NORMALIZE[16]*qK64;
	// num += (fabs(I[64]) >= tol);
	// I[65] += scale*NORMALIZE[15]*NORMALIZE[16]*qK65;
	// num += (fabs(I[65]) >= tol);
	// I[66] += scale*NORMALIZE[16]*NORMALIZE[16]*qK66;
	// num += (fabs(I[66]) >= tol);
	// I[67] += scale*NORMALIZE[17]*NORMALIZE[16]*qK67;
	// num += (fabs(I[67]) >= tol);
	// I[68] += scale*NORMALIZE[18]*NORMALIZE[16]*qK68;
	// num += (fabs(I[68]) >= tol);
	// I[69] += scale*NORMALIZE[19]*NORMALIZE[16]*qK69;
	// num += (fabs(I[69]) >= tol);
    }
    else {
	// I[60] += scale*qK60;
	// num += (fabs(I[60]) >= tol);
	// I[61] += scale*qK61;
	// num += (fabs(I[61]) >= tol);
	// I[62] += scale*qK62;
	// num += (fabs(I[62]) >= tol);
	// I[63] += scale*qK63;
	// num += (fabs(I[63]) >= tol);
	// I[64] += scale*qK64;
	// num += (fabs(I[64]) >= tol);
	// I[65] += scale*qK65;
	// num += (fabs(I[65]) >= tol);
	// I[66] += scale*qK66;
	// num += (fabs(I[66]) >= tol);
	// I[67] += scale*qK67;
	// num += (fabs(I[67]) >= tol);
	// I[68] += scale*qK68;
	// num += (fabs(I[68]) >= tol);
	// I[69] += scale*qK69;
	// num += (fabs(I[69]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK70 = ZERO;
     D128 qK72 = ZERO;
     D128 qK74 = ZERO;
     D128 qK76 = ZERO;
     D128 qK78 = ZERO;
#else
     //double qK70 = 0.0;
     //double qK71 = 0.0;
     //double qK72 = 0.0;
     //double qK73 = 0.0;
     //double qK74 = 0.0;
     //double qK75 = 0.0;
     //double qK76 = 0.0;
     //double qK77 = 0.0;
     //double qK78 = 0.0;
     //double qK79 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q70 = ZERO;
	D128 q71 = ZERO;
	D128 q72 = ZERO;
	D128 q73 = ZERO;
	D128 q74 = ZERO;
	D128 q75 = ZERO;
	D128 q76 = ZERO;
	D128 q77 = ZERO;
	D128 q78 = ZERO;
	D128 q79 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x31 = LOAD(&Ix(a,3,1));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 x21 = LOAD(&Ix(a,2,1));
	    D128 y30 = LOAD(&Iy(a,3,0));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z12 = LOAD(&Iz(a,1,2));
	    D128 z32 = LOAD(&Iz(a,3,2));
	    D128 z02 = LOAD(&Iz(a,0,2));
	    D128 z22 = LOAD(&Iz(a,2,2));
	    q70 = ADD(q70, MUL(MUL(x31, y00), z02));
	    q71 = ADD(q71, MUL(MUL(x01, y30), z02));
	    q72 = ADD(q72, MUL(MUL(x01, y00), z32));
	    q73 = ADD(q73, MUL(MUL(x21, y10), z02));
	    q74 = ADD(q74, MUL(MUL(x21, y00), z12));
	    q75 = ADD(q75, MUL(MUL(x11, y20), z02));
	    q76 = ADD(q76, MUL(MUL(x01, y20), z12));
	    q77 = ADD(q77, MUL(MUL(x11, y00), z22));
	    q78 = ADD(q78, MUL(MUL(x01, y10), z22));
	    q79 = ADD(q79, MUL(MUL(x11, y10), z12));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x31 = LOAD1(&Ix(a,3,1));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 x21 = LOAD1(&Ix(a,2,1));
	    D128 y30 = LOAD1(&Iy(a,3,0));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z12 = LOAD1(&Iz(a,1,2));
	    D128 z32 = LOAD1(&Iz(a,3,2));
	    D128 z02 = LOAD1(&Iz(a,0,2));
	    D128 z22 = LOAD1(&Iz(a,2,2));
	    q70 = ADD1(q70, MUL1(MUL1(x31, y00), z02));
	    q71 = ADD1(q71, MUL1(MUL1(x01, y30), z02));
	    q72 = ADD1(q72, MUL1(MUL1(x01, y00), z32));
	    q73 = ADD1(q73, MUL1(MUL1(x21, y10), z02));
	    q74 = ADD1(q74, MUL1(MUL1(x21, y00), z12));
	    q75 = ADD1(q75, MUL1(MUL1(x11, y20), z02));
	    q76 = ADD1(q76, MUL1(MUL1(x01, y20), z12));
	    q77 = ADD1(q77, MUL1(MUL1(x11, y00), z22));
	    q78 = ADD1(q78, MUL1(MUL1(x01, y10), z22));
	    q79 = ADD1(q79, MUL1(MUL1(x11, y10), z12));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK70 = ADD(qK70, MUL(C00, HADD(q70, q71)));
	qK72 = ADD(qK72, MUL(C00, HADD(q72, q73)));
	qK74 = ADD(qK74, MUL(C00, HADD(q74, q75)));
	qK76 = ADD(qK76, MUL(C00, HADD(q76, q77)));
	qK78 = ADD(qK78, MUL(C00, HADD(q78, q79)));

#else // SSE
	    
	// function registers
	T q70 = 0.0;
	T q71 = 0.0;
	T q72 = 0.0;
	T q73 = 0.0;
	T q74 = 0.0;
	T q75 = 0.0;
	T q76 = 0.0;
	T q77 = 0.0;
	T q78 = 0.0;
	T q79 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q70 += Ix(a,3,1)*Iy(a,0,0)*Iz(a,0,2);
	    q71 += Ix(a,0,1)*Iy(a,3,0)*Iz(a,0,2);
	    q72 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,3,2);
	    q73 += Ix(a,2,1)*Iy(a,1,0)*Iz(a,0,2);
	    q74 += Ix(a,2,1)*Iy(a,0,0)*Iz(a,1,2);
	    q75 += Ix(a,1,1)*Iy(a,2,0)*Iz(a,0,2);
	    q76 += Ix(a,0,1)*Iy(a,2,0)*Iz(a,1,2);
	    q77 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,2,2);
	    q78 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,2,2);
	    q79 += Ix(a,1,1)*Iy(a,1,0)*Iz(a,1,2);
	}
	    
	//contraction coefficients
	//qK70 += q70*C[k+0];
	//I[70] += q70*C[k+0];
	I[70] += q70*C_[0];
	//qK71 += q71*C[k+0];
	//I[71] += q71*C[k+0];
	I[71] += q71*C_[0];
	//qK72 += q72*C[k+0];
	//I[72] += q72*C[k+0];
	I[72] += q72*C_[0];
	//qK73 += q73*C[k+0];
	//I[73] += q73*C[k+0];
	I[73] += q73*C_[0];
	//qK74 += q74*C[k+0];
	//I[74] += q74*C[k+0];
	I[74] += q74*C_[0];
	//qK75 += q75*C[k+0];
	//I[75] += q75*C[k+0];
	I[75] += q75*C_[0];
	//qK76 += q76*C[k+0];
	//I[76] += q76*C[k+0];
	I[76] += q76*C_[0];
	//qK77 += q77*C[k+0];
	//I[77] += q77*C[k+0];
	I[77] += q77*C_[0];
	//qK78 += q78*C[k+0];
	//I[78] += q78*C[k+0];
	I[78] += q78*C_[0];
	//qK79 += q79*C[k+0];
	//I[79] += q79*C[k+0];
	I[79] += q79*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK70 = MUL(q, qK70);
	qK72 = MUL(q, qK72);
	qK74 = MUL(q, qK74);
	qK76 = MUL(q, qK76);
	qK78 = MUL(q, qK78);
	num += 10; //num += (fabs(I[78]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[70]) & 0xF) {
	    // 70
	    STOREU(&I[70], ADD(qK70, LOADU(&I[70])));
	    STOREU(&I[72], ADD(qK72, LOADU(&I[72])));
	    STOREU(&I[74], ADD(qK74, LOADU(&I[74])));
	    STOREU(&I[76], ADD(qK76, LOADU(&I[76])));
	    STOREU(&I[78], ADD(qK78, LOADU(&I[78])));
	}
	else {
	    STORE(&I[70], ADD(qK70, LOADU(&I[70])));
	    STORE(&I[72], ADD(qK72, LOADU(&I[72])));
	    STORE(&I[74], ADD(qK74, LOADU(&I[74])));
	    STORE(&I[76], ADD(qK76, LOADU(&I[76])));
	    STORE(&I[78], ADD(qK78, LOADU(&I[78])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[70] += scale*NORMALIZE[10]*NORMALIZE[17]*qK70;
	// num += (fabs(I[70]) >= tol);
	// I[71] += scale*NORMALIZE[11]*NORMALIZE[17]*qK71;
	// num += (fabs(I[71]) >= tol);
	// I[72] += scale*NORMALIZE[12]*NORMALIZE[17]*qK72;
	// num += (fabs(I[72]) >= tol);
	// I[73] += scale*NORMALIZE[13]*NORMALIZE[17]*qK73;
	// num += (fabs(I[73]) >= tol);
	// I[74] += scale*NORMALIZE[14]*NORMALIZE[17]*qK74;
	// num += (fabs(I[74]) >= tol);
	// I[75] += scale*NORMALIZE[15]*NORMALIZE[17]*qK75;
	// num += (fabs(I[75]) >= tol);
	// I[76] += scale*NORMALIZE[16]*NORMALIZE[17]*qK76;
	// num += (fabs(I[76]) >= tol);
	// I[77] += scale*NORMALIZE[17]*NORMALIZE[17]*qK77;
	// num += (fabs(I[77]) >= tol);
	// I[78] += scale*NORMALIZE[18]*NORMALIZE[17]*qK78;
	// num += (fabs(I[78]) >= tol);
	// I[79] += scale*NORMALIZE[19]*NORMALIZE[17]*qK79;
	// num += (fabs(I[79]) >= tol);
    }
    else {
	// I[70] += scale*qK70;
	// num += (fabs(I[70]) >= tol);
	// I[71] += scale*qK71;
	// num += (fabs(I[71]) >= tol);
	// I[72] += scale*qK72;
	// num += (fabs(I[72]) >= tol);
	// I[73] += scale*qK73;
	// num += (fabs(I[73]) >= tol);
	// I[74] += scale*qK74;
	// num += (fabs(I[74]) >= tol);
	// I[75] += scale*qK75;
	// num += (fabs(I[75]) >= tol);
	// I[76] += scale*qK76;
	// num += (fabs(I[76]) >= tol);
	// I[77] += scale*qK77;
	// num += (fabs(I[77]) >= tol);
	// I[78] += scale*qK78;
	// num += (fabs(I[78]) >= tol);
	// I[79] += scale*qK79;
	// num += (fabs(I[79]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK80 = ZERO;
     D128 qK82 = ZERO;
     D128 qK84 = ZERO;
     D128 qK86 = ZERO;
     D128 qK88 = ZERO;
#else
     //double qK80 = 0.0;
     //double qK81 = 0.0;
     //double qK82 = 0.0;
     //double qK83 = 0.0;
     //double qK84 = 0.0;
     //double qK85 = 0.0;
     //double qK86 = 0.0;
     //double qK87 = 0.0;
     //double qK88 = 0.0;
     //double qK89 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q80 = ZERO;
	D128 q81 = ZERO;
	D128 q82 = ZERO;
	D128 q83 = ZERO;
	D128 q84 = ZERO;
	D128 q85 = ZERO;
	D128 q86 = ZERO;
	D128 q87 = ZERO;
	D128 q88 = ZERO;
	D128 q89 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x30 = LOAD(&Ix(a,3,0));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y31 = LOAD(&Iy(a,3,1));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 y21 = LOAD(&Iy(a,2,1));
	    D128 z12 = LOAD(&Iz(a,1,2));
	    D128 z32 = LOAD(&Iz(a,3,2));
	    D128 z02 = LOAD(&Iz(a,0,2));
	    D128 z22 = LOAD(&Iz(a,2,2));
	    q80 = ADD(q80, MUL(MUL(x30, y01), z02));
	    q81 = ADD(q81, MUL(MUL(x00, y31), z02));
	    q82 = ADD(q82, MUL(MUL(x00, y01), z32));
	    q83 = ADD(q83, MUL(MUL(x20, y11), z02));
	    q84 = ADD(q84, MUL(MUL(x20, y01), z12));
	    q85 = ADD(q85, MUL(MUL(x10, y21), z02));
	    q86 = ADD(q86, MUL(MUL(x00, y21), z12));
	    q87 = ADD(q87, MUL(MUL(x10, y01), z22));
	    q88 = ADD(q88, MUL(MUL(x00, y11), z22));
	    q89 = ADD(q89, MUL(MUL(x10, y11), z12));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x30 = LOAD1(&Ix(a,3,0));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y31 = LOAD1(&Iy(a,3,1));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 y21 = LOAD1(&Iy(a,2,1));
	    D128 z12 = LOAD1(&Iz(a,1,2));
	    D128 z32 = LOAD1(&Iz(a,3,2));
	    D128 z02 = LOAD1(&Iz(a,0,2));
	    D128 z22 = LOAD1(&Iz(a,2,2));
	    q80 = ADD1(q80, MUL1(MUL1(x30, y01), z02));
	    q81 = ADD1(q81, MUL1(MUL1(x00, y31), z02));
	    q82 = ADD1(q82, MUL1(MUL1(x00, y01), z32));
	    q83 = ADD1(q83, MUL1(MUL1(x20, y11), z02));
	    q84 = ADD1(q84, MUL1(MUL1(x20, y01), z12));
	    q85 = ADD1(q85, MUL1(MUL1(x10, y21), z02));
	    q86 = ADD1(q86, MUL1(MUL1(x00, y21), z12));
	    q87 = ADD1(q87, MUL1(MUL1(x10, y01), z22));
	    q88 = ADD1(q88, MUL1(MUL1(x00, y11), z22));
	    q89 = ADD1(q89, MUL1(MUL1(x10, y11), z12));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK80 = ADD(qK80, MUL(C00, HADD(q80, q81)));
	qK82 = ADD(qK82, MUL(C00, HADD(q82, q83)));
	qK84 = ADD(qK84, MUL(C00, HADD(q84, q85)));
	qK86 = ADD(qK86, MUL(C00, HADD(q86, q87)));
	qK88 = ADD(qK88, MUL(C00, HADD(q88, q89)));

#else // SSE
	    
	// function registers
	T q80 = 0.0;
	T q81 = 0.0;
	T q82 = 0.0;
	T q83 = 0.0;
	T q84 = 0.0;
	T q85 = 0.0;
	T q86 = 0.0;
	T q87 = 0.0;
	T q88 = 0.0;
	T q89 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q80 += Ix(a,3,0)*Iy(a,0,1)*Iz(a,0,2);
	    q81 += Ix(a,0,0)*Iy(a,3,1)*Iz(a,0,2);
	    q82 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,3,2);
	    q83 += Ix(a,2,0)*Iy(a,1,1)*Iz(a,0,2);
	    q84 += Ix(a,2,0)*Iy(a,0,1)*Iz(a,1,2);
	    q85 += Ix(a,1,0)*Iy(a,2,1)*Iz(a,0,2);
	    q86 += Ix(a,0,0)*Iy(a,2,1)*Iz(a,1,2);
	    q87 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,2,2);
	    q88 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,2,2);
	    q89 += Ix(a,1,0)*Iy(a,1,1)*Iz(a,1,2);
	}
	    
	//contraction coefficients
	//qK80 += q80*C[k+0];
	//I[80] += q80*C[k+0];
	I[80] += q80*C_[0];
	//qK81 += q81*C[k+0];
	//I[81] += q81*C[k+0];
	I[81] += q81*C_[0];
	//qK82 += q82*C[k+0];
	//I[82] += q82*C[k+0];
	I[82] += q82*C_[0];
	//qK83 += q83*C[k+0];
	//I[83] += q83*C[k+0];
	I[83] += q83*C_[0];
	//qK84 += q84*C[k+0];
	//I[84] += q84*C[k+0];
	I[84] += q84*C_[0];
	//qK85 += q85*C[k+0];
	//I[85] += q85*C[k+0];
	I[85] += q85*C_[0];
	//qK86 += q86*C[k+0];
	//I[86] += q86*C[k+0];
	I[86] += q86*C_[0];
	//qK87 += q87*C[k+0];
	//I[87] += q87*C[k+0];
	I[87] += q87*C_[0];
	//qK88 += q88*C[k+0];
	//I[88] += q88*C[k+0];
	I[88] += q88*C_[0];
	//qK89 += q89*C[k+0];
	//I[89] += q89*C[k+0];
	I[89] += q89*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK80 = MUL(q, qK80);
	qK82 = MUL(q, qK82);
	qK84 = MUL(q, qK84);
	qK86 = MUL(q, qK86);
	qK88 = MUL(q, qK88);
	num += 10; //num += (fabs(I[88]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[80]) & 0xF) {
	    // 80
	    STOREU(&I[80], ADD(qK80, LOADU(&I[80])));
	    STOREU(&I[82], ADD(qK82, LOADU(&I[82])));
	    STOREU(&I[84], ADD(qK84, LOADU(&I[84])));
	    STOREU(&I[86], ADD(qK86, LOADU(&I[86])));
	    STOREU(&I[88], ADD(qK88, LOADU(&I[88])));
	}
	else {
	    STORE(&I[80], ADD(qK80, LOADU(&I[80])));
	    STORE(&I[82], ADD(qK82, LOADU(&I[82])));
	    STORE(&I[84], ADD(qK84, LOADU(&I[84])));
	    STORE(&I[86], ADD(qK86, LOADU(&I[86])));
	    STORE(&I[88], ADD(qK88, LOADU(&I[88])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[80] += scale*NORMALIZE[10]*NORMALIZE[18]*qK80;
	// num += (fabs(I[80]) >= tol);
	// I[81] += scale*NORMALIZE[11]*NORMALIZE[18]*qK81;
	// num += (fabs(I[81]) >= tol);
	// I[82] += scale*NORMALIZE[12]*NORMALIZE[18]*qK82;
	// num += (fabs(I[82]) >= tol);
	// I[83] += scale*NORMALIZE[13]*NORMALIZE[18]*qK83;
	// num += (fabs(I[83]) >= tol);
	// I[84] += scale*NORMALIZE[14]*NORMALIZE[18]*qK84;
	// num += (fabs(I[84]) >= tol);
	// I[85] += scale*NORMALIZE[15]*NORMALIZE[18]*qK85;
	// num += (fabs(I[85]) >= tol);
	// I[86] += scale*NORMALIZE[16]*NORMALIZE[18]*qK86;
	// num += (fabs(I[86]) >= tol);
	// I[87] += scale*NORMALIZE[17]*NORMALIZE[18]*qK87;
	// num += (fabs(I[87]) >= tol);
	// I[88] += scale*NORMALIZE[18]*NORMALIZE[18]*qK88;
	// num += (fabs(I[88]) >= tol);
	// I[89] += scale*NORMALIZE[19]*NORMALIZE[18]*qK89;
	// num += (fabs(I[89]) >= tol);
    }
    else {
	// I[80] += scale*qK80;
	// num += (fabs(I[80]) >= tol);
	// I[81] += scale*qK81;
	// num += (fabs(I[81]) >= tol);
	// I[82] += scale*qK82;
	// num += (fabs(I[82]) >= tol);
	// I[83] += scale*qK83;
	// num += (fabs(I[83]) >= tol);
	// I[84] += scale*qK84;
	// num += (fabs(I[84]) >= tol);
	// I[85] += scale*qK85;
	// num += (fabs(I[85]) >= tol);
	// I[86] += scale*qK86;
	// num += (fabs(I[86]) >= tol);
	// I[87] += scale*qK87;
	// num += (fabs(I[87]) >= tol);
	// I[88] += scale*qK88;
	// num += (fabs(I[88]) >= tol);
	// I[89] += scale*qK89;
	// num += (fabs(I[89]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK90 = ZERO;
     D128 qK92 = ZERO;
     D128 qK94 = ZERO;
     D128 qK96 = ZERO;
     D128 qK98 = ZERO;
#else
     //double qK90 = 0.0;
     //double qK91 = 0.0;
     //double qK92 = 0.0;
     //double qK93 = 0.0;
     //double qK94 = 0.0;
     //double qK95 = 0.0;
     //double qK96 = 0.0;
     //double qK97 = 0.0;
     //double qK98 = 0.0;
     //double qK99 = 0.0;
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q90 = ZERO;
	D128 q91 = ZERO;
	D128 q92 = ZERO;
	D128 q93 = ZERO;
	D128 q94 = ZERO;
	D128 q95 = ZERO;
	D128 q96 = ZERO;
	D128 q97 = ZERO;
	D128 q98 = ZERO;
	D128 q99 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x31 = LOAD(&Ix(a,3,1));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 x21 = LOAD(&Ix(a,2,1));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y31 = LOAD(&Iy(a,3,1));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 y21 = LOAD(&Iy(a,2,1));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z31 = LOAD(&Iz(a,3,1));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    D128 z21 = LOAD(&Iz(a,2,1));
	    q90 = ADD(q90, MUL(MUL(x31, y01), z01));
	    q91 = ADD(q91, MUL(MUL(x01, y31), z01));
	    q92 = ADD(q92, MUL(MUL(x01, y01), z31));
	    q93 = ADD(q93, MUL(MUL(x21, y11), z01));
	    q94 = ADD(q94, MUL(MUL(x21, y01), z11));
	    q95 = ADD(q95, MUL(MUL(x11, y21), z01));
	    q96 = ADD(q96, MUL(MUL(x01, y21), z11));
	    q97 = ADD(q97, MUL(MUL(x11, y01), z21));
	    q98 = ADD(q98, MUL(MUL(x01, y11), z21));
	    q99 = ADD(q99, MUL(MUL(x11, y11), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x31 = LOAD1(&Ix(a,3,1));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 x21 = LOAD1(&Ix(a,2,1));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y31 = LOAD1(&Iy(a,3,1));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 y21 = LOAD1(&Iy(a,2,1));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z31 = LOAD1(&Iz(a,3,1));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    D128 z21 = LOAD1(&Iz(a,2,1));
	    q90 = ADD1(q90, MUL1(MUL1(x31, y01), z01));
	    q91 = ADD1(q91, MUL1(MUL1(x01, y31), z01));
	    q92 = ADD1(q92, MUL1(MUL1(x01, y01), z31));
	    q93 = ADD1(q93, MUL1(MUL1(x21, y11), z01));
	    q94 = ADD1(q94, MUL1(MUL1(x21, y01), z11));
	    q95 = ADD1(q95, MUL1(MUL1(x11, y21), z01));
	    q96 = ADD1(q96, MUL1(MUL1(x01, y21), z11));
	    q97 = ADD1(q97, MUL1(MUL1(x11, y01), z21));
	    q98 = ADD1(q98, MUL1(MUL1(x01, y11), z21));
	    q99 = ADD1(q99, MUL1(MUL1(x11, y11), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK90 = ADD(qK90, MUL(C00, HADD(q90, q91)));
	qK92 = ADD(qK92, MUL(C00, HADD(q92, q93)));
	qK94 = ADD(qK94, MUL(C00, HADD(q94, q95)));
	qK96 = ADD(qK96, MUL(C00, HADD(q96, q97)));
	qK98 = ADD(qK98, MUL(C00, HADD(q98, q99)));

#else // SSE
	    
	// function registers
	T q90 = 0.0;
	T q91 = 0.0;
	T q92 = 0.0;
	T q93 = 0.0;
	T q94 = 0.0;
	T q95 = 0.0;
	T q96 = 0.0;
	T q97 = 0.0;
	T q98 = 0.0;
	T q99 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q90 += Ix(a,3,1)*Iy(a,0,1)*Iz(a,0,1);
	    q91 += Ix(a,0,1)*Iy(a,3,1)*Iz(a,0,1);
	    q92 += Ix(a,0,1)*Iy(a,0,1)*Iz(a,3,1);
	    q93 += Ix(a,2,1)*Iy(a,1,1)*Iz(a,0,1);
	    q94 += Ix(a,2,1)*Iy(a,0,1)*Iz(a,1,1);
	    q95 += Ix(a,1,1)*Iy(a,2,1)*Iz(a,0,1);
	    q96 += Ix(a,0,1)*Iy(a,2,1)*Iz(a,1,1);
	    q97 += Ix(a,1,1)*Iy(a,0,1)*Iz(a,2,1);
	    q98 += Ix(a,0,1)*Iy(a,1,1)*Iz(a,2,1);
	    q99 += Ix(a,1,1)*Iy(a,1,1)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK90 += q90*C[k+0];
	//I[90] += q90*C[k+0];
	I[90] += q90*C_[0];
	//qK91 += q91*C[k+0];
	//I[91] += q91*C[k+0];
	I[91] += q91*C_[0];
	//qK92 += q92*C[k+0];
	//I[92] += q92*C[k+0];
	I[92] += q92*C_[0];
	//qK93 += q93*C[k+0];
	//I[93] += q93*C[k+0];
	I[93] += q93*C_[0];
	//qK94 += q94*C[k+0];
	//I[94] += q94*C[k+0];
	I[94] += q94*C_[0];
	//qK95 += q95*C[k+0];
	//I[95] += q95*C[k+0];
	I[95] += q95*C_[0];
	//qK96 += q96*C[k+0];
	//I[96] += q96*C[k+0];
	I[96] += q96*C_[0];
	//qK97 += q97*C[k+0];
	//I[97] += q97*C[k+0];
	I[97] += q97*C_[0];
	//qK98 += q98*C[k+0];
	//I[98] += q98*C[k+0];
	I[98] += q98*C_[0];
	//qK99 += q99*C[k+0];
	//I[99] += q99*C[k+0];
	I[99] += q99*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK90 = MUL(q, qK90);
	qK92 = MUL(q, qK92);
	qK94 = MUL(q, qK94);
	qK96 = MUL(q, qK96);
	qK98 = MUL(q, qK98);
	num += 10; //num += (fabs(I[98]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[90]) & 0xF) {
	    // 90
	    STOREU(&I[90], ADD(qK90, LOADU(&I[90])));
	    STOREU(&I[92], ADD(qK92, LOADU(&I[92])));
	    STOREU(&I[94], ADD(qK94, LOADU(&I[94])));
	    STOREU(&I[96], ADD(qK96, LOADU(&I[96])));
	    STOREU(&I[98], ADD(qK98, LOADU(&I[98])));
	}
	else {
	    STORE(&I[90], ADD(qK90, LOADU(&I[90])));
	    STORE(&I[92], ADD(qK92, LOADU(&I[92])));
	    STORE(&I[94], ADD(qK94, LOADU(&I[94])));
	    STORE(&I[96], ADD(qK96, LOADU(&I[96])));
	    STORE(&I[98], ADD(qK98, LOADU(&I[98])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[90] += scale*NORMALIZE[10]*NORMALIZE[19]*qK90;
	// num += (fabs(I[90]) >= tol);
	// I[91] += scale*NORMALIZE[11]*NORMALIZE[19]*qK91;
	// num += (fabs(I[91]) >= tol);
	// I[92] += scale*NORMALIZE[12]*NORMALIZE[19]*qK92;
	// num += (fabs(I[92]) >= tol);
	// I[93] += scale*NORMALIZE[13]*NORMALIZE[19]*qK93;
	// num += (fabs(I[93]) >= tol);
	// I[94] += scale*NORMALIZE[14]*NORMALIZE[19]*qK94;
	// num += (fabs(I[94]) >= tol);
	// I[95] += scale*NORMALIZE[15]*NORMALIZE[19]*qK95;
	// num += (fabs(I[95]) >= tol);
	// I[96] += scale*NORMALIZE[16]*NORMALIZE[19]*qK96;
	// num += (fabs(I[96]) >= tol);
	// I[97] += scale*NORMALIZE[17]*NORMALIZE[19]*qK97;
	// num += (fabs(I[97]) >= tol);
	// I[98] += scale*NORMALIZE[18]*NORMALIZE[19]*qK98;
	// num += (fabs(I[98]) >= tol);
	// I[99] += scale*NORMALIZE[19]*NORMALIZE[19]*qK99;
	// num += (fabs(I[99]) >= tol);
    }
    else {
	// I[90] += scale*qK90;
	// num += (fabs(I[90]) >= tol);
	// I[91] += scale*qK91;
	// num += (fabs(I[91]) >= tol);
	// I[92] += scale*qK92;
	// num += (fabs(I[92]) >= tol);
	// I[93] += scale*qK93;
	// num += (fabs(I[93]) >= tol);
	// I[94] += scale*qK94;
	// num += (fabs(I[94]) >= tol);
	// I[95] += scale*qK95;
	// num += (fabs(I[95]) >= tol);
	// I[96] += scale*qK96;
	// num += (fabs(I[96]) >= tol);
	// I[97] += scale*qK97;
	// num += (fabs(I[97]) >= tol);
	// I[98] += scale*qK98;
	// num += (fabs(I[98]) >= tol);
	// I[99] += scale*qK99;
	// num += (fabs(I[99]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <spf| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::SP,rysq::F> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::SP,rysq::F> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 2;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x03 = LOAD(&Ix(a,0,3));
	    D128 x13 = LOAD(&Ix(a,1,3));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y13 = LOAD(&Iy(a,1,3));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y03 = LOAD(&Iy(a,0,3));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z03 = LOAD(&Iz(a,0,3));
	    q0 = ADD(q0, MUL(MUL(x03, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x13, y00), z00));
	    q2 = ADD(q2, MUL(MUL(x03, y10), z00));
	    q3 = ADD(q3, MUL(MUL(x03, y00), z10));
	    q4 = ADD(q4, MUL(MUL(x00, y03), z00));
	    q5 = ADD(q5, MUL(MUL(x10, y03), z00));
	    q6 = ADD(q6, MUL(MUL(x00, y13), z00));
	    q7 = ADD(q7, MUL(MUL(x00, y03), z10));
	    q8 = ADD(q8, MUL(MUL(x00, y00), z03));
	    q9 = ADD(q9, MUL(MUL(x10, y00), z03));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x03 = LOAD1(&Ix(a,0,3));
	    D128 x13 = LOAD1(&Ix(a,1,3));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y13 = LOAD1(&Iy(a,1,3));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y03 = LOAD1(&Iy(a,0,3));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z03 = LOAD1(&Iz(a,0,3));
	    q0 = ADD1(q0, MUL1(MUL1(x03, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x13, y00), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x03, y10), z00));
	    q3 = ADD1(q3, MUL1(MUL1(x03, y00), z10));
	    q4 = ADD1(q4, MUL1(MUL1(x00, y03), z00));
	    q5 = ADD1(q5, MUL1(MUL1(x10, y03), z00));
	    q6 = ADD1(q6, MUL1(MUL1(x00, y13), z00));
	    q7 = ADD1(q7, MUL1(MUL1(x00, y03), z10));
	    q8 = ADD1(q8, MUL1(MUL1(x00, y00), z03));
	    q9 = ADD1(q9, MUL1(MUL1(x10, y00), z03));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C01 = LOADU(&C[k+0]);
	qK0 = ADD(qK0, MUL(C01, HADD(q0, q1)));
	D128 C11 = LOADDUP(&C[k+1]);
	qK2 = ADD(qK2, MUL(C11, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C01, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C11, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C01, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,0,3)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,1,3)*Iy(a,0,0)*Iz(a,0,0);
	    q2 += Ix(a,0,3)*Iy(a,1,0)*Iz(a,0,0);
	    q3 += Ix(a,0,3)*Iy(a,0,0)*Iz(a,1,0);
	    q4 += Ix(a,0,0)*Iy(a,0,3)*Iz(a,0,0);
	    q5 += Ix(a,1,0)*Iy(a,0,3)*Iz(a,0,0);
	    q6 += Ix(a,0,0)*Iy(a,1,3)*Iz(a,0,0);
	    q7 += Ix(a,0,0)*Iy(a,0,3)*Iz(a,1,0);
	    q8 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,0,3);
	    q9 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,0,3);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+1];
	//I[1] += q1*C[k+1];
	I[1] += q1*C_[1];
	//qK2 += q2*C[k+1];
	//I[2] += q2*C[k+1];
	I[2] += q2*C_[1];
	//qK3 += q3*C[k+1];
	//I[3] += q3*C[k+1];
	I[3] += q3*C_[1];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+1];
	//I[5] += q5*C[k+1];
	I[5] += q5*C_[1];
	//qK6 += q6*C[k+1];
	//I[6] += q6*C[k+1];
	I[6] += q6*C_[1];
	//qK7 += q7*C[k+1];
	//I[7] += q7*C[k+1];
	I[7] += q7*C_[1];
	//qK8 += q8*C[k+0];
	//I[8] += q8*C[k+0];
	I[8] += q8*C_[0];
	//qK9 += q9*C[k+1];
	//I[9] += q9*C[k+1];
	I[9] += q9*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[0]*NORMALIZE[10]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[1]*NORMALIZE[10]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[2]*NORMALIZE[10]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[3]*NORMALIZE[10]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[0]*NORMALIZE[11]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[1]*NORMALIZE[11]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[2]*NORMALIZE[11]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[3]*NORMALIZE[11]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[0]*NORMALIZE[12]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[1]*NORMALIZE[12]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK10 = ZERO;
     D128 qK12 = ZERO;
     D128 qK14 = ZERO;
     D128 qK16 = ZERO;
     D128 qK18 = ZERO;
#else
     //double qK10 = 0.0;
     //double qK11 = 0.0;
     //double qK12 = 0.0;
     //double qK13 = 0.0;
     //double qK14 = 0.0;
     //double qK15 = 0.0;
     //double qK16 = 0.0;
     //double qK17 = 0.0;
     //double qK18 = 0.0;
     //double qK19 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q10 = ZERO;
	D128 q11 = ZERO;
	D128 q12 = ZERO;
	D128 q13 = ZERO;
	D128 q14 = ZERO;
	D128 q15 = ZERO;
	D128 q16 = ZERO;
	D128 q17 = ZERO;
	D128 q18 = ZERO;
	D128 q19 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x12 = LOAD(&Ix(a,1,2));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x02 = LOAD(&Ix(a,0,2));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z13 = LOAD(&Iz(a,1,3));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z03 = LOAD(&Iz(a,0,3));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    q10 = ADD(q10, MUL(MUL(x00, y10), z03));
	    q11 = ADD(q11, MUL(MUL(x00, y00), z13));
	    q12 = ADD(q12, MUL(MUL(x02, y01), z00));
	    q13 = ADD(q13, MUL(MUL(x12, y01), z00));
	    q14 = ADD(q14, MUL(MUL(x02, y11), z00));
	    q15 = ADD(q15, MUL(MUL(x02, y01), z10));
	    q16 = ADD(q16, MUL(MUL(x02, y00), z01));
	    q17 = ADD(q17, MUL(MUL(x12, y00), z01));
	    q18 = ADD(q18, MUL(MUL(x02, y10), z01));
	    q19 = ADD(q19, MUL(MUL(x02, y00), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x12 = LOAD1(&Ix(a,1,2));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x02 = LOAD1(&Ix(a,0,2));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z13 = LOAD1(&Iz(a,1,3));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z03 = LOAD1(&Iz(a,0,3));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    q10 = ADD1(q10, MUL1(MUL1(x00, y10), z03));
	    q11 = ADD1(q11, MUL1(MUL1(x00, y00), z13));
	    q12 = ADD1(q12, MUL1(MUL1(x02, y01), z00));
	    q13 = ADD1(q13, MUL1(MUL1(x12, y01), z00));
	    q14 = ADD1(q14, MUL1(MUL1(x02, y11), z00));
	    q15 = ADD1(q15, MUL1(MUL1(x02, y01), z10));
	    q16 = ADD1(q16, MUL1(MUL1(x02, y00), z01));
	    q17 = ADD1(q17, MUL1(MUL1(x12, y00), z01));
	    q18 = ADD1(q18, MUL1(MUL1(x02, y10), z01));
	    q19 = ADD1(q19, MUL1(MUL1(x02, y00), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C11 = LOADDUP(&C[k+1]);
	qK10 = ADD(qK10, MUL(C11, HADD(q10, q11)));
	D128 C01 = LOADU(&C[k+0]);
	qK12 = ADD(qK12, MUL(C01, HADD(q12, q13)));
	qK14 = ADD(qK14, MUL(C11, HADD(q14, q15)));
	qK16 = ADD(qK16, MUL(C01, HADD(q16, q17)));
	qK18 = ADD(qK18, MUL(C11, HADD(q18, q19)));

#else // SSE
	    
	// function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
	T q16 = 0.0;
	T q17 = 0.0;
	T q18 = 0.0;
	T q19 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,0,3);
	    q11 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,1,3);
	    q12 += Ix(a,0,2)*Iy(a,0,1)*Iz(a,0,0);
	    q13 += Ix(a,1,2)*Iy(a,0,1)*Iz(a,0,0);
	    q14 += Ix(a,0,2)*Iy(a,1,1)*Iz(a,0,0);
	    q15 += Ix(a,0,2)*Iy(a,0,1)*Iz(a,1,0);
	    q16 += Ix(a,0,2)*Iy(a,0,0)*Iz(a,0,1);
	    q17 += Ix(a,1,2)*Iy(a,0,0)*Iz(a,0,1);
	    q18 += Ix(a,0,2)*Iy(a,1,0)*Iz(a,0,1);
	    q19 += Ix(a,0,2)*Iy(a,0,0)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK10 += q10*C[k+1];
	//I[10] += q10*C[k+1];
	I[10] += q10*C_[1];
	//qK11 += q11*C[k+1];
	//I[11] += q11*C[k+1];
	I[11] += q11*C_[1];
	//qK12 += q12*C[k+0];
	//I[12] += q12*C[k+0];
	I[12] += q12*C_[0];
	//qK13 += q13*C[k+1];
	//I[13] += q13*C[k+1];
	I[13] += q13*C_[1];
	//qK14 += q14*C[k+1];
	//I[14] += q14*C[k+1];
	I[14] += q14*C_[1];
	//qK15 += q15*C[k+1];
	//I[15] += q15*C[k+1];
	I[15] += q15*C_[1];
	//qK16 += q16*C[k+0];
	//I[16] += q16*C[k+0];
	I[16] += q16*C_[0];
	//qK17 += q17*C[k+1];
	//I[17] += q17*C[k+1];
	I[17] += q17*C_[1];
	//qK18 += q18*C[k+1];
	//I[18] += q18*C[k+1];
	I[18] += q18*C_[1];
	//qK19 += q19*C[k+1];
	//I[19] += q19*C[k+1];
	I[19] += q19*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK10 = MUL(q, qK10);
	qK12 = MUL(q, qK12);
	qK14 = MUL(q, qK14);
	qK16 = MUL(q, qK16);
	qK18 = MUL(q, qK18);
	num += 10; //num += (fabs(I[18]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[10]) & 0xF) {
	    // 10
	    STOREU(&I[10], ADD(qK10, LOADU(&I[10])));
	    STOREU(&I[12], ADD(qK12, LOADU(&I[12])));
	    STOREU(&I[14], ADD(qK14, LOADU(&I[14])));
	    STOREU(&I[16], ADD(qK16, LOADU(&I[16])));
	    STOREU(&I[18], ADD(qK18, LOADU(&I[18])));
	}
	else {
	    STORE(&I[10], ADD(qK10, LOADU(&I[10])));
	    STORE(&I[12], ADD(qK12, LOADU(&I[12])));
	    STORE(&I[14], ADD(qK14, LOADU(&I[14])));
	    STORE(&I[16], ADD(qK16, LOADU(&I[16])));
	    STORE(&I[18], ADD(qK18, LOADU(&I[18])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[10] += scale*NORMALIZE[2]*NORMALIZE[12]*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*NORMALIZE[3]*NORMALIZE[12]*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*NORMALIZE[0]*NORMALIZE[13]*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*NORMALIZE[1]*NORMALIZE[13]*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*NORMALIZE[2]*NORMALIZE[13]*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*NORMALIZE[3]*NORMALIZE[13]*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*NORMALIZE[0]*NORMALIZE[14]*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*NORMALIZE[1]*NORMALIZE[14]*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*NORMALIZE[2]*NORMALIZE[14]*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*NORMALIZE[3]*NORMALIZE[14]*qK19;
	// num += (fabs(I[19]) >= tol);
    }
    else {
	// I[10] += scale*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*qK19;
	// num += (fabs(I[19]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK20 = ZERO;
     D128 qK22 = ZERO;
     D128 qK24 = ZERO;
     D128 qK26 = ZERO;
     D128 qK28 = ZERO;
#else
     //double qK20 = 0.0;
     //double qK21 = 0.0;
     //double qK22 = 0.0;
     //double qK23 = 0.0;
     //double qK24 = 0.0;
     //double qK25 = 0.0;
     //double qK26 = 0.0;
     //double qK27 = 0.0;
     //double qK28 = 0.0;
     //double qK29 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q20 = ZERO;
	D128 q21 = ZERO;
	D128 q22 = ZERO;
	D128 q23 = ZERO;
	D128 q24 = ZERO;
	D128 q25 = ZERO;
	D128 q26 = ZERO;
	D128 q27 = ZERO;
	D128 q28 = ZERO;
	D128 q29 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 y12 = LOAD(&Iy(a,1,2));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y02 = LOAD(&Iy(a,0,2));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    D128 z02 = LOAD(&Iz(a,0,2));
	    q20 = ADD(q20, MUL(MUL(x01, y02), z00));
	    q21 = ADD(q21, MUL(MUL(x11, y02), z00));
	    q22 = ADD(q22, MUL(MUL(x01, y12), z00));
	    q23 = ADD(q23, MUL(MUL(x01, y02), z10));
	    q24 = ADD(q24, MUL(MUL(x00, y02), z01));
	    q25 = ADD(q25, MUL(MUL(x10, y02), z01));
	    q26 = ADD(q26, MUL(MUL(x00, y12), z01));
	    q27 = ADD(q27, MUL(MUL(x00, y02), z11));
	    q28 = ADD(q28, MUL(MUL(x01, y00), z02));
	    q29 = ADD(q29, MUL(MUL(x11, y00), z02));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 y12 = LOAD1(&Iy(a,1,2));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y02 = LOAD1(&Iy(a,0,2));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    D128 z02 = LOAD1(&Iz(a,0,2));
	    q20 = ADD1(q20, MUL1(MUL1(x01, y02), z00));
	    q21 = ADD1(q21, MUL1(MUL1(x11, y02), z00));
	    q22 = ADD1(q22, MUL1(MUL1(x01, y12), z00));
	    q23 = ADD1(q23, MUL1(MUL1(x01, y02), z10));
	    q24 = ADD1(q24, MUL1(MUL1(x00, y02), z01));
	    q25 = ADD1(q25, MUL1(MUL1(x10, y02), z01));
	    q26 = ADD1(q26, MUL1(MUL1(x00, y12), z01));
	    q27 = ADD1(q27, MUL1(MUL1(x00, y02), z11));
	    q28 = ADD1(q28, MUL1(MUL1(x01, y00), z02));
	    q29 = ADD1(q29, MUL1(MUL1(x11, y00), z02));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C01 = LOADU(&C[k+0]);
	qK20 = ADD(qK20, MUL(C01, HADD(q20, q21)));
	D128 C11 = LOADDUP(&C[k+1]);
	qK22 = ADD(qK22, MUL(C11, HADD(q22, q23)));
	qK24 = ADD(qK24, MUL(C01, HADD(q24, q25)));
	qK26 = ADD(qK26, MUL(C11, HADD(q26, q27)));
	qK28 = ADD(qK28, MUL(C01, HADD(q28, q29)));

#else // SSE
	    
	// function registers
	T q20 = 0.0;
	T q21 = 0.0;
	T q22 = 0.0;
	T q23 = 0.0;
	T q24 = 0.0;
	T q25 = 0.0;
	T q26 = 0.0;
	T q27 = 0.0;
	T q28 = 0.0;
	T q29 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q20 += Ix(a,0,1)*Iy(a,0,2)*Iz(a,0,0);
	    q21 += Ix(a,1,1)*Iy(a,0,2)*Iz(a,0,0);
	    q22 += Ix(a,0,1)*Iy(a,1,2)*Iz(a,0,0);
	    q23 += Ix(a,0,1)*Iy(a,0,2)*Iz(a,1,0);
	    q24 += Ix(a,0,0)*Iy(a,0,2)*Iz(a,0,1);
	    q25 += Ix(a,1,0)*Iy(a,0,2)*Iz(a,0,1);
	    q26 += Ix(a,0,0)*Iy(a,1,2)*Iz(a,0,1);
	    q27 += Ix(a,0,0)*Iy(a,0,2)*Iz(a,1,1);
	    q28 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,0,2);
	    q29 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,0,2);
	}
	    
	//contraction coefficients
	//qK20 += q20*C[k+0];
	//I[20] += q20*C[k+0];
	I[20] += q20*C_[0];
	//qK21 += q21*C[k+1];
	//I[21] += q21*C[k+1];
	I[21] += q21*C_[1];
	//qK22 += q22*C[k+1];
	//I[22] += q22*C[k+1];
	I[22] += q22*C_[1];
	//qK23 += q23*C[k+1];
	//I[23] += q23*C[k+1];
	I[23] += q23*C_[1];
	//qK24 += q24*C[k+0];
	//I[24] += q24*C[k+0];
	I[24] += q24*C_[0];
	//qK25 += q25*C[k+1];
	//I[25] += q25*C[k+1];
	I[25] += q25*C_[1];
	//qK26 += q26*C[k+1];
	//I[26] += q26*C[k+1];
	I[26] += q26*C_[1];
	//qK27 += q27*C[k+1];
	//I[27] += q27*C[k+1];
	I[27] += q27*C_[1];
	//qK28 += q28*C[k+0];
	//I[28] += q28*C[k+0];
	I[28] += q28*C_[0];
	//qK29 += q29*C[k+1];
	//I[29] += q29*C[k+1];
	I[29] += q29*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK20 = MUL(q, qK20);
	qK22 = MUL(q, qK22);
	qK24 = MUL(q, qK24);
	qK26 = MUL(q, qK26);
	qK28 = MUL(q, qK28);
	num += 10; //num += (fabs(I[28]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[20]) & 0xF) {
	    // 20
	    STOREU(&I[20], ADD(qK20, LOADU(&I[20])));
	    STOREU(&I[22], ADD(qK22, LOADU(&I[22])));
	    STOREU(&I[24], ADD(qK24, LOADU(&I[24])));
	    STOREU(&I[26], ADD(qK26, LOADU(&I[26])));
	    STOREU(&I[28], ADD(qK28, LOADU(&I[28])));
	}
	else {
	    STORE(&I[20], ADD(qK20, LOADU(&I[20])));
	    STORE(&I[22], ADD(qK22, LOADU(&I[22])));
	    STORE(&I[24], ADD(qK24, LOADU(&I[24])));
	    STORE(&I[26], ADD(qK26, LOADU(&I[26])));
	    STORE(&I[28], ADD(qK28, LOADU(&I[28])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[20] += scale*NORMALIZE[0]*NORMALIZE[15]*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*NORMALIZE[1]*NORMALIZE[15]*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*NORMALIZE[2]*NORMALIZE[15]*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*NORMALIZE[3]*NORMALIZE[15]*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*NORMALIZE[0]*NORMALIZE[16]*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*NORMALIZE[1]*NORMALIZE[16]*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*NORMALIZE[2]*NORMALIZE[16]*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*NORMALIZE[3]*NORMALIZE[16]*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*NORMALIZE[0]*NORMALIZE[17]*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*NORMALIZE[1]*NORMALIZE[17]*qK29;
	// num += (fabs(I[29]) >= tol);
    }
    else {
	// I[20] += scale*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*qK29;
	// num += (fabs(I[29]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK30 = ZERO;
     D128 qK32 = ZERO;
     D128 qK34 = ZERO;
     D128 qK36 = ZERO;
     D128 qK38 = ZERO;
#else
     //double qK30 = 0.0;
     //double qK31 = 0.0;
     //double qK32 = 0.0;
     //double qK33 = 0.0;
     //double qK34 = 0.0;
     //double qK35 = 0.0;
     //double qK36 = 0.0;
     //double qK37 = 0.0;
     //double qK38 = 0.0;
     //double qK39 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q30 = ZERO;
	D128 q31 = ZERO;
	D128 q32 = ZERO;
	D128 q33 = ZERO;
	D128 q34 = ZERO;
	D128 q35 = ZERO;
	D128 q36 = ZERO;
	D128 q37 = ZERO;
	D128 q38 = ZERO;
	D128 q39 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z12 = LOAD(&Iz(a,1,2));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z02 = LOAD(&Iz(a,0,2));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    q30 = ADD(q30, MUL(MUL(x01, y10), z02));
	    q31 = ADD(q31, MUL(MUL(x01, y00), z12));
	    q32 = ADD(q32, MUL(MUL(x00, y01), z02));
	    q33 = ADD(q33, MUL(MUL(x10, y01), z02));
	    q34 = ADD(q34, MUL(MUL(x00, y11), z02));
	    q35 = ADD(q35, MUL(MUL(x00, y01), z12));
	    q36 = ADD(q36, MUL(MUL(x01, y01), z01));
	    q37 = ADD(q37, MUL(MUL(x11, y01), z01));
	    q38 = ADD(q38, MUL(MUL(x01, y11), z01));
	    q39 = ADD(q39, MUL(MUL(x01, y01), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z12 = LOAD1(&Iz(a,1,2));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z02 = LOAD1(&Iz(a,0,2));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    q30 = ADD1(q30, MUL1(MUL1(x01, y10), z02));
	    q31 = ADD1(q31, MUL1(MUL1(x01, y00), z12));
	    q32 = ADD1(q32, MUL1(MUL1(x00, y01), z02));
	    q33 = ADD1(q33, MUL1(MUL1(x10, y01), z02));
	    q34 = ADD1(q34, MUL1(MUL1(x00, y11), z02));
	    q35 = ADD1(q35, MUL1(MUL1(x00, y01), z12));
	    q36 = ADD1(q36, MUL1(MUL1(x01, y01), z01));
	    q37 = ADD1(q37, MUL1(MUL1(x11, y01), z01));
	    q38 = ADD1(q38, MUL1(MUL1(x01, y11), z01));
	    q39 = ADD1(q39, MUL1(MUL1(x01, y01), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C11 = LOADDUP(&C[k+1]);
	qK30 = ADD(qK30, MUL(C11, HADD(q30, q31)));
	D128 C01 = LOADU(&C[k+0]);
	qK32 = ADD(qK32, MUL(C01, HADD(q32, q33)));
	qK34 = ADD(qK34, MUL(C11, HADD(q34, q35)));
	qK36 = ADD(qK36, MUL(C01, HADD(q36, q37)));
	qK38 = ADD(qK38, MUL(C11, HADD(q38, q39)));

#else // SSE
	    
	// function registers
	T q30 = 0.0;
	T q31 = 0.0;
	T q32 = 0.0;
	T q33 = 0.0;
	T q34 = 0.0;
	T q35 = 0.0;
	T q36 = 0.0;
	T q37 = 0.0;
	T q38 = 0.0;
	T q39 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q30 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,0,2);
	    q31 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,1,2);
	    q32 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,0,2);
	    q33 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,0,2);
	    q34 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,0,2);
	    q35 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,1,2);
	    q36 += Ix(a,0,1)*Iy(a,0,1)*Iz(a,0,1);
	    q37 += Ix(a,1,1)*Iy(a,0,1)*Iz(a,0,1);
	    q38 += Ix(a,0,1)*Iy(a,1,1)*Iz(a,0,1);
	    q39 += Ix(a,0,1)*Iy(a,0,1)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK30 += q30*C[k+1];
	//I[30] += q30*C[k+1];
	I[30] += q30*C_[1];
	//qK31 += q31*C[k+1];
	//I[31] += q31*C[k+1];
	I[31] += q31*C_[1];
	//qK32 += q32*C[k+0];
	//I[32] += q32*C[k+0];
	I[32] += q32*C_[0];
	//qK33 += q33*C[k+1];
	//I[33] += q33*C[k+1];
	I[33] += q33*C_[1];
	//qK34 += q34*C[k+1];
	//I[34] += q34*C[k+1];
	I[34] += q34*C_[1];
	//qK35 += q35*C[k+1];
	//I[35] += q35*C[k+1];
	I[35] += q35*C_[1];
	//qK36 += q36*C[k+0];
	//I[36] += q36*C[k+0];
	I[36] += q36*C_[0];
	//qK37 += q37*C[k+1];
	//I[37] += q37*C[k+1];
	I[37] += q37*C_[1];
	//qK38 += q38*C[k+1];
	//I[38] += q38*C[k+1];
	I[38] += q38*C_[1];
	//qK39 += q39*C[k+1];
	//I[39] += q39*C[k+1];
	I[39] += q39*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK30 = MUL(q, qK30);
	qK32 = MUL(q, qK32);
	qK34 = MUL(q, qK34);
	qK36 = MUL(q, qK36);
	qK38 = MUL(q, qK38);
	num += 10; //num += (fabs(I[38]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[30]) & 0xF) {
	    // 30
	    STOREU(&I[30], ADD(qK30, LOADU(&I[30])));
	    STOREU(&I[32], ADD(qK32, LOADU(&I[32])));
	    STOREU(&I[34], ADD(qK34, LOADU(&I[34])));
	    STOREU(&I[36], ADD(qK36, LOADU(&I[36])));
	    STOREU(&I[38], ADD(qK38, LOADU(&I[38])));
	}
	else {
	    STORE(&I[30], ADD(qK30, LOADU(&I[30])));
	    STORE(&I[32], ADD(qK32, LOADU(&I[32])));
	    STORE(&I[34], ADD(qK34, LOADU(&I[34])));
	    STORE(&I[36], ADD(qK36, LOADU(&I[36])));
	    STORE(&I[38], ADD(qK38, LOADU(&I[38])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[30] += scale*NORMALIZE[2]*NORMALIZE[17]*qK30;
	// num += (fabs(I[30]) >= tol);
	// I[31] += scale*NORMALIZE[3]*NORMALIZE[17]*qK31;
	// num += (fabs(I[31]) >= tol);
	// I[32] += scale*NORMALIZE[0]*NORMALIZE[18]*qK32;
	// num += (fabs(I[32]) >= tol);
	// I[33] += scale*NORMALIZE[1]*NORMALIZE[18]*qK33;
	// num += (fabs(I[33]) >= tol);
	// I[34] += scale*NORMALIZE[2]*NORMALIZE[18]*qK34;
	// num += (fabs(I[34]) >= tol);
	// I[35] += scale*NORMALIZE[3]*NORMALIZE[18]*qK35;
	// num += (fabs(I[35]) >= tol);
	// I[36] += scale*NORMALIZE[0]*NORMALIZE[19]*qK36;
	// num += (fabs(I[36]) >= tol);
	// I[37] += scale*NORMALIZE[1]*NORMALIZE[19]*qK37;
	// num += (fabs(I[37]) >= tol);
	// I[38] += scale*NORMALIZE[2]*NORMALIZE[19]*qK38;
	// num += (fabs(I[38]) >= tol);
	// I[39] += scale*NORMALIZE[3]*NORMALIZE[19]*qK39;
	// num += (fabs(I[39]) >= tol);
    }
    else {
	// I[30] += scale*qK30;
	// num += (fabs(I[30]) >= tol);
	// I[31] += scale*qK31;
	// num += (fabs(I[31]) >= tol);
	// I[32] += scale*qK32;
	// num += (fabs(I[32]) >= tol);
	// I[33] += scale*qK33;
	// num += (fabs(I[33]) >= tol);
	// I[34] += scale*qK34;
	// num += (fabs(I[34]) >= tol);
	// I[35] += scale*qK35;
	// num += (fabs(I[35]) >= tol);
	// I[36] += scale*qK36;
	// num += (fabs(I[36]) >= tol);
	// I[37] += scale*qK37;
	// num += (fabs(I[37]) >= tol);
	// I[38] += scale*qK38;
	// num += (fabs(I[38]) >= tol);
	// I[39] += scale*qK39;
	// num += (fabs(I[39]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <ssp| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::S,rysq::SP> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::S,rysq::SP> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 1;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x00, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x01, y00), z00));
	    q2 = ADD(q2, MUL(MUL(x00, y01), z00));
	    q3 = ADD(q3, MUL(MUL(x00, y00), z01));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x00, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x01, y00), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x00, y01), z00));
	    q3 = ADD1(q3, MUL1(MUL1(x00, y00), z01));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C01 = LOADU(&C[k+0]);
	qK0 = ADD(qK0, MUL(C01, HADD(q0, q1)));
	D128 C11 = LOADDUP(&C[k+1]);
	qK2 = ADD(qK2, MUL(C11, HADD(q2, q3)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,0,0);
	    q2 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,0,0);
	    q3 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,0,1);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+1];
	//I[1] += q1*C[k+1];
	I[1] += q1*C_[1];
	//qK2 += q2*C[k+1];
	//I[2] += q2*C[k+1];
	I[2] += q2*C_[1];
	//qK3 += q3*C[k+1];
	//I[3] += q3*C[k+1];
	I[3] += q3*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	num += 4; //num += (fabs(I[2]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[0]*NORMALIZE[0]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[0]*NORMALIZE[1]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[0]*NORMALIZE[2]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[0]*NORMALIZE[3]*qK3;
	// num += (fabs(I[3]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <psp| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::P,rysq::SP> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::P,rysq::SP> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 2;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x10, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x00, y10), z00));
	    q2 = ADD(q2, MUL(MUL(x00, y00), z10));
	    q3 = ADD(q3, MUL(MUL(x11, y00), z00));
	    q4 = ADD(q4, MUL(MUL(x01, y10), z00));
	    q5 = ADD(q5, MUL(MUL(x01, y00), z10));
	    q6 = ADD(q6, MUL(MUL(x10, y01), z00));
	    q7 = ADD(q7, MUL(MUL(x00, y11), z00));
	    q8 = ADD(q8, MUL(MUL(x00, y01), z10));
	    q9 = ADD(q9, MUL(MUL(x10, y00), z01));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x10, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x00, y10), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x00, y00), z10));
	    q3 = ADD1(q3, MUL1(MUL1(x11, y00), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x01, y10), z00));
	    q5 = ADD1(q5, MUL1(MUL1(x01, y00), z10));
	    q6 = ADD1(q6, MUL1(MUL1(x10, y01), z00));
	    q7 = ADD1(q7, MUL1(MUL1(x00, y11), z00));
	    q8 = ADD1(q8, MUL1(MUL1(x00, y01), z10));
	    q9 = ADD1(q9, MUL1(MUL1(x10, y00), z01));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	D128 C01 = LOADU(&C[k+0]);
	qK2 = ADD(qK2, MUL(C01, HADD(q2, q3)));
	D128 C11 = LOADDUP(&C[k+1]);
	qK4 = ADD(qK4, MUL(C11, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C11, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C11, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,0,0);
	    q2 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,1,0);
	    q3 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,0,0);
	    q4 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,0,0);
	    q5 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,1,0);
	    q6 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,0,0);
	    q7 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,0,0);
	    q8 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,1,0);
	    q9 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,0,1);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+1];
	//I[3] += q3*C[k+1];
	I[3] += q3*C_[1];
	//qK4 += q4*C[k+1];
	//I[4] += q4*C[k+1];
	I[4] += q4*C_[1];
	//qK5 += q5*C[k+1];
	//I[5] += q5*C[k+1];
	I[5] += q5*C_[1];
	//qK6 += q6*C[k+1];
	//I[6] += q6*C[k+1];
	I[6] += q6*C_[1];
	//qK7 += q7*C[k+1];
	//I[7] += q7*C[k+1];
	I[7] += q7*C_[1];
	//qK8 += q8*C[k+1];
	//I[8] += q8*C[k+1];
	I[8] += q8*C_[1];
	//qK9 += q9*C[k+1];
	//I[9] += q9*C[k+1];
	I[9] += q9*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[1]*NORMALIZE[0]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[2]*NORMALIZE[0]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[3]*NORMALIZE[0]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[1]*NORMALIZE[1]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[2]*NORMALIZE[1]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[3]*NORMALIZE[1]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[1]*NORMALIZE[2]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[2]*NORMALIZE[2]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[3]*NORMALIZE[2]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[1]*NORMALIZE[3]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK10 = ZERO;
#else
     //double qK10 = 0.0;
     //double qK11 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q10 = ZERO;
	D128 q11 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    q10 = ADD(q10, MUL(MUL(x00, y10), z01));
	    q11 = ADD(q11, MUL(MUL(x00, y00), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    q10 = ADD1(q10, MUL1(MUL1(x00, y10), z01));
	    q11 = ADD1(q11, MUL1(MUL1(x00, y00), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C11 = LOADDUP(&C[k+1]);
	qK10 = ADD(qK10, MUL(C11, HADD(q10, q11)));

#else // SSE
	    
	// function registers
	T q10 = 0.0;
	T q11 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,0,1);
	    q11 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK10 += q10*C[k+1];
	//I[10] += q10*C[k+1];
	I[10] += q10*C_[1];
	//qK11 += q11*C[k+1];
	//I[11] += q11*C[k+1];
	I[11] += q11*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK10 = MUL(q, qK10);
	num += 2; //num += (fabs(I[10]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[10]) & 0xF) {
	    // 10
	    STOREU(&I[10], ADD(qK10, LOADU(&I[10])));
	}
	else {
	    STORE(&I[10], ADD(qK10, LOADU(&I[10])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[10] += scale*NORMALIZE[2]*NORMALIZE[3]*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*NORMALIZE[3]*NORMALIZE[3]*qK11;
	// num += (fabs(I[11]) >= tol);
    }
    else {
	// I[10] += scale*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*qK11;
	// num += (fabs(I[11]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <dsp| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::D,rysq::SP> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::D,rysq::SP> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 3;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x21 = LOAD(&Ix(a,2,1));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x20, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x00, y20), z00));
	    q2 = ADD(q2, MUL(MUL(x00, y00), z20));
	    q3 = ADD(q3, MUL(MUL(x10, y10), z00));
	    q4 = ADD(q4, MUL(MUL(x10, y00), z10));
	    q5 = ADD(q5, MUL(MUL(x00, y10), z10));
	    q6 = ADD(q6, MUL(MUL(x21, y00), z00));
	    q7 = ADD(q7, MUL(MUL(x01, y20), z00));
	    q8 = ADD(q8, MUL(MUL(x01, y00), z20));
	    q9 = ADD(q9, MUL(MUL(x11, y10), z00));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x21 = LOAD1(&Ix(a,2,1));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x20, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x00, y20), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x00, y00), z20));
	    q3 = ADD1(q3, MUL1(MUL1(x10, y10), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x10, y00), z10));
	    q5 = ADD1(q5, MUL1(MUL1(x00, y10), z10));
	    q6 = ADD1(q6, MUL1(MUL1(x21, y00), z00));
	    q7 = ADD1(q7, MUL1(MUL1(x01, y20), z00));
	    q8 = ADD1(q8, MUL1(MUL1(x01, y00), z20));
	    q9 = ADD1(q9, MUL1(MUL1(x11, y10), z00));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD(qK2, MUL(C00, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C00, HADD(q4, q5)));
	D128 C11 = LOADDUP(&C[k+1]);
	qK6 = ADD(qK6, MUL(C11, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C11, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,2,0)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,0)*Iy(a,2,0)*Iz(a,0,0);
	    q2 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,2,0);
	    q3 += Ix(a,1,0)*Iy(a,1,0)*Iz(a,0,0);
	    q4 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,1,0);
	    q5 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,1,0);
	    q6 += Ix(a,2,1)*Iy(a,0,0)*Iz(a,0,0);
	    q7 += Ix(a,0,1)*Iy(a,2,0)*Iz(a,0,0);
	    q8 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,2,0);
	    q9 += Ix(a,1,1)*Iy(a,1,0)*Iz(a,0,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+0];
	//I[3] += q3*C[k+0];
	I[3] += q3*C_[0];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+0];
	//I[5] += q5*C[k+0];
	I[5] += q5*C_[0];
	//qK6 += q6*C[k+1];
	//I[6] += q6*C[k+1];
	I[6] += q6*C_[1];
	//qK7 += q7*C[k+1];
	//I[7] += q7*C[k+1];
	I[7] += q7*C_[1];
	//qK8 += q8*C[k+1];
	//I[8] += q8*C[k+1];
	I[8] += q8*C_[1];
	//qK9 += q9*C[k+1];
	//I[9] += q9*C[k+1];
	I[9] += q9*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[4]*NORMALIZE[0]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[5]*NORMALIZE[0]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[6]*NORMALIZE[0]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[7]*NORMALIZE[0]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[8]*NORMALIZE[0]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[9]*NORMALIZE[0]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[4]*NORMALIZE[1]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[5]*NORMALIZE[1]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[6]*NORMALIZE[1]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[7]*NORMALIZE[1]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK10 = ZERO;
     D128 qK12 = ZERO;
     D128 qK14 = ZERO;
     D128 qK16 = ZERO;
     D128 qK18 = ZERO;
#else
     //double qK10 = 0.0;
     //double qK11 = 0.0;
     //double qK12 = 0.0;
     //double qK13 = 0.0;
     //double qK14 = 0.0;
     //double qK15 = 0.0;
     //double qK16 = 0.0;
     //double qK17 = 0.0;
     //double qK18 = 0.0;
     //double qK19 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q10 = ZERO;
	D128 q11 = ZERO;
	D128 q12 = ZERO;
	D128 q13 = ZERO;
	D128 q14 = ZERO;
	D128 q15 = ZERO;
	D128 q16 = ZERO;
	D128 q17 = ZERO;
	D128 q18 = ZERO;
	D128 q19 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y21 = LOAD(&Iy(a,2,1));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q10 = ADD(q10, MUL(MUL(x11, y00), z10));
	    q11 = ADD(q11, MUL(MUL(x01, y10), z10));
	    q12 = ADD(q12, MUL(MUL(x20, y01), z00));
	    q13 = ADD(q13, MUL(MUL(x00, y21), z00));
	    q14 = ADD(q14, MUL(MUL(x00, y01), z20));
	    q15 = ADD(q15, MUL(MUL(x10, y11), z00));
	    q16 = ADD(q16, MUL(MUL(x10, y01), z10));
	    q17 = ADD(q17, MUL(MUL(x00, y11), z10));
	    q18 = ADD(q18, MUL(MUL(x20, y00), z01));
	    q19 = ADD(q19, MUL(MUL(x00, y20), z01));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y21 = LOAD1(&Iy(a,2,1));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q10 = ADD1(q10, MUL1(MUL1(x11, y00), z10));
	    q11 = ADD1(q11, MUL1(MUL1(x01, y10), z10));
	    q12 = ADD1(q12, MUL1(MUL1(x20, y01), z00));
	    q13 = ADD1(q13, MUL1(MUL1(x00, y21), z00));
	    q14 = ADD1(q14, MUL1(MUL1(x00, y01), z20));
	    q15 = ADD1(q15, MUL1(MUL1(x10, y11), z00));
	    q16 = ADD1(q16, MUL1(MUL1(x10, y01), z10));
	    q17 = ADD1(q17, MUL1(MUL1(x00, y11), z10));
	    q18 = ADD1(q18, MUL1(MUL1(x20, y00), z01));
	    q19 = ADD1(q19, MUL1(MUL1(x00, y20), z01));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C11 = LOADDUP(&C[k+1]);
	qK10 = ADD(qK10, MUL(C11, HADD(q10, q11)));
	qK12 = ADD(qK12, MUL(C11, HADD(q12, q13)));
	qK14 = ADD(qK14, MUL(C11, HADD(q14, q15)));
	qK16 = ADD(qK16, MUL(C11, HADD(q16, q17)));
	qK18 = ADD(qK18, MUL(C11, HADD(q18, q19)));

#else // SSE
	    
	// function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
	T q16 = 0.0;
	T q17 = 0.0;
	T q18 = 0.0;
	T q19 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,1,0);
	    q11 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,1,0);
	    q12 += Ix(a,2,0)*Iy(a,0,1)*Iz(a,0,0);
	    q13 += Ix(a,0,0)*Iy(a,2,1)*Iz(a,0,0);
	    q14 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,2,0);
	    q15 += Ix(a,1,0)*Iy(a,1,1)*Iz(a,0,0);
	    q16 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,1,0);
	    q17 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,1,0);
	    q18 += Ix(a,2,0)*Iy(a,0,0)*Iz(a,0,1);
	    q19 += Ix(a,0,0)*Iy(a,2,0)*Iz(a,0,1);
	}
	    
	//contraction coefficients
	//qK10 += q10*C[k+1];
	//I[10] += q10*C[k+1];
	I[10] += q10*C_[1];
	//qK11 += q11*C[k+1];
	//I[11] += q11*C[k+1];
	I[11] += q11*C_[1];
	//qK12 += q12*C[k+1];
	//I[12] += q12*C[k+1];
	I[12] += q12*C_[1];
	//qK13 += q13*C[k+1];
	//I[13] += q13*C[k+1];
	I[13] += q13*C_[1];
	//qK14 += q14*C[k+1];
	//I[14] += q14*C[k+1];
	I[14] += q14*C_[1];
	//qK15 += q15*C[k+1];
	//I[15] += q15*C[k+1];
	I[15] += q15*C_[1];
	//qK16 += q16*C[k+1];
	//I[16] += q16*C[k+1];
	I[16] += q16*C_[1];
	//qK17 += q17*C[k+1];
	//I[17] += q17*C[k+1];
	I[17] += q17*C_[1];
	//qK18 += q18*C[k+1];
	//I[18] += q18*C[k+1];
	I[18] += q18*C_[1];
	//qK19 += q19*C[k+1];
	//I[19] += q19*C[k+1];
	I[19] += q19*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK10 = MUL(q, qK10);
	qK12 = MUL(q, qK12);
	qK14 = MUL(q, qK14);
	qK16 = MUL(q, qK16);
	qK18 = MUL(q, qK18);
	num += 10; //num += (fabs(I[18]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[10]) & 0xF) {
	    // 10
	    STOREU(&I[10], ADD(qK10, LOADU(&I[10])));
	    STOREU(&I[12], ADD(qK12, LOADU(&I[12])));
	    STOREU(&I[14], ADD(qK14, LOADU(&I[14])));
	    STOREU(&I[16], ADD(qK16, LOADU(&I[16])));
	    STOREU(&I[18], ADD(qK18, LOADU(&I[18])));
	}
	else {
	    STORE(&I[10], ADD(qK10, LOADU(&I[10])));
	    STORE(&I[12], ADD(qK12, LOADU(&I[12])));
	    STORE(&I[14], ADD(qK14, LOADU(&I[14])));
	    STORE(&I[16], ADD(qK16, LOADU(&I[16])));
	    STORE(&I[18], ADD(qK18, LOADU(&I[18])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[10] += scale*NORMALIZE[8]*NORMALIZE[1]*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*NORMALIZE[9]*NORMALIZE[1]*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*NORMALIZE[4]*NORMALIZE[2]*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*NORMALIZE[5]*NORMALIZE[2]*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*NORMALIZE[6]*NORMALIZE[2]*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*NORMALIZE[7]*NORMALIZE[2]*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*NORMALIZE[8]*NORMALIZE[2]*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*NORMALIZE[9]*NORMALIZE[2]*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*NORMALIZE[4]*NORMALIZE[3]*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*NORMALIZE[5]*NORMALIZE[3]*qK19;
	// num += (fabs(I[19]) >= tol);
    }
    else {
	// I[10] += scale*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*qK19;
	// num += (fabs(I[19]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK20 = ZERO;
     D128 qK22 = ZERO;
#else
     //double qK20 = 0.0;
     //double qK21 = 0.0;
     //double qK22 = 0.0;
     //double qK23 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q20 = ZERO;
	D128 q21 = ZERO;
	D128 q22 = ZERO;
	D128 q23 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    D128 z21 = LOAD(&Iz(a,2,1));
	    q20 = ADD(q20, MUL(MUL(x00, y00), z21));
	    q21 = ADD(q21, MUL(MUL(x10, y10), z01));
	    q22 = ADD(q22, MUL(MUL(x10, y00), z11));
	    q23 = ADD(q23, MUL(MUL(x00, y10), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    D128 z21 = LOAD1(&Iz(a,2,1));
	    q20 = ADD1(q20, MUL1(MUL1(x00, y00), z21));
	    q21 = ADD1(q21, MUL1(MUL1(x10, y10), z01));
	    q22 = ADD1(q22, MUL1(MUL1(x10, y00), z11));
	    q23 = ADD1(q23, MUL1(MUL1(x00, y10), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C11 = LOADDUP(&C[k+1]);
	qK20 = ADD(qK20, MUL(C11, HADD(q20, q21)));
	qK22 = ADD(qK22, MUL(C11, HADD(q22, q23)));

#else // SSE
	    
	// function registers
	T q20 = 0.0;
	T q21 = 0.0;
	T q22 = 0.0;
	T q23 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q20 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,2,1);
	    q21 += Ix(a,1,0)*Iy(a,1,0)*Iz(a,0,1);
	    q22 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,1,1);
	    q23 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK20 += q20*C[k+1];
	//I[20] += q20*C[k+1];
	I[20] += q20*C_[1];
	//qK21 += q21*C[k+1];
	//I[21] += q21*C[k+1];
	I[21] += q21*C_[1];
	//qK22 += q22*C[k+1];
	//I[22] += q22*C[k+1];
	I[22] += q22*C_[1];
	//qK23 += q23*C[k+1];
	//I[23] += q23*C[k+1];
	I[23] += q23*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK20 = MUL(q, qK20);
	qK22 = MUL(q, qK22);
	num += 4; //num += (fabs(I[22]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[20]) & 0xF) {
	    // 20
	    STOREU(&I[20], ADD(qK20, LOADU(&I[20])));
	    STOREU(&I[22], ADD(qK22, LOADU(&I[22])));
	}
	else {
	    STORE(&I[20], ADD(qK20, LOADU(&I[20])));
	    STORE(&I[22], ADD(qK22, LOADU(&I[22])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[20] += scale*NORMALIZE[6]*NORMALIZE[3]*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*NORMALIZE[7]*NORMALIZE[3]*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*NORMALIZE[8]*NORMALIZE[3]*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*NORMALIZE[9]*NORMALIZE[3]*qK23;
	// num += (fabs(I[23]) >= tol);
    }
    else {
	// I[20] += scale*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*qK23;
	// num += (fabs(I[23]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <fsp| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::F,rysq::SP> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::F,rysq::SP> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 4;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x30 = LOAD(&Ix(a,3,0));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y30 = LOAD(&Iy(a,3,0));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z30 = LOAD(&Iz(a,3,0));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x30, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x00, y30), z00));
	    q2 = ADD(q2, MUL(MUL(x00, y00), z30));
	    q3 = ADD(q3, MUL(MUL(x20, y10), z00));
	    q4 = ADD(q4, MUL(MUL(x20, y00), z10));
	    q5 = ADD(q5, MUL(MUL(x10, y20), z00));
	    q6 = ADD(q6, MUL(MUL(x00, y20), z10));
	    q7 = ADD(q7, MUL(MUL(x10, y00), z20));
	    q8 = ADD(q8, MUL(MUL(x00, y10), z20));
	    q9 = ADD(q9, MUL(MUL(x10, y10), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x30 = LOAD1(&Ix(a,3,0));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y30 = LOAD1(&Iy(a,3,0));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z30 = LOAD1(&Iz(a,3,0));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x30, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x00, y30), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x00, y00), z30));
	    q3 = ADD1(q3, MUL1(MUL1(x20, y10), z00));
	    q4 = ADD1(q4, MUL1(MUL1(x20, y00), z10));
	    q5 = ADD1(q5, MUL1(MUL1(x10, y20), z00));
	    q6 = ADD1(q6, MUL1(MUL1(x00, y20), z10));
	    q7 = ADD1(q7, MUL1(MUL1(x10, y00), z20));
	    q8 = ADD1(q8, MUL1(MUL1(x00, y10), z20));
	    q9 = ADD1(q9, MUL1(MUL1(x10, y10), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C00 = LOADDUP(&C[k+0]);
	qK0 = ADD(qK0, MUL(C00, HADD(q0, q1)));
	qK2 = ADD(qK2, MUL(C00, HADD(q2, q3)));
	qK4 = ADD(qK4, MUL(C00, HADD(q4, q5)));
	qK6 = ADD(qK6, MUL(C00, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C00, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,3,0)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,0,0)*Iy(a,3,0)*Iz(a,0,0);
	    q2 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,3,0);
	    q3 += Ix(a,2,0)*Iy(a,1,0)*Iz(a,0,0);
	    q4 += Ix(a,2,0)*Iy(a,0,0)*Iz(a,1,0);
	    q5 += Ix(a,1,0)*Iy(a,2,0)*Iz(a,0,0);
	    q6 += Ix(a,0,0)*Iy(a,2,0)*Iz(a,1,0);
	    q7 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,2,0);
	    q8 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,2,0);
	    q9 += Ix(a,1,0)*Iy(a,1,0)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+0];
	//I[1] += q1*C[k+0];
	I[1] += q1*C_[0];
	//qK2 += q2*C[k+0];
	//I[2] += q2*C[k+0];
	I[2] += q2*C_[0];
	//qK3 += q3*C[k+0];
	//I[3] += q3*C[k+0];
	I[3] += q3*C_[0];
	//qK4 += q4*C[k+0];
	//I[4] += q4*C[k+0];
	I[4] += q4*C_[0];
	//qK5 += q5*C[k+0];
	//I[5] += q5*C[k+0];
	I[5] += q5*C_[0];
	//qK6 += q6*C[k+0];
	//I[6] += q6*C[k+0];
	I[6] += q6*C_[0];
	//qK7 += q7*C[k+0];
	//I[7] += q7*C[k+0];
	I[7] += q7*C_[0];
	//qK8 += q8*C[k+0];
	//I[8] += q8*C[k+0];
	I[8] += q8*C_[0];
	//qK9 += q9*C[k+0];
	//I[9] += q9*C[k+0];
	I[9] += q9*C_[0];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[10]*NORMALIZE[0]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[11]*NORMALIZE[0]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[12]*NORMALIZE[0]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[13]*NORMALIZE[0]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[14]*NORMALIZE[0]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[15]*NORMALIZE[0]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[16]*NORMALIZE[0]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[17]*NORMALIZE[0]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[18]*NORMALIZE[0]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[19]*NORMALIZE[0]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK10 = ZERO;
     D128 qK12 = ZERO;
     D128 qK14 = ZERO;
     D128 qK16 = ZERO;
     D128 qK18 = ZERO;
#else
     //double qK10 = 0.0;
     //double qK11 = 0.0;
     //double qK12 = 0.0;
     //double qK13 = 0.0;
     //double qK14 = 0.0;
     //double qK15 = 0.0;
     //double qK16 = 0.0;
     //double qK17 = 0.0;
     //double qK18 = 0.0;
     //double qK19 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q10 = ZERO;
	D128 q11 = ZERO;
	D128 q12 = ZERO;
	D128 q13 = ZERO;
	D128 q14 = ZERO;
	D128 q15 = ZERO;
	D128 q16 = ZERO;
	D128 q17 = ZERO;
	D128 q18 = ZERO;
	D128 q19 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x31 = LOAD(&Ix(a,3,1));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 x21 = LOAD(&Ix(a,2,1));
	    D128 y30 = LOAD(&Iy(a,3,0));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z30 = LOAD(&Iz(a,3,0));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q10 = ADD(q10, MUL(MUL(x31, y00), z00));
	    q11 = ADD(q11, MUL(MUL(x01, y30), z00));
	    q12 = ADD(q12, MUL(MUL(x01, y00), z30));
	    q13 = ADD(q13, MUL(MUL(x21, y10), z00));
	    q14 = ADD(q14, MUL(MUL(x21, y00), z10));
	    q15 = ADD(q15, MUL(MUL(x11, y20), z00));
	    q16 = ADD(q16, MUL(MUL(x01, y20), z10));
	    q17 = ADD(q17, MUL(MUL(x11, y00), z20));
	    q18 = ADD(q18, MUL(MUL(x01, y10), z20));
	    q19 = ADD(q19, MUL(MUL(x11, y10), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x31 = LOAD1(&Ix(a,3,1));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 x21 = LOAD1(&Ix(a,2,1));
	    D128 y30 = LOAD1(&Iy(a,3,0));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z30 = LOAD1(&Iz(a,3,0));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q10 = ADD1(q10, MUL1(MUL1(x31, y00), z00));
	    q11 = ADD1(q11, MUL1(MUL1(x01, y30), z00));
	    q12 = ADD1(q12, MUL1(MUL1(x01, y00), z30));
	    q13 = ADD1(q13, MUL1(MUL1(x21, y10), z00));
	    q14 = ADD1(q14, MUL1(MUL1(x21, y00), z10));
	    q15 = ADD1(q15, MUL1(MUL1(x11, y20), z00));
	    q16 = ADD1(q16, MUL1(MUL1(x01, y20), z10));
	    q17 = ADD1(q17, MUL1(MUL1(x11, y00), z20));
	    q18 = ADD1(q18, MUL1(MUL1(x01, y10), z20));
	    q19 = ADD1(q19, MUL1(MUL1(x11, y10), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C11 = LOADDUP(&C[k+1]);
	qK10 = ADD(qK10, MUL(C11, HADD(q10, q11)));
	qK12 = ADD(qK12, MUL(C11, HADD(q12, q13)));
	qK14 = ADD(qK14, MUL(C11, HADD(q14, q15)));
	qK16 = ADD(qK16, MUL(C11, HADD(q16, q17)));
	qK18 = ADD(qK18, MUL(C11, HADD(q18, q19)));

#else // SSE
	    
	// function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
	T q16 = 0.0;
	T q17 = 0.0;
	T q18 = 0.0;
	T q19 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += Ix(a,3,1)*Iy(a,0,0)*Iz(a,0,0);
	    q11 += Ix(a,0,1)*Iy(a,3,0)*Iz(a,0,0);
	    q12 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,3,0);
	    q13 += Ix(a,2,1)*Iy(a,1,0)*Iz(a,0,0);
	    q14 += Ix(a,2,1)*Iy(a,0,0)*Iz(a,1,0);
	    q15 += Ix(a,1,1)*Iy(a,2,0)*Iz(a,0,0);
	    q16 += Ix(a,0,1)*Iy(a,2,0)*Iz(a,1,0);
	    q17 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,2,0);
	    q18 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,2,0);
	    q19 += Ix(a,1,1)*Iy(a,1,0)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK10 += q10*C[k+1];
	//I[10] += q10*C[k+1];
	I[10] += q10*C_[1];
	//qK11 += q11*C[k+1];
	//I[11] += q11*C[k+1];
	I[11] += q11*C_[1];
	//qK12 += q12*C[k+1];
	//I[12] += q12*C[k+1];
	I[12] += q12*C_[1];
	//qK13 += q13*C[k+1];
	//I[13] += q13*C[k+1];
	I[13] += q13*C_[1];
	//qK14 += q14*C[k+1];
	//I[14] += q14*C[k+1];
	I[14] += q14*C_[1];
	//qK15 += q15*C[k+1];
	//I[15] += q15*C[k+1];
	I[15] += q15*C_[1];
	//qK16 += q16*C[k+1];
	//I[16] += q16*C[k+1];
	I[16] += q16*C_[1];
	//qK17 += q17*C[k+1];
	//I[17] += q17*C[k+1];
	I[17] += q17*C_[1];
	//qK18 += q18*C[k+1];
	//I[18] += q18*C[k+1];
	I[18] += q18*C_[1];
	//qK19 += q19*C[k+1];
	//I[19] += q19*C[k+1];
	I[19] += q19*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK10 = MUL(q, qK10);
	qK12 = MUL(q, qK12);
	qK14 = MUL(q, qK14);
	qK16 = MUL(q, qK16);
	qK18 = MUL(q, qK18);
	num += 10; //num += (fabs(I[18]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[10]) & 0xF) {
	    // 10
	    STOREU(&I[10], ADD(qK10, LOADU(&I[10])));
	    STOREU(&I[12], ADD(qK12, LOADU(&I[12])));
	    STOREU(&I[14], ADD(qK14, LOADU(&I[14])));
	    STOREU(&I[16], ADD(qK16, LOADU(&I[16])));
	    STOREU(&I[18], ADD(qK18, LOADU(&I[18])));
	}
	else {
	    STORE(&I[10], ADD(qK10, LOADU(&I[10])));
	    STORE(&I[12], ADD(qK12, LOADU(&I[12])));
	    STORE(&I[14], ADD(qK14, LOADU(&I[14])));
	    STORE(&I[16], ADD(qK16, LOADU(&I[16])));
	    STORE(&I[18], ADD(qK18, LOADU(&I[18])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[10] += scale*NORMALIZE[10]*NORMALIZE[1]*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*NORMALIZE[11]*NORMALIZE[1]*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*NORMALIZE[12]*NORMALIZE[1]*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*NORMALIZE[13]*NORMALIZE[1]*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*NORMALIZE[14]*NORMALIZE[1]*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*NORMALIZE[15]*NORMALIZE[1]*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*NORMALIZE[16]*NORMALIZE[1]*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*NORMALIZE[17]*NORMALIZE[1]*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*NORMALIZE[18]*NORMALIZE[1]*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*NORMALIZE[19]*NORMALIZE[1]*qK19;
	// num += (fabs(I[19]) >= tol);
    }
    else {
	// I[10] += scale*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*qK15;
	// num += (fabs(I[15]) >= tol);
	// I[16] += scale*qK16;
	// num += (fabs(I[16]) >= tol);
	// I[17] += scale*qK17;
	// num += (fabs(I[17]) >= tol);
	// I[18] += scale*qK18;
	// num += (fabs(I[18]) >= tol);
	// I[19] += scale*qK19;
	// num += (fabs(I[19]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK20 = ZERO;
     D128 qK22 = ZERO;
     D128 qK24 = ZERO;
     D128 qK26 = ZERO;
     D128 qK28 = ZERO;
#else
     //double qK20 = 0.0;
     //double qK21 = 0.0;
     //double qK22 = 0.0;
     //double qK23 = 0.0;
     //double qK24 = 0.0;
     //double qK25 = 0.0;
     //double qK26 = 0.0;
     //double qK27 = 0.0;
     //double qK28 = 0.0;
     //double qK29 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q20 = ZERO;
	D128 q21 = ZERO;
	D128 q22 = ZERO;
	D128 q23 = ZERO;
	D128 q24 = ZERO;
	D128 q25 = ZERO;
	D128 q26 = ZERO;
	D128 q27 = ZERO;
	D128 q28 = ZERO;
	D128 q29 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x30 = LOAD(&Ix(a,3,0));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y31 = LOAD(&Iy(a,3,1));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 y21 = LOAD(&Iy(a,2,1));
	    D128 z30 = LOAD(&Iz(a,3,0));
	    D128 z20 = LOAD(&Iz(a,2,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q20 = ADD(q20, MUL(MUL(x30, y01), z00));
	    q21 = ADD(q21, MUL(MUL(x00, y31), z00));
	    q22 = ADD(q22, MUL(MUL(x00, y01), z30));
	    q23 = ADD(q23, MUL(MUL(x20, y11), z00));
	    q24 = ADD(q24, MUL(MUL(x20, y01), z10));
	    q25 = ADD(q25, MUL(MUL(x10, y21), z00));
	    q26 = ADD(q26, MUL(MUL(x00, y21), z10));
	    q27 = ADD(q27, MUL(MUL(x10, y01), z20));
	    q28 = ADD(q28, MUL(MUL(x00, y11), z20));
	    q29 = ADD(q29, MUL(MUL(x10, y11), z10));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x30 = LOAD1(&Ix(a,3,0));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y31 = LOAD1(&Iy(a,3,1));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 y21 = LOAD1(&Iy(a,2,1));
	    D128 z30 = LOAD1(&Iz(a,3,0));
	    D128 z20 = LOAD1(&Iz(a,2,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q20 = ADD1(q20, MUL1(MUL1(x30, y01), z00));
	    q21 = ADD1(q21, MUL1(MUL1(x00, y31), z00));
	    q22 = ADD1(q22, MUL1(MUL1(x00, y01), z30));
	    q23 = ADD1(q23, MUL1(MUL1(x20, y11), z00));
	    q24 = ADD1(q24, MUL1(MUL1(x20, y01), z10));
	    q25 = ADD1(q25, MUL1(MUL1(x10, y21), z00));
	    q26 = ADD1(q26, MUL1(MUL1(x00, y21), z10));
	    q27 = ADD1(q27, MUL1(MUL1(x10, y01), z20));
	    q28 = ADD1(q28, MUL1(MUL1(x00, y11), z20));
	    q29 = ADD1(q29, MUL1(MUL1(x10, y11), z10));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C11 = LOADDUP(&C[k+1]);
	qK20 = ADD(qK20, MUL(C11, HADD(q20, q21)));
	qK22 = ADD(qK22, MUL(C11, HADD(q22, q23)));
	qK24 = ADD(qK24, MUL(C11, HADD(q24, q25)));
	qK26 = ADD(qK26, MUL(C11, HADD(q26, q27)));
	qK28 = ADD(qK28, MUL(C11, HADD(q28, q29)));

#else // SSE
	    
	// function registers
	T q20 = 0.0;
	T q21 = 0.0;
	T q22 = 0.0;
	T q23 = 0.0;
	T q24 = 0.0;
	T q25 = 0.0;
	T q26 = 0.0;
	T q27 = 0.0;
	T q28 = 0.0;
	T q29 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q20 += Ix(a,3,0)*Iy(a,0,1)*Iz(a,0,0);
	    q21 += Ix(a,0,0)*Iy(a,3,1)*Iz(a,0,0);
	    q22 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,3,0);
	    q23 += Ix(a,2,0)*Iy(a,1,1)*Iz(a,0,0);
	    q24 += Ix(a,2,0)*Iy(a,0,1)*Iz(a,1,0);
	    q25 += Ix(a,1,0)*Iy(a,2,1)*Iz(a,0,0);
	    q26 += Ix(a,0,0)*Iy(a,2,1)*Iz(a,1,0);
	    q27 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,2,0);
	    q28 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,2,0);
	    q29 += Ix(a,1,0)*Iy(a,1,1)*Iz(a,1,0);
	}
	    
	//contraction coefficients
	//qK20 += q20*C[k+1];
	//I[20] += q20*C[k+1];
	I[20] += q20*C_[1];
	//qK21 += q21*C[k+1];
	//I[21] += q21*C[k+1];
	I[21] += q21*C_[1];
	//qK22 += q22*C[k+1];
	//I[22] += q22*C[k+1];
	I[22] += q22*C_[1];
	//qK23 += q23*C[k+1];
	//I[23] += q23*C[k+1];
	I[23] += q23*C_[1];
	//qK24 += q24*C[k+1];
	//I[24] += q24*C[k+1];
	I[24] += q24*C_[1];
	//qK25 += q25*C[k+1];
	//I[25] += q25*C[k+1];
	I[25] += q25*C_[1];
	//qK26 += q26*C[k+1];
	//I[26] += q26*C[k+1];
	I[26] += q26*C_[1];
	//qK27 += q27*C[k+1];
	//I[27] += q27*C[k+1];
	I[27] += q27*C_[1];
	//qK28 += q28*C[k+1];
	//I[28] += q28*C[k+1];
	I[28] += q28*C_[1];
	//qK29 += q29*C[k+1];
	//I[29] += q29*C[k+1];
	I[29] += q29*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK20 = MUL(q, qK20);
	qK22 = MUL(q, qK22);
	qK24 = MUL(q, qK24);
	qK26 = MUL(q, qK26);
	qK28 = MUL(q, qK28);
	num += 10; //num += (fabs(I[28]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[20]) & 0xF) {
	    // 20
	    STOREU(&I[20], ADD(qK20, LOADU(&I[20])));
	    STOREU(&I[22], ADD(qK22, LOADU(&I[22])));
	    STOREU(&I[24], ADD(qK24, LOADU(&I[24])));
	    STOREU(&I[26], ADD(qK26, LOADU(&I[26])));
	    STOREU(&I[28], ADD(qK28, LOADU(&I[28])));
	}
	else {
	    STORE(&I[20], ADD(qK20, LOADU(&I[20])));
	    STORE(&I[22], ADD(qK22, LOADU(&I[22])));
	    STORE(&I[24], ADD(qK24, LOADU(&I[24])));
	    STORE(&I[26], ADD(qK26, LOADU(&I[26])));
	    STORE(&I[28], ADD(qK28, LOADU(&I[28])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[20] += scale*NORMALIZE[10]*NORMALIZE[2]*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*NORMALIZE[11]*NORMALIZE[2]*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*NORMALIZE[12]*NORMALIZE[2]*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*NORMALIZE[13]*NORMALIZE[2]*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*NORMALIZE[14]*NORMALIZE[2]*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*NORMALIZE[15]*NORMALIZE[2]*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*NORMALIZE[16]*NORMALIZE[2]*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*NORMALIZE[17]*NORMALIZE[2]*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*NORMALIZE[18]*NORMALIZE[2]*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*NORMALIZE[19]*NORMALIZE[2]*qK29;
	// num += (fabs(I[29]) >= tol);
    }
    else {
	// I[20] += scale*qK20;
	// num += (fabs(I[20]) >= tol);
	// I[21] += scale*qK21;
	// num += (fabs(I[21]) >= tol);
	// I[22] += scale*qK22;
	// num += (fabs(I[22]) >= tol);
	// I[23] += scale*qK23;
	// num += (fabs(I[23]) >= tol);
	// I[24] += scale*qK24;
	// num += (fabs(I[24]) >= tol);
	// I[25] += scale*qK25;
	// num += (fabs(I[25]) >= tol);
	// I[26] += scale*qK26;
	// num += (fabs(I[26]) >= tol);
	// I[27] += scale*qK27;
	// num += (fabs(I[27]) >= tol);
	// I[28] += scale*qK28;
	// num += (fabs(I[28]) >= tol);
	// I[29] += scale*qK29;
	// num += (fabs(I[29]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK30 = ZERO;
     D128 qK32 = ZERO;
     D128 qK34 = ZERO;
     D128 qK36 = ZERO;
     D128 qK38 = ZERO;
#else
     //double qK30 = 0.0;
     //double qK31 = 0.0;
     //double qK32 = 0.0;
     //double qK33 = 0.0;
     //double qK34 = 0.0;
     //double qK35 = 0.0;
     //double qK36 = 0.0;
     //double qK37 = 0.0;
     //double qK38 = 0.0;
     //double qK39 = 0.0;
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q30 = ZERO;
	D128 q31 = ZERO;
	D128 q32 = ZERO;
	D128 q33 = ZERO;
	D128 q34 = ZERO;
	D128 q35 = ZERO;
	D128 q36 = ZERO;
	D128 q37 = ZERO;
	D128 q38 = ZERO;
	D128 q39 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x30 = LOAD(&Ix(a,3,0));
	    D128 x20 = LOAD(&Ix(a,2,0));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y30 = LOAD(&Iy(a,3,0));
	    D128 y20 = LOAD(&Iy(a,2,0));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z31 = LOAD(&Iz(a,3,1));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    D128 z21 = LOAD(&Iz(a,2,1));
	    q30 = ADD(q30, MUL(MUL(x30, y00), z01));
	    q31 = ADD(q31, MUL(MUL(x00, y30), z01));
	    q32 = ADD(q32, MUL(MUL(x00, y00), z31));
	    q33 = ADD(q33, MUL(MUL(x20, y10), z01));
	    q34 = ADD(q34, MUL(MUL(x20, y00), z11));
	    q35 = ADD(q35, MUL(MUL(x10, y20), z01));
	    q36 = ADD(q36, MUL(MUL(x00, y20), z11));
	    q37 = ADD(q37, MUL(MUL(x10, y00), z21));
	    q38 = ADD(q38, MUL(MUL(x00, y10), z21));
	    q39 = ADD(q39, MUL(MUL(x10, y10), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x30 = LOAD1(&Ix(a,3,0));
	    D128 x20 = LOAD1(&Ix(a,2,0));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y30 = LOAD1(&Iy(a,3,0));
	    D128 y20 = LOAD1(&Iy(a,2,0));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z31 = LOAD1(&Iz(a,3,1));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    D128 z21 = LOAD1(&Iz(a,2,1));
	    q30 = ADD1(q30, MUL1(MUL1(x30, y00), z01));
	    q31 = ADD1(q31, MUL1(MUL1(x00, y30), z01));
	    q32 = ADD1(q32, MUL1(MUL1(x00, y00), z31));
	    q33 = ADD1(q33, MUL1(MUL1(x20, y10), z01));
	    q34 = ADD1(q34, MUL1(MUL1(x20, y00), z11));
	    q35 = ADD1(q35, MUL1(MUL1(x10, y20), z01));
	    q36 = ADD1(q36, MUL1(MUL1(x00, y20), z11));
	    q37 = ADD1(q37, MUL1(MUL1(x10, y00), z21));
	    q38 = ADD1(q38, MUL1(MUL1(x00, y10), z21));
	    q39 = ADD1(q39, MUL1(MUL1(x10, y10), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C11 = LOADDUP(&C[k+1]);
	qK30 = ADD(qK30, MUL(C11, HADD(q30, q31)));
	qK32 = ADD(qK32, MUL(C11, HADD(q32, q33)));
	qK34 = ADD(qK34, MUL(C11, HADD(q34, q35)));
	qK36 = ADD(qK36, MUL(C11, HADD(q36, q37)));
	qK38 = ADD(qK38, MUL(C11, HADD(q38, q39)));

#else // SSE
	    
	// function registers
	T q30 = 0.0;
	T q31 = 0.0;
	T q32 = 0.0;
	T q33 = 0.0;
	T q34 = 0.0;
	T q35 = 0.0;
	T q36 = 0.0;
	T q37 = 0.0;
	T q38 = 0.0;
	T q39 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q30 += Ix(a,3,0)*Iy(a,0,0)*Iz(a,0,1);
	    q31 += Ix(a,0,0)*Iy(a,3,0)*Iz(a,0,1);
	    q32 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,3,1);
	    q33 += Ix(a,2,0)*Iy(a,1,0)*Iz(a,0,1);
	    q34 += Ix(a,2,0)*Iy(a,0,0)*Iz(a,1,1);
	    q35 += Ix(a,1,0)*Iy(a,2,0)*Iz(a,0,1);
	    q36 += Ix(a,0,0)*Iy(a,2,0)*Iz(a,1,1);
	    q37 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,2,1);
	    q38 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,2,1);
	    q39 += Ix(a,1,0)*Iy(a,1,0)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK30 += q30*C[k+1];
	//I[30] += q30*C[k+1];
	I[30] += q30*C_[1];
	//qK31 += q31*C[k+1];
	//I[31] += q31*C[k+1];
	I[31] += q31*C_[1];
	//qK32 += q32*C[k+1];
	//I[32] += q32*C[k+1];
	I[32] += q32*C_[1];
	//qK33 += q33*C[k+1];
	//I[33] += q33*C[k+1];
	I[33] += q33*C_[1];
	//qK34 += q34*C[k+1];
	//I[34] += q34*C[k+1];
	I[34] += q34*C_[1];
	//qK35 += q35*C[k+1];
	//I[35] += q35*C[k+1];
	I[35] += q35*C_[1];
	//qK36 += q36*C[k+1];
	//I[36] += q36*C[k+1];
	I[36] += q36*C_[1];
	//qK37 += q37*C[k+1];
	//I[37] += q37*C[k+1];
	I[37] += q37*C_[1];
	//qK38 += q38*C[k+1];
	//I[38] += q38*C[k+1];
	I[38] += q38*C_[1];
	//qK39 += q39*C[k+1];
	//I[39] += q39*C[k+1];
	I[39] += q39*C_[1];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK30 = MUL(q, qK30);
	qK32 = MUL(q, qK32);
	qK34 = MUL(q, qK34);
	qK36 = MUL(q, qK36);
	qK38 = MUL(q, qK38);
	num += 10; //num += (fabs(I[38]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[30]) & 0xF) {
	    // 30
	    STOREU(&I[30], ADD(qK30, LOADU(&I[30])));
	    STOREU(&I[32], ADD(qK32, LOADU(&I[32])));
	    STOREU(&I[34], ADD(qK34, LOADU(&I[34])));
	    STOREU(&I[36], ADD(qK36, LOADU(&I[36])));
	    STOREU(&I[38], ADD(qK38, LOADU(&I[38])));
	}
	else {
	    STORE(&I[30], ADD(qK30, LOADU(&I[30])));
	    STORE(&I[32], ADD(qK32, LOADU(&I[32])));
	    STORE(&I[34], ADD(qK34, LOADU(&I[34])));
	    STORE(&I[36], ADD(qK36, LOADU(&I[36])));
	    STORE(&I[38], ADD(qK38, LOADU(&I[38])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[30] += scale*NORMALIZE[10]*NORMALIZE[3]*qK30;
	// num += (fabs(I[30]) >= tol);
	// I[31] += scale*NORMALIZE[11]*NORMALIZE[3]*qK31;
	// num += (fabs(I[31]) >= tol);
	// I[32] += scale*NORMALIZE[12]*NORMALIZE[3]*qK32;
	// num += (fabs(I[32]) >= tol);
	// I[33] += scale*NORMALIZE[13]*NORMALIZE[3]*qK33;
	// num += (fabs(I[33]) >= tol);
	// I[34] += scale*NORMALIZE[14]*NORMALIZE[3]*qK34;
	// num += (fabs(I[34]) >= tol);
	// I[35] += scale*NORMALIZE[15]*NORMALIZE[3]*qK35;
	// num += (fabs(I[35]) >= tol);
	// I[36] += scale*NORMALIZE[16]*NORMALIZE[3]*qK36;
	// num += (fabs(I[36]) >= tol);
	// I[37] += scale*NORMALIZE[17]*NORMALIZE[3]*qK37;
	// num += (fabs(I[37]) >= tol);
	// I[38] += scale*NORMALIZE[18]*NORMALIZE[3]*qK38;
	// num += (fabs(I[38]) >= tol);
	// I[39] += scale*NORMALIZE[19]*NORMALIZE[3]*qK39;
	// num += (fabs(I[39]) >= tol);
    }
    else {
	// I[30] += scale*qK30;
	// num += (fabs(I[30]) >= tol);
	// I[31] += scale*qK31;
	// num += (fabs(I[31]) >= tol);
	// I[32] += scale*qK32;
	// num += (fabs(I[32]) >= tol);
	// I[33] += scale*qK33;
	// num += (fabs(I[33]) >= tol);
	// I[34] += scale*qK34;
	// num += (fabs(I[34]) >= tol);
	// I[35] += scale*qK35;
	// num += (fabs(I[35]) >= tol);
	// I[36] += scale*qK36;
	// num += (fabs(I[36]) >= tol);
	// I[37] += scale*qK37;
	// num += (fabs(I[37]) >= tol);
	// I[38] += scale*qK38;
	// num += (fabs(I[38]) >= tol);
	// I[39] += scale*qK39;
	// num += (fabs(I[39]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}



/** 
    @brief <spsp| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param[out] I integral batch
    @return number of screened integrals
*/


template<>
struct impl< meta::state<rysq::SP,rysq::SP> > {

    template<typename T>
    struct aligned {
#if ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3))  
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::state<rysq::SP,rysq::SP> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I) {

    const int Li1 = 2;
    
    int num = 0;
    
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK0 = ZERO;
     D128 qK2 = ZERO;
     D128 qK4 = ZERO;
     D128 qK6 = ZERO;
     D128 qK8 = ZERO;
#else
     //double qK0 = 0.0;
     //double qK1 = 0.0;
     //double qK2 = 0.0;
     //double qK3 = 0.0;
     //double qK4 = 0.0;
     //double qK5 = 0.0;
     //double qK6 = 0.0;
     //double qK7 = 0.0;
     //double qK8 = 0.0;
     //double qK9 = 0.0;
#endif

    for(int k = 0; k < K*4; k += 4) {
	double C_[4];
	for (int i = 0; i < 4; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q0 = ZERO;
	D128 q1 = ZERO;
	D128 q2 = ZERO;
	D128 q3 = ZERO;
	D128 q4 = ZERO;
	D128 q5 = ZERO;
	D128 q6 = ZERO;
	D128 q7 = ZERO;
	D128 q8 = ZERO;
	D128 q9 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x01 = LOAD(&Ix(a,0,1));
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 x11 = LOAD(&Ix(a,1,1));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    q0 = ADD(q0, MUL(MUL(x00, y00), z00));
	    q1 = ADD(q1, MUL(MUL(x10, y00), z00));
	    q2 = ADD(q2, MUL(MUL(x00, y10), z00));
	    q3 = ADD(q3, MUL(MUL(x00, y00), z10));
	    q4 = ADD(q4, MUL(MUL(x01, y00), z00));
	    q5 = ADD(q5, MUL(MUL(x11, y00), z00));
	    q6 = ADD(q6, MUL(MUL(x01, y10), z00));
	    q7 = ADD(q7, MUL(MUL(x01, y00), z10));
	    q8 = ADD(q8, MUL(MUL(x00, y01), z00));
	    q9 = ADD(q9, MUL(MUL(x10, y01), z00));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x01 = LOAD1(&Ix(a,0,1));
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 x11 = LOAD1(&Ix(a,1,1));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    q0 = ADD1(q0, MUL1(MUL1(x00, y00), z00));
	    q1 = ADD1(q1, MUL1(MUL1(x10, y00), z00));
	    q2 = ADD1(q2, MUL1(MUL1(x00, y10), z00));
	    q3 = ADD1(q3, MUL1(MUL1(x00, y00), z10));
	    q4 = ADD1(q4, MUL1(MUL1(x01, y00), z00));
	    q5 = ADD1(q5, MUL1(MUL1(x11, y00), z00));
	    q6 = ADD1(q6, MUL1(MUL1(x01, y10), z00));
	    q7 = ADD1(q7, MUL1(MUL1(x01, y00), z10));
	    q8 = ADD1(q8, MUL1(MUL1(x00, y01), z00));
	    q9 = ADD1(q9, MUL1(MUL1(x10, y01), z00));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C01 = LOADU(&C[k+0]);
	qK0 = ADD(qK0, MUL(C01, HADD(q0, q1)));
	D128 C11 = LOADDUP(&C[k+1]);
	qK2 = ADD(qK2, MUL(C11, HADD(q2, q3)));
	D128 C23 = LOADU(&C[k+2]);
	qK4 = ADD(qK4, MUL(C23, HADD(q4, q5)));
	D128 C33 = LOADDUP(&C[k+3]);
	qK6 = ADD(qK6, MUL(C33, HADD(q6, q7)));
	qK8 = ADD(qK8, MUL(C23, HADD(q8, q9)));

#else // SSE
	    
	// function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,0,0);
	    q1 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,0,0);
	    q2 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,0,0);
	    q3 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,1,0);
	    q4 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,0,0);
	    q5 += Ix(a,1,1)*Iy(a,0,0)*Iz(a,0,0);
	    q6 += Ix(a,0,1)*Iy(a,1,0)*Iz(a,0,0);
	    q7 += Ix(a,0,1)*Iy(a,0,0)*Iz(a,1,0);
	    q8 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,0,0);
	    q9 += Ix(a,1,0)*Iy(a,0,1)*Iz(a,0,0);
	}
	    
	//contraction coefficients
	//qK0 += q0*C[k+0];
	//I[0] += q0*C[k+0];
	I[0] += q0*C_[0];
	//qK1 += q1*C[k+1];
	//I[1] += q1*C[k+1];
	I[1] += q1*C_[1];
	//qK2 += q2*C[k+1];
	//I[2] += q2*C[k+1];
	I[2] += q2*C_[1];
	//qK3 += q3*C[k+1];
	//I[3] += q3*C[k+1];
	I[3] += q3*C_[1];
	//qK4 += q4*C[k+2];
	//I[4] += q4*C[k+2];
	I[4] += q4*C_[2];
	//qK5 += q5*C[k+3];
	//I[5] += q5*C[k+3];
	I[5] += q5*C_[3];
	//qK6 += q6*C[k+3];
	//I[6] += q6*C[k+3];
	I[6] += q6*C_[3];
	//qK7 += q7*C[k+3];
	//I[7] += q7*C[k+3];
	I[7] += q7*C_[3];
	//qK8 += q8*C[k+2];
	//I[8] += q8*C[k+2];
	I[8] += q8*C_[2];
	//qK9 += q9*C[k+3];
	//I[9] += q9*C[k+3];
	I[9] += q9*C_[3];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK0 = MUL(q, qK0);
	qK2 = MUL(q, qK2);
	qK4 = MUL(q, qK4);
	qK6 = MUL(q, qK6);
	qK8 = MUL(q, qK8);
	num += 10; //num += (fabs(I[8]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[0]) & 0xF) {
	    // 0
	    STOREU(&I[0], ADD(qK0, LOADU(&I[0])));
	    STOREU(&I[2], ADD(qK2, LOADU(&I[2])));
	    STOREU(&I[4], ADD(qK4, LOADU(&I[4])));
	    STOREU(&I[6], ADD(qK6, LOADU(&I[6])));
	    STOREU(&I[8], ADD(qK8, LOADU(&I[8])));
	}
	else {
	    STORE(&I[0], ADD(qK0, LOADU(&I[0])));
	    STORE(&I[2], ADD(qK2, LOADU(&I[2])));
	    STORE(&I[4], ADD(qK4, LOADU(&I[4])));
	    STORE(&I[6], ADD(qK6, LOADU(&I[6])));
	    STORE(&I[8], ADD(qK8, LOADU(&I[8])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[0] += scale*NORMALIZE[0]*NORMALIZE[0]*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*NORMALIZE[1]*NORMALIZE[0]*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*NORMALIZE[2]*NORMALIZE[0]*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*NORMALIZE[3]*NORMALIZE[0]*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*NORMALIZE[0]*NORMALIZE[1]*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*NORMALIZE[1]*NORMALIZE[1]*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*NORMALIZE[2]*NORMALIZE[1]*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*NORMALIZE[3]*NORMALIZE[1]*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*NORMALIZE[0]*NORMALIZE[2]*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*NORMALIZE[1]*NORMALIZE[2]*qK9;
	// num += (fabs(I[9]) >= tol);
    }
    else {
	// I[0] += scale*qK0;
	// num += (fabs(I[0]) >= tol);
	// I[1] += scale*qK1;
	// num += (fabs(I[1]) >= tol);
	// I[2] += scale*qK2;
	// num += (fabs(I[2]) >= tol);
	// I[3] += scale*qK3;
	// num += (fabs(I[3]) >= tol);
	// I[4] += scale*qK4;
	// num += (fabs(I[4]) >= tol);
	// I[5] += scale*qK5;
	// num += (fabs(I[5]) >= tol);
	// I[6] += scale*qK6;
	// num += (fabs(I[6]) >= tol);
	// I[7] += scale*qK7;
	// num += (fabs(I[7]) >= tol);
	// I[8] += scale*qK8;
	// num += (fabs(I[8]) >= tol);
	// I[9] += scale*qK9;
	// num += (fabs(I[9]) >= tol);
    }

#endif // !SSE
    
    
    //contracted registers
#ifdef RYSQ_WITH_SSE
    // std::cout <<  "SSE" << std::endl;
     D128 qK10 = ZERO;
     D128 qK12 = ZERO;
     D128 qK14 = ZERO;
#else
     //double qK10 = 0.0;
     //double qK11 = 0.0;
     //double qK12 = 0.0;
     //double qK13 = 0.0;
     //double qK14 = 0.0;
     //double qK15 = 0.0;
#endif

    for(int k = 0; k < K*4; k += 4) {
	double C_[4];
	for (int i = 0; i < 4; ++ i) C_[i] = C[k + i];

	// // prefetch next contraction
	// _mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	// _mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

#ifdef RYSQ_WITH_SSE
	
	// prefetch next contraction
	_mm_prefetch(Ix + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iy + 3*dim2d, _MM_HINT_T0);
	_mm_prefetch(Iz + 3*dim2d, _MM_HINT_T0);

	D128 q10 = ZERO;
	D128 q11 = ZERO;
	D128 q12 = ZERO;
	D128 q13 = ZERO;
	D128 q14 = ZERO;
	D128 q15 = ZERO;
 	for (int a = 0; a < (int(N)/2)*2; a += 2) {
	    D128 x10 = LOAD(&Ix(a,1,0));
	    D128 x00 = LOAD(&Ix(a,0,0));
	    D128 y01 = LOAD(&Iy(a,0,1));
	    D128 y10 = LOAD(&Iy(a,1,0));
	    D128 y00 = LOAD(&Iy(a,0,0));
	    D128 y11 = LOAD(&Iy(a,1,1));
	    D128 z01 = LOAD(&Iz(a,0,1));
	    D128 z10 = LOAD(&Iz(a,1,0));
	    D128 z00 = LOAD(&Iz(a,0,0));
	    D128 z11 = LOAD(&Iz(a,1,1));
	    q10 = ADD(q10, MUL(MUL(x00, y11), z00));
	    q11 = ADD(q11, MUL(MUL(x00, y01), z10));
	    q12 = ADD(q12, MUL(MUL(x00, y00), z01));
	    q13 = ADD(q13, MUL(MUL(x10, y00), z01));
	    q14 = ADD(q14, MUL(MUL(x00, y10), z01));
	    q15 = ADD(q15, MUL(MUL(x00, y00), z11));
	}
 	for (int a = (int(N)/2)*2; a < int(N); a += 1) {
	    D128 x10 = LOAD1(&Ix(a,1,0));
	    D128 x00 = LOAD1(&Ix(a,0,0));
	    D128 y01 = LOAD1(&Iy(a,0,1));
	    D128 y10 = LOAD1(&Iy(a,1,0));
	    D128 y00 = LOAD1(&Iy(a,0,0));
	    D128 y11 = LOAD1(&Iy(a,1,1));
	    D128 z01 = LOAD1(&Iz(a,0,1));
	    D128 z10 = LOAD1(&Iz(a,1,0));
	    D128 z00 = LOAD1(&Iz(a,0,0));
	    D128 z11 = LOAD1(&Iz(a,1,1));
	    q10 = ADD1(q10, MUL1(MUL1(x00, y11), z00));
	    q11 = ADD1(q11, MUL1(MUL1(x00, y01), z10));
	    q12 = ADD1(q12, MUL1(MUL1(x00, y00), z01));
	    q13 = ADD1(q13, MUL1(MUL1(x10, y00), z01));
	    q14 = ADD1(q14, MUL1(MUL1(x00, y10), z01));
	    q15 = ADD1(q15, MUL1(MUL1(x00, y00), z11));
	}

	//contraction coefficients
	//D128 C0 = _mm_loaddup_pd(&C[k]);
	D128 C33 = LOADDUP(&C[k+3]);
	qK10 = ADD(qK10, MUL(C33, HADD(q10, q11)));
	D128 C23 = LOADU(&C[k+2]);
	qK12 = ADD(qK12, MUL(C23, HADD(q12, q13)));
	qK14 = ADD(qK14, MUL(C33, HADD(q14, q15)));

#else // SSE
	    
	// function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += Ix(a,0,0)*Iy(a,1,1)*Iz(a,0,0);
	    q11 += Ix(a,0,0)*Iy(a,0,1)*Iz(a,1,0);
	    q12 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,0,1);
	    q13 += Ix(a,1,0)*Iy(a,0,0)*Iz(a,0,1);
	    q14 += Ix(a,0,0)*Iy(a,1,0)*Iz(a,0,1);
	    q15 += Ix(a,0,0)*Iy(a,0,0)*Iz(a,1,1);
	}
	    
	//contraction coefficients
	//qK10 += q10*C[k+3];
	//I[10] += q10*C[k+3];
	I[10] += q10*C_[3];
	//qK11 += q11*C[k+3];
	//I[11] += q11*C[k+3];
	I[11] += q11*C_[3];
	//qK12 += q12*C[k+2];
	//I[12] += q12*C[k+2];
	I[12] += q12*C_[2];
	//qK13 += q13*C[k+3];
	//I[13] += q13*C[k+3];
	I[13] += q13*C_[3];
	//qK14 += q14*C[k+3];
	//I[14] += q14*C[k+3];
	I[14] += q14*C_[3];
	//qK15 += q15*C[k+3];
	//I[15] += q15*C[k+3];
	I[15] += q15*C_[3];

#endif // not SSE
	    
	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    // normalization, scaling, and storage

#ifdef RYSQ_WITH_SSE

    {
	D128 q = SET1(scale);
	qK10 = MUL(q, qK10);
	qK12 = MUL(q, qK12);
	qK14 = MUL(q, qK14);
	num += 6; //num += (fabs(I[14]) >= tol);
    }

    if (normalize) {
	throw std::runtime_error("not implemented");
    }
    else {
	if (size_t(&I[10]) & 0xF) {
	    // 10
	    STOREU(&I[10], ADD(qK10, LOADU(&I[10])));
	    STOREU(&I[12], ADD(qK12, LOADU(&I[12])));
	    STOREU(&I[14], ADD(qK14, LOADU(&I[14])));
	}
	else {
	    STORE(&I[10], ADD(qK10, LOADU(&I[10])));
	    STORE(&I[12], ADD(qK12, LOADU(&I[12])));
	    STORE(&I[14], ADD(qK14, LOADU(&I[14])));
	}
    }

#else // !SSE

    if(normalize) {
	// I[10] += scale*NORMALIZE[2]*NORMALIZE[2]*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*NORMALIZE[3]*NORMALIZE[2]*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*NORMALIZE[0]*NORMALIZE[3]*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*NORMALIZE[1]*NORMALIZE[3]*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*NORMALIZE[2]*NORMALIZE[3]*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*NORMALIZE[3]*NORMALIZE[3]*qK15;
	// num += (fabs(I[15]) >= tol);
    }
    else {
	// I[10] += scale*qK10;
	// num += (fabs(I[10]) >= tol);
	// I[11] += scale*qK11;
	// num += (fabs(I[11]) >= tol);
	// I[12] += scale*qK12;
	// num += (fabs(I[12]) >= tol);
	// I[13] += scale*qK13;
	// num += (fabs(I[13]) >= tol);
	// I[14] += scale*qK14;
	// num += (fabs(I[14]) >= tol);
	// I[15] += scale*qK15;
	// num += (fabs(I[15]) >= tol);
    }

#endif // !SSE
    
    
    return num;
}





#undef Ix
#undef Iy
#undef Iz

END_NAMESPACE(rysq, kernel, quadrature)

#undef D128
#undef ZERO
#undef SET1

#undef LOAD
#undef LOADU
#undef LOAD1
#undef LOADDUP

#undef STORE
#undef STOREU
#undef STORE1

#undef MUL
#undef ADD
#undef HADD

#undef MUL1
#undef ADD1

#endif // _RYSQ_KERNEL_QUADRATURE1_IMPL_HPP_

