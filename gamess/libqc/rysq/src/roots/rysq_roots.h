#ifndef RYSQ_ROOTS_H
#define RYSQ_ROOTS_H

#ifndef __device__
#define __device__
#endif

#ifndef __constant__
#define __constant__
#endif


#ifdef __CUDACC__
#define _constant_ /*__constant__*/  __device__ const
#else
#define _constant_ static const
#endif

#include "rysq_roots0.h"
#include "rysq_roots1.h"
#include "rysq_roots2.h"
#include "rysq_roots3.h"


#ifndef __CUDACC__
#include "roots/asymptotic.hpp"
#include "roots/opq.h"
#include "externals/cxx/pretty_function.hpp"
#include <boost/math/special_functions/pow.hpp>
#endif

#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/at.hpp>

#include <iostream>
#include <stdexcept>

/**
   @brief Initializes Rys quadrature root finding algorithms
*/
int Rysq_roots_init();

/**
   @brief Finializes Rys quadrature root finding algorithms
*/
int Rysq_roots_finalize();


/**
   @brief Finds the Rys quadrature roots and weights
   @param N Number of roots and weights to find, must be at least 0
   @param X The X value, cf. DRK
   @param[out] t2 Roots
   @param[out] W Weights
*/
void Rysq_roots(int N, double X, double *t2, double *W);

namespace rysq {


    template<int N, typename T>
    static void roots(T X, T *t2, T *W);

	     // template<int N> __device__
    // static inline void roots(double X, double *t2, double *W) {
    // 	stieltjes<N>(X, t2, W);
    // }

    template<int N> __device__
    static inline void roots(double X, double *t2, double *W,
			     const unsigned short thread);	

    template<> __device__
    inline void roots<0>(double X, double *t2, double *W) {
	W[0] = Rysq_roots0(X);
    }

    template<> __device__
    inline void roots<1>(double X, double *t2, double *W) {
	Rysq_roots1(X, t2, W);
    }





// #include "rysq_roots2.h"
    // #include "rysq_roots3.h"
#include "rysq_roots4.h"
#include "rysq_roots5.h"

    template<> __device__
    inline void roots<2>(double X, double *t2, double *W) {
	rysq::roots2::evaluate(X, t2, W);
    }

    template<> __device__
    inline void roots<3>(double X, double *t2, double *W) {
    	rysq::roots3::evaluate(X, t2, W);
    }

    template<> __device__
    inline void roots<4>(double X, double *t2, double *W) {
    	roots4::evaluate(X, t2, W);
    }

    template<> __device__
    inline void roots<5>(double X, double *t2, double *W) {
    	roots5::evaluate(X, t2, W);
    }

    template<> __device__
    inline void roots<2>(double X, double *t2, double *W,
    			 const unsigned short thread) {
    	roots2::evaluate(X, t2, W, thread);
    }

    template<> __device__
    inline void roots<3>(double X, double *t2, double *W,
    			 const unsigned short thread) {
    	roots3::evaluate(X, t2, W, thread);
    }

    template<> __device__
    inline void roots<4>(double X, double *t2, double *W,
			 const unsigned short thread) {
	roots4::evaluate(X, t2, W, thread);
    }

    template<> __device__
    inline void roots<5>(double X, double *t2, double *W,
			 const unsigned short thread) {
	roots5::evaluate(X, t2, W, thread);
    }

#ifndef __CUDACC__

    /**
       @brief Auxiliary quadrature size corresponding to a number of roots
    */
    template<size_t n>
    struct STIELTJES {
    	typedef typename boost::mpl::vector_c
    	<int,0, 1, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 7>::type index_vector;
    	static const size_t index = boost::mpl::at_c<index_vector,n-1>::type::value;
    	static const size_t N = 20 + index*5;
    };

    static const int STIELTJES_N[] = { 20, 25, 30, 35, 40, 45, 50, 55 }; 

    /**
       @brief Auxiliary quadrature index corresponding to a number of roots
    */
    static const int STIELTJES_INDEX[] = { 0, 1, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 7 };

    /**
       @brief Auxiliary quadrature roots
    */ 
    extern double *STIELTJES_X[8]; 

    /**
       @brief Auxiliary quadrature weights
    */
    extern double *STIELTJES_W[8]; 

    void stieltjes_init();

    void stieltjes_finalize();

    /**
       @brief Evaluates Rys quadrature roots and weights using Stieltjes method
       @param n Number of roots
       @param X X value
       @param[out] x Roots
       @param[out] w Weights 
    */   
    template<int n, typename T>
    void stieltjes(T X, T *x, T *w) {
    
	// int is = STIELTJES_INDEX[n-1];
	// int N = STIELTJES_N[is];
	static const int N = STIELTJES<n>::N;
	static const int is = STIELTJES<n>::index;

	T rg[N], wg[N];
	T a[n], b[n];
	
	using boost::math::pow;
	// if underflow is likely, scale grid weights (shift exponent)
	for(int i = 0; i < N; ++i) {
	    T t2 = pow<2>(STIELTJES_X[is][i]);
	    rg[i] = t2;
	    wg[i] = STIELTJES_W[is][i]*exp(-X*t2);
	}

	int status = 0;
	status = opq::stieltjes(n, N, rg, wg, a, b);
	if (status != 0) {
	    throw std::runtime_error(PRETTY_FUNCTION("opq::stieltjes",
						     " returned ", status));
	}

	status = opq::coefficients(n, a, b, x, w);
	if (status != 0) {
	    throw std::runtime_error(PRETTY_FUNCTION("opq::coefficients",
						     " returned ", status));
	}
		     
	//    for(int i = 0; i < n; ++i) {
	//	x[i] = x[i]/(1.0 - x[i]);
	// }
 
    }

    template<int N, typename T>
    static void roots(T X, T *t2, T *W) {
	if (X > 500) asymptotic::roots<N>(X, t2, W);
	else stieltjes<N>(X, t2, W);
    }

#endif // __CUDACC__

}

#undef _constant_

#ifndef __CUDACC__
#undef __constant__
#undef __device__
#endif

#endif //  RYSQ_ROOTS_H
