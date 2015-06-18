#ifndef _RYSQ_TRANSFER_HPP_
#define _RYSQ_TRANSFER_HPP_


#include <assert.h>
#include <math.h>

namespace rysq {

    template<typename T, size_t N, size_t A = 16>
    struct align {
	static const size_t value = N%2;//(A/sizeof(T));
    };

    template<int A, int B, int N, typename S, typename T>
    static inline void transfer(int CD, S dq, S *__restrict G, T *__restrict I) {

	static const int ldN = N + N%2;
#define G(n,ab) G[(n) + ldN*((ab))]
#define I(n,a,b) I[(n) + ldN*((a)+(A+1)*((b)))]

	for(int kl = 0; kl <= CD; ++kl) {
	
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector aligned
#endif
	    for(int a = 0; a < N; ++a) {
		for (int i = 0; i <= A; ++i) {
		    I(a,i,0) = G(a,i);
		}

		for (int j = 1; j <= B; ++j) {
		    for (int i = 0; i <= A; ++i) {
			G(a,i) = dq*G(a,i) + G(a,i+1);
			I(a,i,j) = G(a,i);
		    }
		    for (int i = A+1; i <= A+B-j; ++i) {
			G(a,i) = dq*G(a,i) + G(a,i+1);
		    }
		}
	    }
	    G += ldN*(A+B+1);
	    I += ldN*(A+1)*(B+1);
	}

#undef G
#undef I
    }

    template<int N, typename S, typename T>
    static inline void transfer(int A, int B, S dq, S *__restrict G, T *__restrict I) {

	static const int ldN = N + N%2;
#define G(n,ij) G[(n) + ldN*((ij))]
#define I(n,i,j) I[(n) + ldN*((i)+(A+1)*((j)))]



	for (int i = 0; i <= A; ++i) {
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector aligned
#endif
	    for(int a = 0; a < N; ++a) {
		I(a,i,0) = G(a,i);
	    }
	}

	for (int j = 1; j <= B; ++j) {
	    for (int i = 0; i <= A; ++i) {
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector aligned
#endif
		for(int a = 0; a < N; ++a) {
		    G(a,i) = dq*G(a,i) + G(a,i+1);
		    I(a,i,j) = G(a,i);
		}
	    }

	    for (int i = A+1; i <= A+B-j; ++i) {
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector aligned
#endif
		for(int a = 0; a < N; ++a) {
		    G(a,i) = dq*G(a,i) + G(a,i+1);	    
		}
	    }
	}

#undef G
#undef I
    }

    template<int A, int B, int N, typename T, typename U>
    static void transfer(int C, int D, U dij, U dkl,
			 U *G, T *I, U *mem16) {

	static const int ldN = N + align<T,N>::value;

	if (B + D == 0) return;
	else if (D == 0) transfer<A,B,N>(C, dij, G, I);
	else if (B == 0) transfer<(A+1)*ldN>(C, D, dkl, G, I);
	else {
	    transfer<A,B,N>(C+D, dij, G, mem16);
	    transfer<(A+1)*(B+1)*ldN>(C, D, dkl, mem16, I);
	}
    }

}


#endif /* _RYSQ_TRANSFER_HPP_ */
