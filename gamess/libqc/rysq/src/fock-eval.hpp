#ifndef _RYSQ_FOCK_EVAL_HPP_
#define _RYSQ_FOCK_EVAL_HPP_

#include <rysq/core.hpp>
#include <rysq/fock.hpp>
#include "meta.hpp"

BEGIN_NAMESPACE(rysq, fock)

typedef Fock::density_set density_set;
typedef Fock::fock_set fock_set;

template<class braket>
inline void eval(const density_set &D, const fock_set &F,
		 const double *__restrict Eri, const Parameters &parameters) {

    using util::add;
    using util::copy;

    double scale = parameters.scale;
    double scale2 = parameters.scale2;

    static const int ni = braket::A::size;
    static const int nj = braket::B::size;
    static const int nk = braket::C::size;
    static const int nl = braket::D::size;


    double Dij[ni*nj] __attribute__ ((aligned(16)));
    double Dkl[nk*nl] __attribute__ ((aligned(16)));
    double Dik[ni*nk] __attribute__ ((aligned(16)));
    double Dil[ni*nl] __attribute__ ((aligned(16)));
    double Djk[nj*nk] __attribute__ ((aligned(16)));
    double Djl[nj*nl] __attribute__ ((aligned(16)));

    copy(ni*nj, Dij, D[0]);
    copy(nk*nl, Dkl, D[1]);
    copy(ni*nk, Dik, D[2]);
    copy(ni*nl, Dil, D[3]);
    copy(nj*nk, Djk, D[4]);
    copy(nj*nl, Djl, D[5]);

    double Fij[ni*nj] __attribute__ ((aligned(16))) = { 0.0 };
    double Fkl[nk*nl] __attribute__ ((aligned(16))) = { 0.0 };
    double Fik[ni*nk] __attribute__ ((aligned(16))) = { 0.0 };
    double Fil[ni*nl] __attribute__ ((aligned(16))) = { 0.0 };
    double Fjk[nj*nk] __attribute__ ((aligned(16))) = { 0.0 };
    double Fjl[nj*nl] __attribute__ ((aligned(16))) = { 0.0 };

    for (int l = 0, kl = 0, ijkl = 0; l < nl; ++l) {

	for (int k = 0; k < nk; ++k, ++kl) {
	    int jk = k*nj;
	    int jl = l*nj;

	    for (int j = 0, ij = 0; j < nj; ++j, ++jk, ++jl) {
		int ik = k*ni;
		int il = l*ni;

		for (int i = 0; i < ni; ++i, ++ij, ++ik, ++il, ++ijkl) {
		    //double q = scale*Eri[ijkl];
		    double q = Eri[ijkl];
		    //if(fabs(q) < 1e-11) continue;

		    Fij[ij] += Dkl[kl]*q;
		    Fkl[kl] += Dij[ij]*q;
		    Fik[ik] += Djl[jl]*q;
		    Fil[il] += Djk[jk]*q;
		    Fjk[jk] += Dil[il]*q;
		    Fjl[jl] += Dik[ik]*q;
		}
	    }

	}
    }
    add(ni*nj, F[0], scale*scale2, Fij);
    add(nk*nl, F[1], scale*scale2, Fkl);
    add(ni*nk, F[2], -scale, Fik);
    add(ni*nl, F[3], -scale, Fil);
    add(nj*nk, F[4], -scale, Fjk);
    add(nj*nl, F[5], -scale, Fjl);
}

template<int ni, int nj, int nk, int nl>
inline void eval(const Data::density_type &D, const Data::fock_type &F,
		 const double *__restrict Q, double scale) {

    const double * __restrict Dij = D[0];
    const double * __restrict Dkl = D[1];
    const double * __restrict Dik = D[2];
    const double * __restrict Dil = D[3];
    const double * __restrict Djk = D[4];
    const double * __restrict Djl = D[5];

    double * __restrict Fij = F[0];
    double * __restrict Fkl = F[1];
    double * __restrict Fik = F[2];
    double * __restrict Fil = F[3];
    double * __restrict Fjk = F[4];
    double * __restrict Fjl = F[5];

    for (int l = 0, kl = 0, ijkl = 0; l < nl; ++l) {

	for (int k = 0; k < nk; ++k, ++kl) {
	    int jk = k*nj;
	    int jl = l*nj;

	    for (int j = 0, ij = 0; j < nj; ++j, ++jk, ++jl) {
		int ik = k*ni;
		int il = l*ni;

		for (int i = 0; i < ni; ++i, ++ij, ++ik, ++il, ++ijkl) {
		    double q = scale*Q[ijkl];
		    // std::cout << Dkl[kl] << std::endl;
		    Fij[ij] += Dkl[kl]*q;
		    Fkl[kl] += Dij[ij]*q;
		    Fik[ik] += Djl[jl]*q;
		    Fil[il] += Djk[jk]*q;
		    Fjk[jk] += Dil[il]*q;
		    Fjl[jl] += Dik[ik]*q;
		}
	    }

	}
    }
}



template<int M,int N>
void eval(int k, int l, int kl,
	  const Data::density_type &D, const Data::fock_type &F,
	  const double * __restrict I, double scale = 1) {

    // bool test = false;
    // for (int i = 0; i < M*N; ++i) {
    // 	test += (scale*I[i] > 1e-10);
    // 	if (test) break;
    // }
    // if (!test) return;

    const double * __restrict Dij = D[0];
    const double * __restrict Dkl = D[1];
    const double * __restrict Dik = D[2];
    const double * __restrict Dil = D[3];
    const double * __restrict Djk = D[4];
    const double * __restrict Djl = D[5];

    double * __restrict Fij = F[0];
    double * __restrict Fkl = F[1];
    double * __restrict Fik = F[2];
    double * __restrict Fil = F[3];
    double * __restrict Fjk = F[4];
    double * __restrict Fjl = F[5];

    int jk = k*N;
    int jl = l*N;
    for (int j = 0, ij = 0; j < N; ++j, ++jk, ++jl) {
	    // std::cout << Fjk + jk << std::endl;

	// bool test = false;
	// for (int i = 0; i < M; ++i) {
	//     test += (scale*I[ij+i] > 1e-10);
	// }
	// if (!test) continue;

    	int ik = k*M;
    	int il = l*M;
    	for (int i = 0; i < M; ++i, ++ij, ++ik, ++il) {
    	    double q = scale*I[ij];


    	    Fij[ij] += Dkl[kl]*q;
    	    Fkl[kl] += Dij[ij]*q;
    	    Fik[ik] += Djl[jl]*q;
    	    Fil[il] += Djk[jk]*q;
    	    Fjk[jk] += Dil[il]*q;
	    Fjl[jl] += Dik[ik]*q;
    	}
    }

}



template<rysq::type A, rysq::type B>
void eval(const rysq::type &c, const rysq::type &d,
	  const density_set &D, const fock_set &F,
	  double scale, const double *eri) {

    static const int ni = rysq::meta::shell<A>::size;
    static const int nj = rysq::meta::shell<B>::size;
    const int nk = rysq::shell(c).size;
    const int nl = rysq::shell(d).size;

    bool test = false;
    for (int i = 0; i < ni*nj; ++i) {
	test += (scale*eri[i] > 1e-10);
    }
    if (!test) return;

    for (int l = 0, kl = 0, ijkl = 0; l < nl; ++l) {

	for (int k = 0; k < nk; ++k, ++kl) {
	    int jk = k*nj;
	    int jl = l*nj;

	    for (int j = 0, ij = 0; j < nj; ++j, ++jk, ++jl) {
		int ik = k*ni;
		int il = l*ni;

		for (int i = 0; i < ni; ++i, ++ij, ++ik, ++il, ++ijkl) {
		    //double q = scale*Eri[ijkl];
		    double q = scale*eri[ijkl];
		    //if(fabs(q) < 1e-11) continue;

		    F[0][ij] += 4.0*D[1][kl]*q;
		    F[1][kl] += 4.0*D[0][ij]*q;
		    F[2][ik] -= D[5][jl]*q;
		    F[3][il] -= D[4][jk]*q;
		    F[4][jk] -= D[3][il]*q;
		    F[5][jl] -= D[2][ik]*q;
		}
	    }

	}
    }
}

// template<class bra, int N>
// inline void eval(const Shell &c, const Shell &d,
// 		 int K, const double* const *C, const double *Ints2d,
// 		 const density_set  &D,
// 		 const fock_set &F) {
     
//     typedef double T;
//     static const int ldN = (sizeof(T) > 8) ? N : N + N%2;
//     static const int Nij = ldN*(bra::A::L+1)*(bra::B::L+1);
//     const int dim2d = Nij*(c.L+1)*(d.L+1);

//     const T *Ix = &Ints2d[0*dim2d];
//     const T *Iy = &Ints2d[1*dim2d];
//     const T *Iz = &Ints2d[2*dim2d];

//     int spk = (c.type < 0);
//     int spl = (d.type < 0);

//     for(int l = d.begin(), kl = 0; l < d.end(); ++l) {
// 	const int lx = (c.L+1)*LX[l];
// 	const int ly = (c.L+1)*LY[l];
// 	const int lz = (c.L+1)*LZ[l];

// 	const int lsp = (spl && l) << spk;

// 	for(int k = c.begin(); k < c.end(); ++k, ++kl) {
// 	    const double *Ckl = C[(spk && k) + lsp];

// 	    const int klx = Nij*(lx + LX[k]);
// 	    const int kly = Nij*(ly + LY[k]);
// 	    const int klz = Nij*(lz + LZ[k]);

// 	    static const int ni = bra::A::size;
// 	    static const int nj = bra::B::size;
// 	    double I[ni*nj] __attribute__((aligned(16))) = { 0.0 };
	    
// 	    int flags = 0;
// 	    double screen = 0.0;
// 	    typedef kernels::quadrature<bra> quadrature;
// 	    quadrature::template apply<N>(flags, screen, K, Ckl, dim2d,
// 				       &Ix[klx], &Iy[kly], &Iz[klz], 1.0, I);

// 	    fock::eval<ni,nj>(k - c.begin(), l - d.begin(), kl, D, F, I);
// 			      // D[0], D[1], D[2], D[3], D[4], D[5],
// 			      // F[0], F[1], F[2], F[3], F[4], F[5], I);
// 	}
//     }

// }

END_NAMESPACE(rysq, fock)


#endif /* _RYSQ_FOCK_EVAL_HPP_ */
