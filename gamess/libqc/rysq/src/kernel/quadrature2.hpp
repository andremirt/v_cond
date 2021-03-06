#ifndef _RYSQ_KERNEL_QUADRATURE2_HPP_
#define _RYSQ_KERNEL_QUADRATURE2_HPP_

#include "kernel/quadrature.hpp"
#include "kernel/quadrature2-impl.hpp"

#include "roots/rysq_roots.h"

BEGIN_NAMESPACE(rysq, kernel, quadrature)


struct coefficients {
    template<int N>
    static void load(const Shell &s, double C[][N], double Cmax[]) {
	for(int k = 0; k < s.K; ++k) {
	    Cmax[k] = 0.0;
	    for(int i = 0; i < N; ++i) {
		C[k][i] = s(k,i);
		Cmax[k] = std::max(Cmax[k], fabs(C[k][i]));
	    }
	}
    }	
};

template<class braket>
void apply(const Quartet<Shell> &quartet,
	   const Vector<3> &ri, const Vector<3> &rj,
	   const Vector<3> &rk, const Vector<3> &rl,
	   double (&Q)[braket::size], int transpose,
	   double cutoff) {

    typedef quadrature::impl<braket> impl;

    const Shell &a = quartet[0];
    const Shell &b = quartet[1];
    const Shell &c = quartet[2];
    const Shell &d = quartet[3];

    static const int L = braket::L;
    static const int N = (L > 0) ? L/2 + 1 : 0;

    Vector<3> rij = ri - rj;
    Vector<3> rkl = rk - rl;

    double rij2 = rij.inner();
    double rkl2 = rkl.inner();

    double eij[a.K*b.K] __attribute__ ((aligned(16)));
    for(int Kj = 0, Kij = 0; Kj < b.K; ++Kj) {
	for(int Ki = 0; Ki < a.K; ++Ki, ++Kij) {
	    double ai = a(Ki);
	    double aj = b(Kj);
	    eij[Kij] = exp(-ai*aj*rij2/(ai + aj));
	}
    }

    double Ci[a.K][braket::nci] __attribute__ ((aligned(16)));
    double Cj[b.K][braket::ncj] __attribute__ ((aligned(16)));
    double Ci_max[a.K] __attribute__ ((aligned(16)));
    double Cj_max[b.K] __attribute__ ((aligned(16)));
    coefficients::load(a, Ci, Ci_max);
    coefficients::load(b, Cj, Cj_max);

    for (int Kl = 0; Kl < d.K; ++Kl) {
	for (int Kk = 0; Kk < c.K; ++Kk) {
	    double ak = c(Kk);
	    double al = d(Kl);
	    double B = ak + al;

	    Vector<3> rB = Vector<3>::center(ak, rk, al, rl);
	    Vector<3> rBk = rB - rk;

	    double ekl = exp(-ak*al*rkl2/B);

	    static const int nckl = braket::nck*braket::ncl;
	    double Ckl[nckl] __attribute__ ((aligned(16)));
	    double Ckl_max = 0.0;

	    for(int l = 0; l < braket::ncl; ++l) {
		for(int k = 0; k < braket::nck; ++k) {
		    double q = c(Kk,k)*d(Kl,l);
		    Ckl[k+l*braket::nck] = q;
		    Ckl_max = std::max(fabs(q), Ckl_max);
		}
	    }
		    
	    for(int Kj = 0; Kj < b.K; ++Kj) {
		for(int Ki = 0; Ki < a.K; ++Ki) {
		    double ai = a(Ki);
		    double aj = b(Kj); 
		    double A = ai + aj;
		    double e = eij[Ki+Kj*a.K]*ekl;
		    double eAB = e/(A*B*sqrt(A+B));

		    if (eAB*Ci_max[Ki]*Cj_max[Kj]*Ckl_max < cutoff) continue;

		    double C[braket::nc] __attribute__ ((aligned(16)));
		    for(int kl = 0, ijkl = 0; kl < nckl; ++kl) {
			double q = eAB*Ckl[kl];
			for(int j = 0; j < braket::ncj; ++j) {
			    for(int i = 0; i < braket::nci; ++i, ++ijkl) {
				C[ijkl] = q*Ci[Ki][i]*Cj[Kj][j];
			    }
			}
		    }

		    Vector<3> rA = Vector<3>::center(ai, ri, aj, rj);
		    Vector<3> rAi = rA - ri;
		    Vector<3> rAB = rA - rB;

		    double rho = (A*B)/(A + B);
		    double X =  rho*rAB.inner();
		    Vector<N> t2;
		    Vector<N+(L == 0)> W;
		    rysq::roots<N,double>(X, t2, W);

		    t2 /= (A + B);
		    impl::eval(A, B, rAi, rBk, rAB, rij, rkl, t2, W, C, Q);

		}
	    }
	}
    }
    impl::reorder(Q);

    if (transpose) {
	rysq::transpose<braket::A::size, braket::B::size,
	    braket::C::size, braket::D::size>(transpose, Q);
    }

}

END_NAMESPACE(rysq, kernel, quadrature)

#endif /* _RYSQ_KERNEL_QUADRATURE2_HPP_ */
