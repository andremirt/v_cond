/**  
@file 
@warning Automatically Generated
*/
/**  
 @warning AUTOMATICALLY GENERATED
*/




#ifndef _RYSQ_KERNEL_QUADRATURE2_IMPL_HPP_
#define _RYSQ_KERNEL_QUADRATURE2_IMPL_HPP_

#include "rysq/config.hpp"
#include "kernel/eri.hpp"
#include "meta.hpp"
#include "vector.hpp"

BEGIN_NAMESPACE(rysq, kernel, quadrature)

namespace recurrence {

    static inline double coefficient(double A1, double B, double t2) {
	return 0.5*A1*(1.0 - B*t2);
    }
    
    template<int q>
    static inline double coefficient(const Vector<3> &rAi, double B,
				     const Vector<3> &rAB, double t2) {
	return rAi[q] - B*rAB[q]*t2;
    }

    template<int q, int N>
    static inline Vector<N> coefficient(const Vector<3> &rAi, double B,
					const Vector<3> &rAB, double (&t2)[N]) {
	Vector<N> C;
	for (int a = 0; a < N; ++a) {
	    C[a] = coefficient<q>(rAi, B, rAB, t2[a]);
	}
    }

}

namespace mpl = boost::mpl;


template<>
struct impl<meta::braket<rysq::P, rysq::P, rysq::P, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[1], double (&I)[27]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f1 = (B00 + Dz*Iz);
	    double f11 = (B10 + Cz*Iz);
	    double f3 = (Cy*Iy + B10);
	    double f7 = (Cx*Ix + B10);

	    I[0] += C[0]*W[a]*((Dz*(B10 + Cz*Iz) + B00*(2*Cz + zij)));
	    I[1] += C[0]*W[a]*(Cy*(Dx*xij + Qx));
	    I[2] += C[0]*W[a]*(Cz*(Dx*xij + Qx));
	    I[3] += C[0]*W[a]*((Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx)));
	    I[4] += C[0]*W[a]*((B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10)));
	    I[5] += C[0]*W[a]*(Cx*(Dy*yij + Qy));
	    I[6] += C[0]*W[a]*(Cz*(Dy*yij + Qy));
	    I[7] += C[0]*W[a]*(Cz*Dy*Ix);
	    I[8] += C[0]*W[a]*(Cx*Dy*Iz);
	    I[9] += C[0]*W[a]*(Cy*Dx*Iz);
	    I[10] += C[0]*W[a]*(Cy*Dz*Ix);
	    I[11] += C[0]*W[a]*(Cx*Dz*Iy);
	    I[12] += C[0]*W[a]*(Cz*Dx*Iy);
	    I[13] += C[0]*W[a]*(Dx*f3);
	    I[14] += C[0]*W[a]*(Dx*f11);
	    I[15] += C[0]*W[a]*(Dy*f11);
	    I[16] += C[0]*W[a]*(Dy*f7);
	    I[17] += C[0]*W[a]*(Dz*f7);
	    I[18] += C[0]*W[a]*(Dz*f3);
	    I[19] += C[0]*W[a]*(Ix*Qz);
	    I[20] += C[0]*W[a]*(Ix*Qy);
	    I[21] += C[0]*W[a]*(Iz*Qy);
	    I[22] += C[0]*W[a]*(Iz*Qx);
	    I[23] += C[0]*W[a]*(Iy*Qx);
	    I[24] += C[0]*W[a]*(Iy*Qz);
	    I[25] += C[0]*W[a]*(Cx*f1);
	    I[26] += C[0]*W[a]*(Cy*f1);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[27]) {
	double T[27];
	for (int i = 0; i < 27; ++i) {
	    T[i] = I[i];
	}
	I[26] = T[0];
	I[1] = T[1];
	I[2] = T[2];
	I[0] = T[3];
	I[13] = T[4];
	I[12] = T[5];
	I[14] = T[6];
	I[11] = T[7];
	I[15] = T[8];
	I[7] = T[9];
	I[19] = T[10];
	I[21] = T[11];
	I[5] = T[12];
	I[4] = T[13];
	I[8] = T[14];
	I[17] = T[15];
	I[9] = T[16];
	I[18] = T[17];
	I[22] = T[18];
	I[20] = T[19];
	I[10] = T[20];
	I[16] = T[21];
	I[6] = T[22];
	I[3] = T[23];
	I[23] = T[24];
	I[24] = T[25];
	I[25] = T[26];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::P, rysq::S, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[1], double (&I)[30]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 

// 	for (int a = 0; a < N; ++a) {
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);


	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double f0 = (3*B10 + pow(Cx,2));
	    double f11 = (3*B10 + pow(Cz,2));
	    double f12 = (3*B10 + pow(Cy,2));
	    double f13 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f2 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));
	    double f5 = (Cy*Iy + B10);
	    double f6 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f8 = 3*pow(B10,2);

	    I[0] += C[0]*W[a]*(Cx*Cy*(Cz*zij + Pz));
	    I[1] += C[0]*W[a]*(Cy*Cz*(Px + Cx*xij));
	    I[2] += C[0]*W[a]*((3*B10*Cx*(xij + 2*Cx) + Ix*pow(Cx,3) + f8));
	    I[3] += C[0]*W[a]*((Iy*pow(Cy,3) + f8 + 3*B10*Cy*(yij + 2*Cy)));
	    I[4] += C[0]*W[a]*((3*B10*Cz*(2*Cz + zij) + f8 + Iz*pow(Cz,3)));
	    I[5] += C[0]*W[a]*(Px*(Cz*zij + Pz));
	    I[6] += C[0]*W[a]*(Py*(Cz*zij + Pz));
	    I[7] += C[0]*W[a]*(Cx*Cz*f5);
	    I[8] += C[0]*W[a]*(Cy*Ix*Pz);
	    I[9] += C[0]*W[a]*(Cy*Ix*f12);
	    I[10] += C[0]*W[a]*(Cy*Iz*f12);
	    I[11] += C[0]*W[a]*(Cy*Iz*Px);
	    I[12] += C[0]*W[a]*(Py*(Px + Cx*xij));
	    I[13] += C[0]*W[a]*(Pz*(Px + Cx*xij));
	    I[14] += C[0]*W[a]*(Cx*Iy*Pz);
	    I[15] += C[0]*W[a]*(Cx*Iy*f0);
	    I[16] += C[0]*W[a]*(Cx*Iz*f0);
	    I[17] += C[0]*W[a]*(Cx*Iz*Py);
	    I[18] += C[0]*W[a]*(Cz*Ix*Py);
	    I[19] += C[0]*W[a]*(Cz*Ix*f11);
	    I[20] += C[0]*W[a]*(Cz*Iy*f11);
	    I[21] += C[0]*W[a]*(Cz*Iy*Px);
	    I[22] += C[0]*W[a]*(Px*f5);
	    I[23] += C[0]*W[a]*(Pz*f5);
	    I[24] += C[0]*W[a]*(Cy*f6);
	    I[25] += C[0]*W[a]*(Cz*f6);
	    I[26] += C[0]*W[a]*(Cz*f2);
	    I[27] += C[0]*W[a]*(Cx*f2);
	    I[28] += C[0]*W[a]*(Cx*f13);
	    I[29] += C[0]*W[a]*(Cy*f13);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[30]) {
	double T[30];
	for (int i = 0; i < 30; ++i) {
	    T[i] = I[i];
	}
	I[29] = T[0];
	I[9] = T[1];
	I[0] = T[2];
	I[11] = T[3];
	I[22] = T[4];
	I[24] = T[5];
	I[26] = T[6];
	I[19] = T[7];
	I[8] = T[8];
	I[1] = T[9];
	I[21] = T[10];
	I[23] = T[11];
	I[5] = T[12];
	I[7] = T[13];
	I[17] = T[14];
	I[10] = T[15];
	I[20] = T[16];
	I[25] = T[17];
	I[6] = T[18];
	I[2] = T[19];
	I[12] = T[20];
	I[14] = T[21];
	I[13] = T[22];
	I[18] = T[23];
	I[3] = T[24];
	I[4] = T[25];
	I[16] = T[26];
	I[15] = T[27];
	I[27] = T[28];
	I[28] = T[29];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::S, rysq::S, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[1], double (&I)[10]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;



// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 

// 	for (int a = 0; a < N; ++a) {
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);


	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));

	    I[0] += C[0]*W[a]*(Cx*Cy*Cz);
	    I[1] += C[0]*W[a]*(Cy*Px);
	    I[2] += C[0]*W[a]*(Cz*Px);
	    I[3] += C[0]*W[a]*(Cx*Py);
	    I[4] += C[0]*W[a]*(Cz*Py);
	    I[5] += C[0]*W[a]*((3*B10*Cz + pow(Cz,3)));
	    I[6] += C[0]*W[a]*((3*B10*Cy + pow(Cy,3)));
	    I[7] += C[0]*W[a]*((3*B10*Cx + pow(Cx,3)));
	    I[8] += C[0]*W[a]*(Cx*Pz);
	    I[9] += C[0]*W[a]*(Cy*Pz);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[10]) {
	double T[10];
	for (int i = 0; i < 10; ++i) {
	    T[i] = I[i];
	}
	I[9] = T[0];
	I[3] = T[1];
	I[4] = T[2];
	I[5] = T[3];
	I[6] = T[4];
	I[2] = T[5];
	I[1] = T[6];
	I[0] = T[7];
	I[7] = T[8];
	I[8] = T[9];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::P, rysq::P, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[1], double (&I)[54]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f1 = (B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10));
	    double f10 = (Dy*Iy + B00);
	    double f11 = (Cy*Iy + B10);
	    double f12 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f14 = (Dx*Px + 2*B00*Cx);
	    double f15 = (Dx*Ix + B00);
	    double f22 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f3 = 3*B00*B10;
	    double f4 = (Dy*Py + 2*B00*Cy);
	    double f5 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));
	    double f6 = (Dz*Pz + 2*B00*Cz);
	    double f8 = (Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx));

	    I[0] += C[0]*W[a]*((Dy*Iy*pow(Cy,2) + f3 + B00*Cy*(3*Cy + 2*yij) + B10*Dy*(3*Cy + yij)));
	    I[1] += C[0]*W[a]*(Cx*Cy*(Dz*zij + Qz));
	    I[2] += C[0]*W[a]*(Cy*Dz*(Px + Cx*xij));
	    I[3] += C[0]*W[a]*(Cz*Dy*(Px + Cx*xij));
	    I[4] += C[0]*W[a]*(Qy*(Px + Cx*xij));
	    I[5] += C[0]*W[a]*(Qz*(Px + Cx*xij));
	    I[6] += C[0]*W[a]*((Dx*Ix*pow(Cx,2) + B00*Cx*(3*Cx + 2*xij) + B10*Dx*(3*Cx + xij) + f3));
	    I[7] += C[0]*W[a]*((B00*Cz*(3*Cz + 2*zij) + Dz*Iz*pow(Cz,2) + B10*Dz*(3*Cz + zij) + f3));
	    I[8] += C[0]*W[a]*(Cx*(f6 + Qz*zij));
	    I[9] += C[0]*W[a]*(Cy*(f6 + Qz*zij));
	    I[10] += C[0]*W[a]*(Cy*Dx*(Cz*zij + Pz));
	    I[11] += C[0]*W[a]*(Qy*(Cz*zij + Pz));
	    I[12] += C[0]*W[a]*(Qx*(Cz*zij + Pz));
	    I[13] += C[0]*W[a]*(Cx*Dy*(Cz*zij + Pz));
	    I[14] += C[0]*W[a]*(Dy*Ix*Pz);
	    I[15] += C[0]*W[a]*(Cz*Ix*Qy);
	    I[16] += C[0]*W[a]*(Cx*Cz*f10);
	    I[17] += C[0]*W[a]*(Cz*Dx*f11);
	    I[18] += C[0]*W[a]*(Dx*Iy*Pz);
	    I[19] += C[0]*W[a]*(Dx*Iz*Py);
	    I[20] += C[0]*W[a]*(Dy*Iz*Px);
	    I[21] += C[0]*W[a]*(Cx*Iz*Qy);
	    I[22] += C[0]*W[a]*(Cy*Iz*Qx);
	    I[23] += C[0]*W[a]*(Cy*Ix*Qz);
	    I[24] += C[0]*W[a]*(Dz*Ix*Py);
	    I[25] += C[0]*W[a]*(Py*(Dz*zij + Qz));
	    I[26] += C[0]*W[a]*(Px*(Dz*zij + Qz));
	    I[27] += C[0]*W[a]*(Dz*Iy*Px);
	    I[28] += C[0]*W[a]*(Cx*Dz*f11);
	    I[29] += C[0]*W[a]*(Cx*Iy*Qz);
	    I[30] += C[0]*W[a]*(Cz*Iy*Qx);
	    I[31] += C[0]*W[a]*(Cy*Cz*f15);
	    I[32] += C[0]*W[a]*(Py*f15);
	    I[33] += C[0]*W[a]*(Pz*f15);
	    I[34] += C[0]*W[a]*(Cy*f8);
	    I[35] += C[0]*W[a]*(Cz*f8);
	    I[36] += C[0]*W[a]*(Px*f10);
	    I[37] += C[0]*W[a]*(Pz*f10);
	    I[38] += C[0]*W[a]*(Cx*f1);
	    I[39] += C[0]*W[a]*(Cz*f1);
	    I[40] += C[0]*W[a]*(Dx*f5);
	    I[41] += C[0]*W[a]*(Dx*f22);
	    I[42] += C[0]*W[a]*(Dy*f22);
	    I[43] += C[0]*W[a]*(Dy*f12);
	    I[44] += C[0]*W[a]*(Dz*f12);
	    I[45] += C[0]*W[a]*(Dz*f5);
	    I[46] += C[0]*W[a]*(Ix*f6);
	    I[47] += C[0]*W[a]*(Ix*f4);
	    I[48] += C[0]*W[a]*(Iz*f4);
	    I[49] += C[0]*W[a]*(Iz*f14);
	    I[50] += C[0]*W[a]*(Iy*f14);
	    I[51] += C[0]*W[a]*(Iy*f6);
	    I[52] += C[0]*W[a]*(Qx*f11);
	    I[53] += C[0]*W[a]*(Qz*f11);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[54]) {
	double T[54];
	for (int i = 0; i < 54; ++i) {
	    T[i] = I[i];
	}
	I[25] = T[0];
	I[51] = T[1];
	I[39] = T[2];
	I[22] = T[3];
	I[21] = T[4];
	I[40] = T[5];
	I[0] = T[6];
	I[50] = T[7];
	I[52] = T[8];
	I[53] = T[9];
	I[17] = T[10];
	I[35] = T[11];
	I[16] = T[12];
	I[34] = T[13];
	I[20] = T[14];
	I[23] = T[15];
	I[28] = T[16];
	I[11] = T[17];
	I[8] = T[18];
	I[13] = T[19];
	I[30] = T[20];
	I[33] = T[21];
	I[15] = T[22];
	I[41] = T[23];
	I[37] = T[24];
	I[49] = T[25];
	I[48] = T[26];
	I[42] = T[27];
	I[45] = T[28];
	I[46] = T[29];
	I[10] = T[30];
	I[5] = T[31];
	I[1] = T[32];
	I[2] = T[33];
	I[3] = T[34];
	I[4] = T[35];
	I[24] = T[36];
	I[26] = T[37];
	I[27] = T[38];
	I[29] = T[39];
	I[7] = T[40];
	I[14] = T[41];
	I[32] = T[42];
	I[18] = T[43];
	I[36] = T[44];
	I[43] = T[45];
	I[38] = T[46];
	I[19] = T[47];
	I[31] = T[48];
	I[12] = T[49];
	I[6] = T[50];
	I[44] = T[51];
	I[9] = T[52];
	I[47] = T[53];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::S, rysq::D, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[1], double (&I)[60]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;



// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double Rx = (B01 + pow(Dx,2));
	    double Ry = (B01 + pow(Dy,2));
	    double Rz = (pow(Dz,2) + B01);
	    double f0 = (3*B10 + pow(Cx,2));
	    double f1 = (Dy*Py + 2*B00*Cy);
	    double f11 = (2*B00*Dx + Cx*Rx);
	    double f12 = (Dx*Px + 2*B00*Cx);
	    double f14 = (Px*Rx + 2*pow(B00,2) + 4*B00*Cx*Dx);
	    double f17 = (Cx*Dx*(3*B10 + pow(Cx,2)) + 3*B00*Px);
	    double f19 = (3*B10 + pow(Cz,2));
	    double f2 = (Dz*Pz + 2*B00*Cz);
	    double f20 = (3*B10 + pow(Cy,2));
	    double f21 = (2*pow(B00,2) + Pz*Rz + 4*B00*Cz*Dz);
	    double f22 = (2*B00*Dy + Cy*Ry);
	    double f3 = (4*B00*Cy*Dy + Py*Ry + 2*pow(B00,2));
	    double f4 = (Cy*Dy*(3*B10 + pow(Cy,2)) + 3*B00*Py);
	    double f7 = (2*B00*Dz + Cz*Rz);
	    double f8 = (3*B00*Pz + Cz*Dz*(3*B10 + pow(Cz,2)));

	    I[0] += C[0]*W[a]*((6*Cx*pow(B00,2) + Cx*Rx*f0 + 6*B00*Dx*Px));
	    I[1] += C[0]*W[a]*((6*Cy*pow(B00,2) + Cy*Ry*f20 + 6*B00*Dy*Py));
	    I[2] += C[0]*W[a]*((6*B00*Dz*Pz + 6*Cz*pow(B00,2) + Cz*Rz*f19));
	    I[3] += C[0]*W[a]*(Cz*Dx*Dy*f19);
	    I[4] += C[0]*W[a]*(Cz*Dx*f1);
	    I[5] += C[0]*W[a]*(Cz*Dy*f12);
	    I[6] += C[0]*W[a]*(Dy*Pz*Qx);
	    I[7] += C[0]*W[a]*(Dx*Pz*Qy);
	    I[8] += C[0]*W[a]*(Cz*Qx*Qy);
	    I[9] += C[0]*W[a]*(Cx*Cz*f22);
	    I[10] += C[0]*W[a]*(Cx*Cy*f7);
	    I[11] += C[0]*W[a]*(Cy*Dx*Dz*f20);
	    I[12] += C[0]*W[a]*(Cy*Dz*f12);
	    I[13] += C[0]*W[a]*(Dz*Py*Qx);
	    I[14] += C[0]*W[a]*(Dx*Py*Qz);
	    I[15] += C[0]*W[a]*(Cy*Dx*f2);
	    I[16] += C[0]*W[a]*(Cy*Qx*Qz);
	    I[17] += C[0]*W[a]*(Cy*Rz*f20);
	    I[18] += C[0]*W[a]*(Cy*Px*Rz);
	    I[19] += C[0]*W[a]*(Dz*Px*Qy);
	    I[20] += C[0]*W[a]*(Dy*Px*Qz);
	    I[21] += C[0]*W[a]*(Cx*Dy*Dz*f0);
	    I[22] += C[0]*W[a]*(Cx*Dz*f1);
	    I[23] += C[0]*W[a]*(Cx*Dy*f2);
	    I[24] += C[0]*W[a]*(Cx*Qy*Qz);
	    I[25] += C[0]*W[a]*(Cx*Py*Rz);
	    I[26] += C[0]*W[a]*(Cx*Rz*f0);
	    I[27] += C[0]*W[a]*(Cx*Ry*f0);
	    I[28] += C[0]*W[a]*(Cx*Pz*Ry);
	    I[29] += C[0]*W[a]*(Cy*Pz*Rx);
	    I[30] += C[0]*W[a]*(Cy*Rx*f20);
	    I[31] += C[0]*W[a]*(Cy*Cz*f11);
	    I[32] += C[0]*W[a]*(Cz*Px*Ry);
	    I[33] += C[0]*W[a]*(Cz*Ry*f19);
	    I[34] += C[0]*W[a]*(Cz*Rx*f19);
	    I[35] += C[0]*W[a]*(Cz*Py*Rx);
	    I[36] += C[0]*W[a]*(Py*f11);
	    I[37] += C[0]*W[a]*(Pz*f11);
	    I[38] += C[0]*W[a]*(Pz*f22);
	    I[39] += C[0]*W[a]*(Px*f22);
	    I[40] += C[0]*W[a]*(Px*f7);
	    I[41] += C[0]*W[a]*(Py*f7);
	    I[42] += C[0]*W[a]*(Cy*f14);
	    I[43] += C[0]*W[a]*(Cz*f14);
	    I[44] += C[0]*W[a]*(Cz*f3);
	    I[45] += C[0]*W[a]*(Cx*f3);
	    I[46] += C[0]*W[a]*(Cx*f21);
	    I[47] += C[0]*W[a]*(Cy*f21);
	    I[48] += C[0]*W[a]*(Dx*f8);
	    I[49] += C[0]*W[a]*(Dx*f4);
	    I[50] += C[0]*W[a]*(Dz*f4);
	    I[51] += C[0]*W[a]*(Dz*f17);
	    I[52] += C[0]*W[a]*(Dy*f17);
	    I[53] += C[0]*W[a]*(Dy*f8);
	    I[54] += C[0]*W[a]*(Qx*f2);
	    I[55] += C[0]*W[a]*(Qx*f1);
	    I[56] += C[0]*W[a]*(Qz*f1);
	    I[57] += C[0]*W[a]*(Qz*f12);
	    I[58] += C[0]*W[a]*(Qy*f12);
	    I[59] += C[0]*W[a]*(Qy*f2);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[60]) {
	double T[60];
	for (int i = 0; i < 60; ++i) {
	    T[i] = I[i];
	}
	I[0] = T[0];
	I[11] = T[1];
	I[22] = T[2];
	I[32] = T[3];
	I[36] = T[4];
	I[34] = T[5];
	I[37] = T[6];
	I[38] = T[7];
	I[39] = T[8];
	I[19] = T[9];
	I[29] = T[10];
	I[41] = T[11];
	I[43] = T[12];
	I[45] = T[13];
	I[46] = T[14];
	I[48] = T[15];
	I[49] = T[16];
	I[21] = T[17];
	I[23] = T[18];
	I[53] = T[19];
	I[54] = T[20];
	I[50] = T[21];
	I[55] = T[22];
	I[57] = T[23];
	I[59] = T[24];
	I[25] = T[25];
	I[20] = T[26];
	I[10] = T[27];
	I[17] = T[28];
	I[8] = T[29];
	I[1] = T[30];
	I[9] = T[31];
	I[14] = T[32];
	I[12] = T[33];
	I[2] = T[34];
	I[6] = T[35];
	I[5] = T[36];
	I[7] = T[37];
	I[18] = T[38];
	I[13] = T[39];
	I[24] = T[40];
	I[26] = T[41];
	I[3] = T[42];
	I[4] = T[43];
	I[16] = T[44];
	I[15] = T[45];
	I[27] = T[46];
	I[28] = T[47];
	I[42] = T[48];
	I[31] = T[49];
	I[51] = T[50];
	I[40] = T[51];
	I[30] = T[52];
	I[52] = T[53];
	I[47] = T[54];
	I[35] = T[55];
	I[56] = T[56];
	I[44] = T[57];
	I[33] = T[58];
	I[58] = T[59];
    }

};

template<>
struct impl<meta::braket<rysq::SP, rysq::SP, rysq::SP, rysq::SP> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[16], double (&I)[256]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];

	double xkl = rkl[0];
	double ykl = rkl[1];
	double zkl = rkl[2];

// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Kx = (xkl + Dx);
	    double Ky = (ykl + Dy);
	    double Kz = (Dz + zkl);
	    double Px = (B10 + pow(Cx,2));
	    double Qy = (Cy*Dy + B00);
	    double Ry = (B01 + pow(Dy,2));
	    double f1 = Cx*Dx*xij;
	    double f10 = Cz*Dz*zij;
	    double f11 = B00*ykl;
	    double f12 = Cy*Dy;
	    double f13 = (Dz*Kz + B01);
	    double f17 = 2*B00*Dx;
	    double f2 = B00*xij;
	    double f20 = Cz*Dz;
	    double f21 = Cz*pow(Dz,2);
	    double f26 = Cz*zkl;
	    double f3 = B00*zkl;
	    double f33 = Dy*yij;
	    double f34 = B10*Dx;
	    double f35 = B01*Cy;
	    double f38 = Cx*xkl;
	    double f4 = Cx*Dx;
	    double f41 = 2*B00*Cy;
	    double f44 = 2*B00*Dy;
	    double f45 = B10*Dy;
	    double f46 = Dz*pow(Cz,2);
	    double f48 = (B01 + Dx*Kx);
	    double f49 = Dz*zij;
	    double f5 = B01*B10;
	    double f51 = Cy*ykl;
	    double f54 = 2*B00*Dz;
	    double f55 = B10*Dz;
	    double f56 = (B01 + Dy*Ky);
	    double f58 = B01*Cx;
	    double f59 = Dy*pow(Cy,2);
	    double f6 = 2*B00*Cz;
	    double f63 = Cx*pow(Dx,2);
	    double f64 = B00*zij;
	    double f67 = 2*B00*Cx;
	    double f68 = Cy*pow(Dy,2);
	    double f7 = Dx*pow(Cx,2);
	    double f76 = Dx*xij;
	    double f77 = B01*Cz;

	    I[0] += C[15]*W[a]*((f58*(Cx + xij) + f1*xkl + f5 + f63*xij + Dx*(2*f2 + Dx*Px) + xkl*(f34 + f2 + f7) + 2*B00*(f38 + B00 + 2*f4)));
	    I[1] += C[15]*W[a]*(Cy*(xij*pow(Dx,2) + f17 + f63 + B01*Ix + xkl*(f76 + B00 + f4)));
	    I[2] += C[15]*W[a]*(Cz*(xij*pow(Dx,2) + f17 + f63 + B01*Ix + xkl*(f76 + B00 + f4)));
	    I[3] += C[14]*W[a]*((xij*pow(Dx,2) + f17 + f63 + B01*Ix + xkl*(f76 + B00 + f4)));
	    I[4] += C[15]*W[a]*(Iy*(f17 + f63 + f58 + xkl*(B00 + f4)));
	    I[5] += C[15]*W[a]*(Iz*(f17 + f63 + f58 + xkl*(B00 + f4)));
	    I[6] += C[13]*W[a]*((f17 + f63 + f58 + xkl*(B00 + f4)));
	    I[7] += C[15]*W[a]*((Cy*(f68 + 2*f11) + 2*B00*(f33 + B00 + 2*f12) + yij*(f11 + f68 + f35) + B01*pow(Cy,2) + f12*yij*ykl + f59*ykl + f5 + f45*(ykl + Dy)));
	    I[8] += C[15]*W[a]*(Cz*(f11 + f68 + f44 + f35 + Ry*yij + ykl*(f12 + f33)));
	    I[9] += C[15]*W[a]*(Cx*(f11 + f68 + f44 + f35 + Ry*yij + ykl*(f12 + f33)));
	    I[10] += C[14]*W[a]*((f11 + f68 + f44 + f35 + Ry*yij + ykl*(f12 + f33)));
	    I[11] += C[15]*W[a]*(Iz*(f11 + f68 + f44 + f35 + f12*ykl));
	    I[12] += C[15]*W[a]*(Ix*(f11 + f68 + f44 + f35 + f12*ykl));
	    I[13] += C[13]*W[a]*((f11 + f68 + f44 + f35 + f12*ykl));
	    I[14] += C[15]*W[a]*(Dz*(f59 + B00*(yij + 2*Cy) + B10*Ky + yij*(f12 + f51) + ykl*pow(Cy,2)));
	    I[15] += C[15]*W[a]*(Dx*(f59 + B00*(yij + 2*Cy) + B10*Ky + yij*(f12 + f51) + ykl*pow(Cy,2)));
	    I[16] += C[11]*W[a]*((f59 + B00*(yij + 2*Cy) + B10*Ky + yij*(f12 + f51) + ykl*pow(Cy,2)));
	    I[17] += C[15]*W[a]*(Kz*(f41 + f45 + f59 + yij*(f12 + B00)));
	    I[18] += C[15]*W[a]*(Kx*(f41 + f45 + f59 + yij*(f12 + B00)));
	    I[19] += C[7]*W[a]*((f41 + f45 + f59 + yij*(f12 + B00)));
	    I[20] += C[15]*W[a]*((2*B00*(f49 + B00) + zkl*(f10 + f46) + f3*zij + f21*zij + f5 + 2*Cz*f3 + pow(Cz,2)*pow(Dz,2) + 2*Dz*f6 + f77*(Cz + zij) + f55*(Dz + zkl)));
	    I[21] += C[15]*W[a]*(Cy*(f54 + f21 + zij*pow(Dz,2) + f3 + zkl*(f49 + f20) + B01*Iz));
	    I[22] += C[15]*W[a]*(Cx*(f54 + f21 + zij*pow(Dz,2) + f3 + zkl*(f49 + f20) + B01*Iz));
	    I[23] += C[14]*W[a]*((f54 + f21 + zij*pow(Dz,2) + f3 + zkl*(f49 + f20) + B01*Iz));
	    I[24] += C[15]*W[a]*(Iy*(f77 + f54 + f21 + f3 + f20*zkl));
	    I[25] += C[15]*W[a]*(Ix*(f77 + f54 + f21 + f3 + f20*zkl));
	    I[26] += C[13]*W[a]*((f77 + f54 + f21 + f3 + f20*zkl));
	    I[27] += C[15]*W[a]*(Dy*(zkl*pow(Cz,2) + f10 + f46 + zij*(f26 + B00) + B10*Kz + f6));
	    I[28] += C[15]*W[a]*(Dx*(zkl*pow(Cz,2) + f10 + f46 + zij*(f26 + B00) + B10*Kz + f6));
	    I[29] += C[11]*W[a]*((zkl*pow(Cz,2) + f10 + f46 + zij*(f26 + B00) + B10*Kz + f6));
	    I[30] += C[15]*W[a]*(Ky*(f10 + f64 + f46 + f55 + f6));
	    I[31] += C[15]*W[a]*(Kx*(f10 + f64 + f46 + f55 + f6));
	    I[32] += C[7]*W[a]*((f10 + f64 + f46 + f55 + f6));
	    I[33] += C[15]*W[a]*(Dz*(f67 + f34 + Px*xkl + f1 + f2 + f7 + f38*xij));
	    I[34] += C[15]*W[a]*(Dy*(f67 + f34 + Px*xkl + f1 + f2 + f7 + f38*xij));
	    I[35] += C[11]*W[a]*((f67 + f34 + Px*xkl + f1 + f2 + f7 + f38*xij));
	    I[36] += C[15]*W[a]*(Ky*(f67 + f34 + f1 + f2 + f7));
	    I[37] += C[15]*W[a]*(Kz*(f67 + f34 + f1 + f2 + f7));
	    I[38] += C[7]*W[a]*((f67 + f34 + f1 + f2 + f7));
	    I[39] += C[15]*W[a]*(Qy*(f49 + f20 + f26 + zij*zkl + B00));
	    I[40] += C[15]*W[a]*((B00 + f4)*(f49 + f20 + f26 + zij*zkl + B00));
	    I[41] += C[14]*W[a]*(Dx*(f49 + f20 + f26 + zij*zkl + B00));
	    I[42] += C[15]*W[a]*(Cy*Dx*(f49 + f20 + f26 + zij*zkl + B00));
	    I[43] += C[11]*W[a]*(Cy*(f49 + f20 + f26 + zij*zkl + B00));
	    I[44] += C[14]*W[a]*(Dy*(f49 + f20 + f26 + zij*zkl + B00));
	    I[45] += C[15]*W[a]*(Cx*Dy*(f49 + f20 + f26 + zij*zkl + B00));
	    I[46] += C[11]*W[a]*(Cx*(f49 + f20 + f26 + zij*zkl + B00));
	    I[47] += C[10]*W[a]*((f49 + f20 + f26 + zij*zkl + B00));
	    I[48] += C[15]*W[a]*((f20 + B00)*(f76 + f38 + xij*xkl + B00 + f4));
	    I[49] += C[15]*W[a]*(Qy*(f76 + f38 + xij*xkl + B00 + f4));
	    I[50] += C[14]*W[a]*(Dy*(f76 + f38 + xij*xkl + B00 + f4));
	    I[51] += C[15]*W[a]*(Cz*Dy*(f76 + f38 + xij*xkl + B00 + f4));
	    I[52] += C[11]*W[a]*(Cz*(f76 + f38 + xij*xkl + B00 + f4));
	    I[53] += C[15]*W[a]*(Cy*Dz*(f76 + f38 + xij*xkl + B00 + f4));
	    I[54] += C[11]*W[a]*(Cy*(f76 + f38 + xij*xkl + B00 + f4));
	    I[55] += C[10]*W[a]*((f76 + f38 + xij*xkl + B00 + f4));
	    I[56] += C[14]*W[a]*(Dz*(f76 + f38 + xij*xkl + B00 + f4));
	    I[57] += C[15]*W[a]*(Dz*Iy*(f38 + B00 + f4));
	    I[58] += C[15]*W[a]*((f20 + B00)*(f51 + f33 + yij*ykl + Qy));
	    I[59] += C[15]*W[a]*((B00 + f4)*(f51 + f33 + yij*ykl + Qy));
	    I[60] += C[14]*W[a]*(Dx*(f51 + f33 + yij*ykl + Qy));
	    I[61] += C[15]*W[a]*(Cz*Dx*(f51 + f33 + yij*ykl + Qy));
	    I[62] += C[11]*W[a]*(Cz*(f51 + f33 + yij*ykl + Qy));
	    I[63] += C[14]*W[a]*(Dz*(f51 + f33 + yij*ykl + Qy));
	    I[64] += C[15]*W[a]*(Cx*Dz*(f51 + f33 + yij*ykl + Qy));
	    I[65] += C[11]*W[a]*(Cx*(f51 + f33 + yij*ykl + Qy));
	    I[66] += C[10]*W[a]*((f51 + f33 + yij*ykl + Qy));
	    I[67] += C[15]*W[a]*((f51 + Qy)*(f49 + f20 + B00));
	    I[68] += C[15]*W[a]*((f38 + B00 + f4)*(f49 + f20 + B00));
	    I[69] += C[15]*W[a]*((f20 + f26 + B00)*(f76 + B00 + f4));
	    I[70] += C[15]*W[a]*((f33 + Qy)*(f20 + f26 + B00));
	    I[71] += C[15]*W[a]*((f33 + Qy)*(f38 + B00 + f4));
	    I[72] += C[15]*W[a]*(Dy*Iz*(f38 + B00 + f4));
	    I[73] += C[13]*W[a]*(Dy*(f38 + B00 + f4));
	    I[74] += C[13]*W[a]*(Dz*(f38 + B00 + f4));
	    I[75] += C[15]*W[a]*(Cz*Ix*f56);
	    I[76] += C[13]*W[a]*(Cy*Dz*Kx);
	    I[77] += C[14]*W[a]*(Dz*Iy*Kx);
	    I[78] += C[15]*W[a]*(Iy*Kx*(f20 + B00));
	    I[79] += C[13]*W[a]*(Kx*(f20 + B00));
	    I[80] += C[14]*W[a]*(Kx*(f49 + f20 + B00));
	    I[81] += C[15]*W[a]*(Cy*Kx*(f49 + f20 + B00));
	    I[82] += C[7]*W[a]*(Cy*(f49 + f20 + B00));
	    I[83] += C[14]*W[a]*(Ky*(f49 + f20 + B00));
	    I[84] += C[15]*W[a]*(Cx*Ky*(f49 + f20 + B00));
	    I[85] += C[7]*W[a]*(Cx*(f49 + f20 + B00));
	    I[86] += C[6]*W[a]*((f49 + f20 + B00));
	    I[87] += C[13]*W[a]*(Dx*(f20 + f26 + B00));
	    I[88] += C[15]*W[a]*(Dx*Iy*(f20 + f26 + B00));
	    I[89] += C[11]*W[a]*(Iy*(f20 + f26 + B00));
	    I[90] += C[13]*W[a]*(Dy*(f20 + f26 + B00));
	    I[91] += C[15]*W[a]*(Dy*Ix*(f20 + f26 + B00));
	    I[92] += C[11]*W[a]*(Ix*(f20 + f26 + B00));
	    I[93] += C[7]*W[a]*(Ix*(f20 + B00));
	    I[94] += C[15]*W[a]*(Ix*Ky*(f20 + B00));
	    I[95] += C[13]*W[a]*(Ky*(f20 + B00));
	    I[96] += C[9]*W[a]*((f20 + f26 + B00));
	    I[97] += C[7]*W[a]*(Iy*(f20 + B00));
	    I[98] += C[15]*W[a]*(Cz*Iy*f48);
	    I[99] += C[15]*W[a]*(f48*(B10 + Cz*Iz));
	    I[100] += C[11]*W[a]*(Kx*(B10 + Cz*Iz));
	    I[101] += C[15]*W[a]*(Dy*Kx*(B10 + Cz*Iz));
	    I[102] += C[7]*W[a]*(Dy*(B10 + Cz*Iz));
	    I[103] += C[11]*W[a]*(Ky*(B10 + Cz*Iz));
	    I[104] += C[15]*W[a]*(Dx*Ky*(B10 + Cz*Iz));
	    I[105] += C[7]*W[a]*(Dx*(B10 + Cz*Iz));
	    I[106] += C[3]*W[a]*((B10 + Cz*Iz));
	    I[107] += C[15]*W[a]*(f56*(B10 + Cz*Iz));
	    I[108] += C[15]*W[a]*(f56*(Cx*Ix + B10));
	    I[109] += C[11]*W[a]*(Ky*(Cx*Ix + B10));
	    I[110] += C[15]*W[a]*(Dz*Ky*(Cx*Ix + B10));
	    I[111] += C[7]*W[a]*(Dz*(Cx*Ix + B10));
	    I[112] += C[11]*W[a]*(Kz*(Cx*Ix + B10));
	    I[113] += C[15]*W[a]*(Dy*Kz*(Cx*Ix + B10));
	    I[114] += C[7]*W[a]*(Dy*(Cx*Ix + B10));
	    I[115] += C[3]*W[a]*((Cx*Ix + B10));
	    I[116] += C[15]*W[a]*(f13*(Cx*Ix + B10));
	    I[117] += C[15]*W[a]*(f13*(Cy*Iy + B10));
	    I[118] += C[11]*W[a]*(Kz*(Cy*Iy + B10));
	    I[119] += C[15]*W[a]*(f48*(Cy*Iy + B10));
	    I[120] += C[11]*W[a]*(Kx*(Cy*Iy + B10));
	    I[121] += C[15]*W[a]*(Dz*Kx*(Cy*Iy + B10));
	    I[122] += C[7]*W[a]*(Dz*(Cy*Iy + B10));
	    I[123] += C[3]*W[a]*((Cy*Iy + B10));
	    I[124] += C[7]*W[a]*(Dx*(Cy*Iy + B10));
	    I[125] += C[15]*W[a]*(Dx*Kz*(Cy*Iy + B10));
	    I[126] += C[13]*W[a]*(Cy*Dx*Kz);
	    I[127] += C[15]*W[a]*(Cy*Ix*f13);
	    I[128] += C[15]*W[a]*(Dz*Ix*(f51 + Qy));
	    I[129] += C[13]*W[a]*(Dz*(f51 + Qy));
	    I[130] += C[15]*W[a]*(Dx*Iz*(f51 + Qy));
	    I[131] += C[13]*W[a]*(Dx*(f51 + Qy));
	    I[132] += C[15]*W[a]*((f51 + Qy)*(f76 + B00 + f4));
	    I[133] += C[7]*W[a]*(Cz*(f76 + B00 + f4));
	    I[134] += C[15]*W[a]*(Cz*Ky*(f76 + B00 + f4));
	    I[135] += C[14]*W[a]*(Ky*(f76 + B00 + f4));
	    I[136] += C[13]*W[a]*(Ky*(B00 + f4));
	    I[137] += C[15]*W[a]*(Iz*Ky*(B00 + f4));
	    I[138] += C[7]*W[a]*(Iz*(B00 + f4));
	    I[139] += C[11]*W[a]*(Iz*(f38 + B00 + f4));
	    I[140] += C[9]*W[a]*((f38 + B00 + f4));
	    I[141] += C[11]*W[a]*(Iy*(f38 + B00 + f4));
	    I[142] += C[7]*W[a]*(Iy*(B00 + f4));
	    I[143] += C[15]*W[a]*(Iy*Kz*(B00 + f4));
	    I[144] += C[13]*W[a]*(Kz*(B00 + f4));
	    I[145] += C[14]*W[a]*(Kz*(f76 + B00 + f4));
	    I[146] += C[15]*W[a]*(Cy*Kz*(f76 + B00 + f4));
	    I[147] += C[7]*W[a]*(Cy*(f76 + B00 + f4));
	    I[148] += C[6]*W[a]*((f76 + B00 + f4));
	    I[149] += C[5]*W[a]*((B00 + f4));
	    I[150] += C[5]*W[a]*((f20 + B00));
	    I[151] += C[15]*W[a]*(Cz*Kx*(f33 + Qy));
	    I[152] += C[7]*W[a]*(Cz*(f33 + Qy));
	    I[153] += C[14]*W[a]*(Kz*(f33 + Qy));
	    I[154] += C[15]*W[a]*(Cx*Kz*(f33 + Qy));
	    I[155] += C[7]*W[a]*(Cx*(f33 + Qy));
	    I[156] += C[6]*W[a]*((f33 + Qy));
	    I[157] += C[14]*W[a]*(Kx*(f33 + Qy));
	    I[158] += C[15]*W[a]*(Iz*Kx*Qy);
	    I[159] += C[11]*W[a]*(Cy*Iz*Kx);
	    I[160] += C[9]*W[a]*(Cy*Kx);
	    I[161] += C[10]*W[a]*(Iy*Kx);
	    I[162] += C[10]*W[a]*(Iz*Kx);
	    I[163] += C[15]*W[a]*(Cy*Iz*f48);
	    I[164] += C[13]*W[a]*(Cy*f48);
	    I[165] += C[13]*W[a]*(Cz*f48);
	    I[166] += C[14]*W[a]*(Iy*f48);
	    I[167] += C[14]*W[a]*(Iz*f48);
	    I[168] += C[6]*W[a]*(Dy*Iz);
	    I[169] += C[14]*W[a]*(Dy*Iz*Kx);
	    I[170] += C[12]*W[a]*(Dy*Kx);
	    I[171] += C[13]*W[a]*(Kx*Qy);
	    I[172] += C[12]*W[a]*(Dz*Kx);
	    I[173] += C[6]*W[a]*(Dz*Iy);
	    I[174] += C[11]*W[a]*(Cx*Iy*Kz);
	    I[175] += C[14]*W[a]*(Dx*Iy*Kz);
	    I[176] += C[6]*W[a]*(Dx*Iy);
	    I[177] += C[6]*W[a]*(Dx*Iz);
	    I[178] += C[14]*W[a]*(Dx*Iz*Ky);
	    I[179] += C[10]*W[a]*(Iz*Ky);
	    I[180] += C[12]*W[a]*(Dx*Ky);
	    I[181] += C[5]*W[a]*(Cy*Dx);
	    I[182] += C[7]*W[a]*(Cy*Dx*Iz);
	    I[183] += C[3]*W[a]*(Cy*Iz);
	    I[184] += C[15]*W[a]*(Cx*Iz*f56);
	    I[185] += C[13]*W[a]*(Cx*f56);
	    I[186] += C[13]*W[a]*(Cz*f56);
	    I[187] += C[14]*W[a]*(Ix*f56);
	    I[188] += C[14]*W[a]*(Iz*f56);
	    I[189] += C[7]*W[a]*(Iz*Qy);
	    I[190] += C[11]*W[a]*(Iz*(f51 + Qy));
	    I[191] += C[9]*W[a]*((f51 + Qy));
	    I[192] += C[11]*W[a]*(Ix*(f51 + Qy));
	    I[193] += C[7]*W[a]*(Ix*Qy);
	    I[194] += C[15]*W[a]*(Ix*Kz*Qy);
	    I[195] += C[14]*W[a]*(Dy*Ix*Kz);
	    I[196] += C[6]*W[a]*(Dy*Ix);
	    I[197] += C[10]*W[a]*(Ix*Ky);
	    I[198] += C[11]*W[a]*(Cz*Ix*Ky);
	    I[199] += C[9]*W[a]*(Cz*Ky);
	    I[200] += C[13]*W[a]*(Cz*Dx*Ky);
	    I[201] += C[5]*W[a]*(Cz*Dx);
	    I[202] += C[7]*W[a]*(Cz*Dx*Iy);
	    I[203] += C[3]*W[a]*(Cz*Iy);
	    I[204] += C[11]*W[a]*(Cz*Iy*Kx);
	    I[205] += C[9]*W[a]*(Cz*Kx);
	    I[206] += C[13]*W[a]*(Cz*Dy*Kx);
	    I[207] += C[5]*W[a]*(Cz*Dy);
	    I[208] += C[7]*W[a]*(Cz*Dy*Ix);
	    I[209] += C[3]*W[a]*(Cz*Ix);
	    I[210] += C[6]*W[a]*(Dz*Ix);
	    I[211] += C[14]*W[a]*(Dz*Ix*Ky);
	    I[212] += C[12]*W[a]*(Dz*Ky);
	    I[213] += C[5]*W[a]*(Cy*Dz);
	    I[214] += C[7]*W[a]*(Cy*Dz*Ix);
	    I[215] += C[3]*W[a]*(Cy*Ix);
	    I[216] += C[11]*W[a]*(Cy*Ix*Kz);
	    I[217] += C[9]*W[a]*(Cy*Kz);
	    I[218] += C[10]*W[a]*(Ix*Kz);
	    I[219] += C[10]*W[a]*(Iy*Kz);
	    I[220] += C[12]*W[a]*(Dx*Kz);
	    I[221] += C[12]*W[a]*(Dy*Kz);
	    I[222] += C[13]*W[a]*(Kz*Qy);
	    I[223] += C[9]*W[a]*(Cx*Kz);
	    I[224] += C[13]*W[a]*(Cx*Dy*Kz);
	    I[225] += C[5]*W[a]*(Cx*Dy);
	    I[226] += C[7]*W[a]*(Cx*Dy*Iz);
	    I[227] += C[3]*W[a]*(Cx*Iz);
	    I[228] += C[11]*W[a]*(Cx*Iz*Ky);
	    I[229] += C[9]*W[a]*(Cx*Ky);
	    I[230] += C[13]*W[a]*(Cx*Dz*Ky);
	    I[231] += C[5]*W[a]*(Cx*Dz);
	    I[232] += C[7]*W[a]*(Cx*Dz*Iy);
	    I[233] += C[15]*W[a]*(Cx*Iy*f13);
	    I[234] += C[13]*W[a]*(Cx*f13);
	    I[235] += C[13]*W[a]*(Cy*f13);
	    I[236] += C[14]*W[a]*(Ix*f13);
	    I[237] += C[14]*W[a]*(Iy*f13);
	    I[238] += C[3]*W[a]*(Cx*Iy);
	    I[239] += C[1]*W[a]*(Cx);
	    I[240] += C[1]*W[a]*(Cy);
	    I[241] += C[1]*W[a]*(Cz);
	    I[242] += C[2]*W[a]*(Ix);
	    I[243] += C[2]*W[a]*(Iy);
	    I[244] += C[2]*W[a]*(Iz);
	    I[245] += C[4]*W[a]*(Dx);
	    I[246] += C[4]*W[a]*(Dy);
	    I[247] += C[5]*W[a]*(Qy);
	    I[248] += C[4]*W[a]*(Dz);
	    I[249] += C[8]*W[a]*(Kx);
	    I[250] += C[12]*W[a]*(f48);
	    I[251] += C[8]*W[a]*(Ky);
	    I[252] += C[12]*W[a]*(f56);
	    I[253] += C[8]*W[a]*(Kz);
	    I[254] += C[12]*W[a]*(f13);
	    I[255] += C[0]*W[a]*(1);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[256]) {
	double T[256];
	for (int i = 0; i < 256; ++i) {
	    T[i] = I[i];
	}
	I[85] = T[0];
	I[86] = T[1];
	I[87] = T[2];
	I[84] = T[3];
	I[89] = T[4];
	I[93] = T[5];
	I[81] = T[6];
	I[170] = T[7];
	I[171] = T[8];
	I[169] = T[9];
	I[168] = T[10];
	I[174] = T[11];
	I[166] = T[12];
	I[162] = T[13];
	I[186] = T[14];
	I[154] = T[15];
	I[138] = T[16];
	I[234] = T[17];
	I[106] = T[18];
	I[42] = T[19];
	I[255] = T[20];
	I[254] = T[21];
	I[253] = T[22];
	I[252] = T[23];
	I[251] = T[24];
	I[247] = T[25];
	I[243] = T[26];
	I[239] = T[27];
	I[223] = T[28];
	I[207] = T[29];
	I[191] = T[30];
	I[127] = T[31];
	I[63] = T[32];
	I[117] = T[33];
	I[101] = T[34];
	I[69] = T[35];
	I[149] = T[36];
	I[213] = T[37];
	I[21] = T[38];
	I[238] = T[39];
	I[221] = T[40];
	I[220] = T[41];
	I[222] = T[42];
	I[206] = T[43];
	I[236] = T[44];
	I[237] = T[45];
	I[205] = T[46];
	I[204] = T[47];
	I[119] = T[48];
	I[102] = T[49];
	I[100] = T[50];
	I[103] = T[51];
	I[71] = T[52];
	I[118] = T[53];
	I[70] = T[54];
	I[68] = T[55];
	I[116] = T[56];
	I[121] = T[57];
	I[187] = T[58];
	I[153] = T[59];
	I[152] = T[60];
	I[155] = T[61];
	I[139] = T[62];
	I[184] = T[63];
	I[185] = T[64];
	I[137] = T[65];
	I[136] = T[66];
	I[190] = T[67];
	I[125] = T[68];
	I[215] = T[69];
	I[235] = T[70];
	I[105] = T[71];
	I[109] = T[72];
	I[97] = T[73];
	I[113] = T[74];
	I[167] = T[75];
	I[114] = T[76];
	I[120] = T[77];
	I[123] = T[78];
	I[115] = T[79];
	I[124] = T[80];
	I[126] = T[81];
	I[62] = T[82];
	I[188] = T[83];
	I[189] = T[84];
	I[61] = T[85];
	I[60] = T[86];
	I[211] = T[87];
	I[219] = T[88];
	I[203] = T[89];
	I[227] = T[90];
	I[231] = T[91];
	I[199] = T[92];
	I[55] = T[93];
	I[183] = T[94];
	I[179] = T[95];
	I[195] = T[96];
	I[59] = T[97];
	I[91] = T[98];
	I[95] = T[99];
	I[79] = T[100];
	I[111] = T[101];
	I[47] = T[102];
	I[143] = T[103];
	I[159] = T[104];
	I[31] = T[105];
	I[15] = T[106];
	I[175] = T[107];
	I[165] = T[108];
	I[133] = T[109];
	I[181] = T[110];
	I[53] = T[111];
	I[197] = T[112];
	I[229] = T[113];
	I[37] = T[114];
	I[5] = T[115];
	I[245] = T[116];
	I[250] = T[117];
	I[202] = T[118];
	I[90] = T[119];
	I[74] = T[120];
	I[122] = T[121];
	I[58] = T[122];
	I[10] = T[123];
	I[26] = T[124];
	I[218] = T[125];
	I[210] = T[126];
	I[246] = T[127];
	I[182] = T[128];
	I[178] = T[129];
	I[158] = T[130];
	I[146] = T[131];
	I[150] = T[132];
	I[23] = T[133];
	I[151] = T[134];
	I[148] = T[135];
	I[145] = T[136];
	I[157] = T[137];
	I[29] = T[138];
	I[77] = T[139];
	I[65] = T[140];
	I[73] = T[141];
	I[25] = T[142];
	I[217] = T[143];
	I[209] = T[144];
	I[212] = T[145];
	I[214] = T[146];
	I[22] = T[147];
	I[20] = T[148];
	I[17] = T[149];
	I[51] = T[150];
	I[107] = T[151];
	I[43] = T[152];
	I[232] = T[153];
	I[233] = T[154];
	I[41] = T[155];
	I[40] = T[156];
	I[104] = T[157];
	I[110] = T[158];
	I[78] = T[159];
	I[66] = T[160];
	I[72] = T[161];
	I[76] = T[162];
	I[94] = T[163];
	I[82] = T[164];
	I[83] = T[165];
	I[88] = T[166];
	I[92] = T[167];
	I[44] = T[168];
	I[108] = T[169];
	I[96] = T[170];
	I[98] = T[171];
	I[112] = T[172];
	I[56] = T[173];
	I[201] = T[174];
	I[216] = T[175];
	I[24] = T[176];
	I[28] = T[177];
	I[156] = T[178];
	I[140] = T[179];
	I[144] = T[180];
	I[18] = T[181];
	I[30] = T[182];
	I[14] = T[183];
	I[173] = T[184];
	I[161] = T[185];
	I[163] = T[186];
	I[164] = T[187];
	I[172] = T[188];
	I[46] = T[189];
	I[142] = T[190];
	I[130] = T[191];
	I[134] = T[192];
	I[38] = T[193];
	I[230] = T[194];
	I[228] = T[195];
	I[36] = T[196];
	I[132] = T[197];
	I[135] = T[198];
	I[131] = T[199];
	I[147] = T[200];
	I[19] = T[201];
	I[27] = T[202];
	I[11] = T[203];
	I[75] = T[204];
	I[67] = T[205];
	I[99] = T[206];
	I[35] = T[207];
	I[39] = T[208];
	I[7] = T[209];
	I[52] = T[210];
	I[180] = T[211];
	I[176] = T[212];
	I[50] = T[213];
	I[54] = T[214];
	I[6] = T[215];
	I[198] = T[216];
	I[194] = T[217];
	I[196] = T[218];
	I[200] = T[219];
	I[208] = T[220];
	I[224] = T[221];
	I[226] = T[222];
	I[193] = T[223];
	I[225] = T[224];
	I[33] = T[225];
	I[45] = T[226];
	I[13] = T[227];
	I[141] = T[228];
	I[129] = T[229];
	I[177] = T[230];
	I[49] = T[231];
	I[57] = T[232];
	I[249] = T[233];
	I[241] = T[234];
	I[242] = T[235];
	I[244] = T[236];
	I[248] = T[237];
	I[9] = T[238];
	I[1] = T[239];
	I[2] = T[240];
	I[3] = T[241];
	I[4] = T[242];
	I[8] = T[243];
	I[12] = T[244];
	I[16] = T[245];
	I[32] = T[246];
	I[34] = T[247];
	I[48] = T[248];
	I[64] = T[249];
	I[80] = T[250];
	I[128] = T[251];
	I[160] = T[252];
	I[192] = T[253];
	I[240] = T[254];
	I[0] = T[255];
    }

};

template<>
struct impl<meta::braket<rysq::SP, rysq::P, rysq::SP, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[4], double (&I)[48]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);

	    I[0] += C[3]*W[a]*((B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10)));
	    I[1] += C[3]*W[a]*((Dz*(B10 + Cz*Iz) + B00*(2*Cz + zij)));
	    I[2] += C[3]*W[a]*((Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx)));
	    I[3] += C[3]*W[a]*(Dy*(Cx*Ix + B10));
	    I[4] += C[3]*W[a]*(Dz*(Cx*Ix + B10));
	    I[5] += C[1]*W[a]*((Cx*Ix + B10));
	    I[6] += C[3]*W[a]*(Cz*(Dx*xij + Qx));
	    I[7] += C[3]*W[a]*(Cy*(Dx*xij + Qx));
	    I[8] += C[2]*W[a]*((Dx*xij + Qx));
	    I[9] += C[3]*W[a]*(Dx*(B10 + Cz*Iz));
	    I[10] += C[1]*W[a]*((B10 + Cz*Iz));
	    I[11] += C[3]*W[a]*(Dy*(B10 + Cz*Iz));
	    I[12] += C[3]*W[a]*(Cz*(Dy*yij + Qy));
	    I[13] += C[3]*W[a]*(Cx*(Dy*yij + Qy));
	    I[14] += C[2]*W[a]*((Dy*yij + Qy));
	    I[15] += C[3]*W[a]*(Dx*(Cy*Iy + B10));
	    I[16] += C[1]*W[a]*((Cy*Iy + B10));
	    I[17] += C[3]*W[a]*(Dz*(Cy*Iy + B10));
	    I[18] += C[3]*W[a]*(Cy*(Dz*zij + Qz));
	    I[19] += C[3]*W[a]*(Cx*(Dz*zij + Qz));
	    I[20] += C[2]*W[a]*((Dz*zij + Qz));
	    I[21] += C[3]*W[a]*(Ix*Qz);
	    I[22] += C[3]*W[a]*(Iz*Qx);
	    I[23] += C[3]*W[a]*(Iy*Qx);
	    I[24] += C[3]*W[a]*(Cz*Dx*Iy);
	    I[25] += C[1]*W[a]*(Cz*Iy);
	    I[26] += C[1]*W[a]*(Cx*Iy);
	    I[27] += C[3]*W[a]*(Cx*Dz*Iy);
	    I[28] += C[2]*W[a]*(Dz*Iy);
	    I[29] += C[3]*W[a]*(Iy*Qz);
	    I[30] += C[2]*W[a]*(Dx*Iy);
	    I[31] += C[2]*W[a]*(Dx*Iz);
	    I[32] += C[3]*W[a]*(Cy*Dx*Iz);
	    I[33] += C[1]*W[a]*(Cy*Iz);
	    I[34] += C[2]*W[a]*(Dy*Iz);
	    I[35] += C[3]*W[a]*(Cx*Dy*Iz);
	    I[36] += C[1]*W[a]*(Cx*Iz);
	    I[37] += C[3]*W[a]*(Iz*Qy);
	    I[38] += C[3]*W[a]*(Ix*Qy);
	    I[39] += C[2]*W[a]*(Dy*Ix);
	    I[40] += C[3]*W[a]*(Cz*Dy*Ix);
	    I[41] += C[1]*W[a]*(Cz*Ix);
	    I[42] += C[2]*W[a]*(Dz*Ix);
	    I[43] += C[3]*W[a]*(Cy*Dz*Ix);
	    I[44] += C[1]*W[a]*(Cy*Ix);
	    I[45] += C[0]*W[a]*(Ix);
	    I[46] += C[0]*W[a]*(Iy);
	    I[47] += C[0]*W[a]*(Iz);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[48]) {
	double T[48];
	for (int i = 0; i < 48; ++i) {
	    T[i] = I[i];
	}
	I[30] = T[0];
	I[47] = T[1];
	I[13] = T[2];
	I[25] = T[3];
	I[37] = T[4];
	I[1] = T[5];
	I[15] = T[6];
	I[14] = T[7];
	I[12] = T[8];
	I[23] = T[9];
	I[11] = T[10];
	I[35] = T[11];
	I[31] = T[12];
	I[29] = T[13];
	I[28] = T[14];
	I[18] = T[15];
	I[6] = T[16];
	I[42] = T[17];
	I[46] = T[18];
	I[45] = T[19];
	I[44] = T[20];
	I[39] = T[21];
	I[21] = T[22];
	I[17] = T[23];
	I[19] = T[24];
	I[7] = T[25];
	I[5] = T[26];
	I[41] = T[27];
	I[40] = T[28];
	I[43] = T[29];
	I[16] = T[30];
	I[20] = T[31];
	I[22] = T[32];
	I[10] = T[33];
	I[32] = T[34];
	I[33] = T[35];
	I[9] = T[36];
	I[34] = T[37];
	I[26] = T[38];
	I[24] = T[39];
	I[27] = T[40];
	I[3] = T[41];
	I[36] = T[42];
	I[38] = T[43];
	I[2] = T[44];
	I[0] = T[45];
	I[4] = T[46];
	I[8] = T[47];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::F, rysq::S, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<4> &t2, const Vector<4> &W,
			    const double (&C)[1], double (&I)[100]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 4;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 

// 	for (int a = 0; a < N; ++a) {
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);


	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double f0 = (3*B10 + pow(Cx,2));
	    double f10 = (3*pow(B10,2) + B10*(6*Cx*xij + pow(xij,2) + 6*pow(Cx,2)) + pow(Cx,2)*pow(Ix,2));
	    double f12 = (3*B10 + pow(Iy,2));
	    double f13 = (B10*(3*Cz + 2*zij) + Cz*pow(Iz,2));
	    double f14 = (3*pow(B10,2) + B10*(6*Cy*yij + pow(yij,2) + 6*pow(Cy,2)) + pow(Cy,2)*pow(Iy,2));
	    double f15 = (Cy*Iy + B10);
	    double f18 = (B10*(3*Cx + 2*xij) + Cx*pow(Ix,2));
	    double f19 = (Cy*pow(Iy,2) + B10*(3*Cy + 2*yij));
	    double f2 = (3*pow(B10,2)*(5*Cz + 2*zij) + B10*Cz*(3*pow(zij,2) + 10*pow(Cz,2) + 12*Cz*zij) + pow(Cz,3)*pow(Iz,2));
	    double f20 = (B10 + pow(Iy,2));
	    double f21 = (pow(Cy,2)*pow(Iy,3) + B10*Iy*(8*Cy*yij + 10*pow(Cy,2) + pow(yij,2)) + 3*pow(B10,2)*(5*Cy + 3*yij));
	    double f22 = (3*pow(B10,2)*(3*zij + 5*Cz) + B10*Iz*(8*Cz*zij + 10*pow(Cz,2) + pow(zij,2)) + pow(Cz,2)*pow(Iz,3));
	    double f23 = (3*pow(B10,2) + Cy*pow(Iy,3) + 3*B10*Iy*(yij + 2*Cy));
	    double f25 = (pow(Cx,2)*pow(Ix,3) + B10*Ix*(8*Cx*xij + pow(xij,2) + 10*pow(Cx,2)) + 3*pow(B10,2)*(5*Cx + 3*xij));
	    double f27 = (3*pow(B10,2) + 3*B10*Cz*(2*Cz + zij) + Iz*pow(Cz,3));
	    double f28 = (B10 + pow(Iz,2));
	    double f29 = (3*B10 + pow(Iz,2));
	    double f3 = (B10*Cx*(12*Cx*xij + 10*pow(Cx,2) + 3*pow(xij,2)) + 3*pow(B10,2)*(5*Cx + 2*xij) + pow(Cx,3)*pow(Ix,2));
	    double f30 = (3*pow(B10,2) + 3*B10*Ix*(xij + 2*Cx) + Cx*pow(Ix,3));
	    double f31 = (3*pow(B10,2) + Cz*pow(Iz,3) + 3*B10*Iz*(2*Cz + zij));
	    double f32 = (3*B10 + pow(Cz,2));
	    double f33 = (3*B10 + pow(Cy,2));
	    double f34 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f35 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f36 = 15*pow(B10,3);
	    double f39 = (3*pow(B10,2) + Iy*pow(Cy,3) + 3*B10*Cy*(yij + 2*Cy));
	    double f6 = (3*pow(B10,2) + pow(Cz,2)*pow(Iz,2) + B10*(6*Cz*zij + pow(zij,2) + 6*pow(Cz,2)));
	    double f7 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));
	    double f8 = (3*pow(B10,2) + 3*B10*Cx*(xij + 2*Cx) + Ix*pow(Cx,3));
	    double f9 = (3*pow(B10,2)*(5*Cy + 2*yij) + B10*Cy*(12*Cy*yij + 3*pow(yij,2) + 10*pow(Cy,2)) + pow(Cy,3)*pow(Iy,2));

	    I[0] += C[0]*W[a]*(Py*(Cz*zij + Pz)*(xij*(xij + 2*Cx) + Px));
	    I[1] += C[0]*W[a]*(f15*(Px + Cx*xij)*(Cz*zij + Pz));
	    I[2] += C[0]*W[a]*(Cy*f18*(Cz*zij + Pz));
	    I[3] += C[0]*W[a]*(Cx*f19*(Cz*zij + Pz));
	    I[4] += C[0]*W[a]*(Iy*f34*(Cz*zij + Pz));
	    I[5] += C[0]*W[a]*(Ix*f7*(Cz*zij + Pz));
	    I[6] += C[0]*W[a]*(Cz*Ix*Py*(xij*(xij + 2*Cx) + f0));
	    I[7] += C[0]*W[a]*(Cz*Ix*f32*(xij*(xij + 2*Cx) + f0));
	    I[8] += C[0]*W[a]*(Cy*Ix*Pz*(xij*(xij + 2*Cx) + f0));
	    I[9] += C[0]*W[a]*(Cy*Ix*f33*(xij*(xij + 2*Cx) + f0));
	    I[10] += C[0]*W[a]*((3*B10*Cx*Ix*(5*Cx*xij + 5*pow(Cx,2) + pow(xij,2)) + f36 + 9*pow(B10,2)*(5*Cx*xij + 5*pow(Cx,2) + pow(xij,2)) + pow(Cx,3)*pow(Ix,3)));
	    I[11] += C[0]*W[a]*((9*pow(B10,2)*(5*Cy*yij + 5*pow(Cy,2) + pow(yij,2)) + f36 + pow(Cy,3)*pow(Iy,3) + 3*B10*Cy*Iy*(5*Cy*yij + 5*pow(Cy,2) + pow(yij,2))));
	    I[12] += C[0]*W[a]*((f36 + 3*B10*Cz*Iz*(5*pow(Cz,2) + 5*Cz*zij + pow(zij,2)) + pow(Cz,3)*pow(Iz,3) + 9*pow(B10,2)*(5*pow(Cz,2) + 5*Cz*zij + pow(zij,2))));
	    I[13] += C[0]*W[a]*(f10*(Cz*zij + Pz));
	    I[14] += C[0]*W[a]*(f14*(Cz*zij + Pz));
	    I[15] += C[0]*W[a]*(Px*f20*(Cz*zij + Pz));
	    I[16] += C[0]*W[a]*(Pz*f20*(Px + Cx*xij));
	    I[17] += C[0]*W[a]*(Pz*f15*(xij*(xij + 2*Cx) + Px));
	    I[18] += C[0]*W[a]*(Iz*f7*(Px + Cx*xij));
	    I[19] += C[0]*W[a]*(Cy*Iz*f33*(xij*(xij + 2*Cx) + Px));
	    I[20] += C[0]*W[a]*(Cy*f13*(Px + Cx*xij));
	    I[21] += C[0]*W[a]*(Cy*f35*(xij*(xij + 2*Cx) + Px));
	    I[22] += C[0]*W[a]*(Iy*f35*(Px + Cx*xij));
	    I[23] += C[0]*W[a]*(Cz*Iy*f32*(xij*(xij + 2*Cx) + Px));
	    I[24] += C[0]*W[a]*(Cz*f19*(Px + Cx*xij));
	    I[25] += C[0]*W[a]*(Cz*f7*(xij*(xij + 2*Cx) + Px));
	    I[26] += C[0]*W[a]*(Py*f28*(Px + Cx*xij));
	    I[27] += C[0]*W[a]*(Cx*Iy*f0*f28);
	    I[28] += C[0]*W[a]*(Cx*f20*f35);
	    I[29] += C[0]*W[a]*(Ix*Pz*f19);
	    I[30] += C[0]*W[a]*(Cy*Ix*f28*f33);
	    I[31] += C[0]*W[a]*(Cy*Ix*f6);
	    I[32] += C[0]*W[a]*(Cy*f28*f34);
	    I[33] += C[0]*W[a]*(Px*f15*f28);
	    I[34] += C[0]*W[a]*(Cx*f28*f7);
	    I[35] += C[0]*W[a]*(Cx*f13*f15);
	    I[36] += C[0]*W[a]*(Iy*Iz*f8);
	    I[37] += C[0]*W[a]*(Ix*Iz*f39);
	    I[38] += C[0]*W[a]*(Ix*Iy*f27);
	    I[39] += C[0]*W[a]*(Cz*Ix*f14);
	    I[40] += C[0]*W[a]*(Cz*Ix*f20*f32);
	    I[41] += C[0]*W[a]*(Cz*f20*f34);
	    I[42] += C[0]*W[a]*(Cz*f15*f18);
	    I[43] += C[0]*W[a]*(Iz*f15*f34);
	    I[44] += C[0]*W[a]*(Ix*f15*f35);
	    I[45] += C[0]*W[a]*(Ix*Py*f13);
	    I[46] += C[0]*W[a]*(Iy*Px*f13);
	    I[47] += C[0]*W[a]*(Cz*Iy*Px*f12);
	    I[48] += C[0]*W[a]*(Cz*Iy*f12*f32);
	    I[49] += C[0]*W[a]*(Cz*Iy*f10);
	    I[50] += C[0]*W[a]*(Iy*Pz*f18);
	    I[51] += C[0]*W[a]*(Cx*Iy*Pz*f12);
	    I[52] += C[0]*W[a]*(Cx*Iy*f0*f12);
	    I[53] += C[0]*W[a]*(Cx*Iy*f6);
	    I[54] += C[0]*W[a]*(f6*(Px + Cx*xij));
	    I[55] += C[0]*W[a]*(f27*(xij*(xij + 2*Cx) + Px));
	    I[56] += C[0]*W[a]*(f39*(xij*(xij + 2*Cx) + Px));
	    I[57] += C[0]*W[a]*(f14*(Px + Cx*xij));
	    I[58] += C[0]*W[a]*(Cx*Iz*f14);
	    I[59] += C[0]*W[a]*(Cx*Iz*f0*f20);
	    I[60] += C[0]*W[a]*(Cx*Iz*f0*f29);
	    I[61] += C[0]*W[a]*(Cx*Iz*Py*f29);
	    I[62] += C[0]*W[a]*(Iz*Py*f18);
	    I[63] += C[0]*W[a]*(Iz*Px*f19);
	    I[64] += C[0]*W[a]*(Cy*Iz*Px*f29);
	    I[65] += C[0]*W[a]*(Cy*Iz*f29*f33);
	    I[66] += C[0]*W[a]*(Cy*Iz*f10);
	    I[67] += C[0]*W[a]*(Cx*Cy*f31);
	    I[68] += C[0]*W[a]*(Cx*Cz*f23);
	    I[69] += C[0]*W[a]*(Cy*Cz*f30);
	    I[70] += C[0]*W[a]*(Py*f30);
	    I[71] += C[0]*W[a]*(Pz*f30);
	    I[72] += C[0]*W[a]*(Pz*f23);
	    I[73] += C[0]*W[a]*(Px*f23);
	    I[74] += C[0]*W[a]*(Px*f31);
	    I[75] += C[0]*W[a]*(Py*f31);
	    I[76] += C[0]*W[a]*(Cy*f25);
	    I[77] += C[0]*W[a]*(Cz*f25);
	    I[78] += C[0]*W[a]*(Cz*f21);
	    I[79] += C[0]*W[a]*(Cx*f21);
	    I[80] += C[0]*W[a]*(Cx*f22);
	    I[81] += C[0]*W[a]*(Cy*f22);
	    I[82] += C[0]*W[a]*(f20*f27);
	    I[83] += C[0]*W[a]*(f20*f8);
	    I[84] += C[0]*W[a]*(f28*f8);
	    I[85] += C[0]*W[a]*(f28*f39);
	    I[86] += C[0]*W[a]*(Ix*f2);
	    I[87] += C[0]*W[a]*(Ix*f9);
	    I[88] += C[0]*W[a]*(Iz*f9);
	    I[89] += C[0]*W[a]*(Iz*f3);
	    I[90] += C[0]*W[a]*(Iy*f3);
	    I[91] += C[0]*W[a]*(Iy*f2);
	    I[92] += C[0]*W[a]*(f13*f34);
	    I[93] += C[0]*W[a]*(f19*f34);
	    I[94] += C[0]*W[a]*(f19*f35);
	    I[95] += C[0]*W[a]*(f18*f35);
	    I[96] += C[0]*W[a]*(f18*f7);
	    I[97] += C[0]*W[a]*(f13*f7);
	    I[98] += C[0]*W[a]*(f10*f15);
	    I[99] += C[0]*W[a]*(f15*f6);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[100]) {
	double T[100];
	for (int i = 0; i < 100; ++i) {
	    T[i] = I[i];
	}
	I[46] = T[0];
	I[99] = T[1];
	I[49] = T[2];
	I[69] = T[3];
	I[94] = T[4];
	I[96] = T[5];
	I[6] = T[6];
	I[2] = T[7];
	I[8] = T[8];
	I[1] = T[9];
	I[0] = T[10];
	I[11] = T[11];
	I[22] = T[12];
	I[44] = T[13];
	I[66] = T[14];
	I[64] = T[15];
	I[57] = T[16];
	I[38] = T[17];
	I[95] = T[18];
	I[41] = T[19];
	I[79] = T[20];
	I[48] = T[21];
	I[97] = T[22];
	I[32] = T[23];
	I[59] = T[24];
	I[36] = T[25];
	I[75] = T[26];
	I[80] = T[27];
	I[67] = T[28];
	I[58] = T[29];
	I[71] = T[30];
	I[78] = T[31];
	I[73] = T[32];
	I[83] = T[33];
	I[85] = T[34];
	I[89] = T[35];
	I[90] = T[36];
	I[91] = T[37];
	I[92] = T[38];
	I[56] = T[39];
	I[52] = T[40];
	I[54] = T[41];
	I[39] = T[42];
	I[93] = T[43];
	I[98] = T[44];
	I[76] = T[45];
	I[84] = T[46];
	I[14] = T[47];
	I[12] = T[48];
	I[34] = T[49];
	I[37] = T[50];
	I[17] = T[51];
	I[10] = T[52];
	I[87] = T[53];
	I[77] = T[54];
	I[42] = T[55];
	I[31] = T[56];
	I[55] = T[57];
	I[65] = T[58];
	I[60] = T[59];
	I[20] = T[60];
	I[25] = T[61];
	I[45] = T[62];
	I[63] = T[63];
	I[23] = T[64];
	I[21] = T[65];
	I[43] = T[66];
	I[29] = T[67];
	I[19] = T[68];
	I[9] = T[69];
	I[5] = T[70];
	I[7] = T[71];
	I[18] = T[72];
	I[13] = T[73];
	I[24] = T[74];
	I[26] = T[75];
	I[3] = T[76];
	I[4] = T[77];
	I[16] = T[78];
	I[15] = T[79];
	I[27] = T[80];
	I[28] = T[81];
	I[62] = T[82];
	I[50] = T[83];
	I[70] = T[84];
	I[81] = T[85];
	I[72] = T[86];
	I[51] = T[87];
	I[61] = T[88];
	I[40] = T[89];
	I[30] = T[90];
	I[82] = T[91];
	I[74] = T[92];
	I[53] = T[93];
	I[68] = T[94];
	I[47] = T[95];
	I[35] = T[96];
	I[86] = T[97];
	I[33] = T[98];
	I[88] = T[99];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::D, rysq::P, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[1], double (&I)[108]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f1 = (Dz*Pz + 2*B00*Cz);
	    double f10 = (3*pow(B10,2) + pow(Cz,2)*pow(Iz,2) + B10*(6*Cz*zij + pow(zij,2) + 6*pow(Cz,2)));
	    double f11 = (B00*(Iz*(3*Cz + zij) + 3*B10) + Dz*(B10*(3*Cz + 2*zij) + Cz*pow(Iz,2)));
	    double f12 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));
	    double f13 = (3*pow(B10,2) + B10*(6*Cx*xij + pow(xij,2) + 6*pow(Cx,2)) + pow(Cx,2)*pow(Ix,2));
	    double f14 = (Dy*(Cy*pow(Iy,2) + B10*(3*Cy + 2*yij)) + B00*(3*B10 + Iy*(3*Cy + yij)));
	    double f2 = (Dz*(Iz*pow(Cz,2) + B10*(3*Cz + zij)) + 3*B00*Pz + 2*B00*Cz*zij);
	    double f20 = (B10*(3*Cz + 2*zij) + Cz*pow(Iz,2));
	    double f21 = (3*pow(B10,2) + B10*(6*Cy*yij + pow(yij,2) + 6*pow(Cy,2)) + pow(Cy,2)*pow(Iy,2));
	    double f22 = (Cy*Iy + B10);
	    double f23 = (Cy*pow(Iy,2) + B10*(3*Cy + 2*yij));
	    double f24 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f25 = (B00*(3*B10 + Ix*(3*Cx + xij)) + Dx*(B10*(3*Cx + 2*xij) + Cx*pow(Ix,2)));
	    double f28 = (B10 + pow(Iy,2));
	    double f29 = (B10 + pow(Iz,2));
	    double f37 = (2*B00*Cy*yij + 3*B00*Py + Dy*(B10*(3*Cy + yij) + Iy*pow(Cy,2)));
	    double f38 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f39 = (Dx*(B10*(3*Cx + xij) + Ix*pow(Cx,2)) + 2*B00*Cx*xij + 3*B00*Px);
	    double f4 = (B10*(3*Cx + 2*xij) + Cx*pow(Ix,2));
	    double f5 = (B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10));
	    double f6 = (Dx*Px + 2*B00*Cx);
	    double f9 = (Dy*Py + 2*B00*Cy);

	    I[0] += C[0]*W[a]*(Cy*Iz*(Qx*xij + f6));
	    I[1] += C[0]*W[a]*(Cy*(Cz*zij + Pz)*(Dx*xij + Qx));
	    I[2] += C[0]*W[a]*(Dy*(Px + Cx*xij)*(Cz*zij + Pz));
	    I[3] += C[0]*W[a]*(Cx*Iy*(f1 + Qz*zij));
	    I[4] += C[0]*W[a]*((Px + Cx*xij)*(f1 + Qz*zij));
	    I[5] += C[0]*W[a]*(Cy*(Px + Cx*xij)*(Dz*zij + Qz));
	    I[6] += C[0]*W[a]*(Cy*Qz*(xij*(xij + 2*Cx) + Px));
	    I[7] += C[0]*W[a]*(Iy*Qz*(Px + Cx*xij));
	    I[8] += C[0]*W[a]*(Dz*f22*(Px + Cx*xij));
	    I[9] += C[0]*W[a]*(Dz*Py*(xij*(xij + 2*Cx) + Px));
	    I[10] += C[0]*W[a]*(Ix*Py*(Dz*zij + Qz));
	    I[11] += C[0]*W[a]*(Iy*Px*(Dz*zij + Qz));
	    I[12] += C[0]*W[a]*(Cx*f22*(Dz*zij + Qz));
	    I[13] += C[0]*W[a]*((Dz*(3*pow(B10,2) + pow(Cz,2)*pow(Iz,2) + B10*(6*Cz*zij + pow(zij,2) + 6*pow(Cz,2))) + 2*B00*(2*Cz + zij)*(3*B10 + Cz*Iz)));
	    I[14] += C[0]*W[a]*(f22*(f1 + Qz*zij));
	    I[15] += C[0]*W[a]*(Py*(zij*(2*B00 + Dz*(2*Cz + zij)) + f1));
	    I[16] += C[0]*W[a]*(Px*(zij*(2*B00 + Dz*(2*Cz + zij)) + f1));
	    I[17] += C[0]*W[a]*(Cy*Ix*(f1 + Qz*zij));
	    I[18] += C[0]*W[a]*(Cx*Cy*(zij*(2*B00 + Dz*(2*Cz + zij)) + f1));
	    I[19] += C[0]*W[a]*(Cx*Cz*(yij*(2*B00 + Dy*(yij + 2*Cy)) + f9));
	    I[20] += C[0]*W[a]*((2*B00*(yij + 2*Cy)*(3*B10 + Cy*Iy) + Dy*(3*pow(B10,2) + B10*(6*Cy*yij + pow(yij,2) + 6*pow(Cy,2)) + pow(Cy,2)*pow(Iy,2))));
	    I[21] += C[0]*W[a]*(Pz*(yij*(2*B00 + Dy*(yij + 2*Cy)) + f9));
	    I[22] += C[0]*W[a]*(Px*(yij*(2*B00 + Dy*(yij + 2*Cy)) + f9));
	    I[23] += C[0]*W[a]*(Iz*Px*(Dy*yij + Qy));
	    I[24] += C[0]*W[a]*(Iz*Qy*(Px + Cx*xij));
	    I[25] += C[0]*W[a]*(Cz*Qy*(xij*(xij + 2*Cx) + Px));
	    I[26] += C[0]*W[a]*(Cz*(Px + Cx*xij)*(Dy*yij + Qy));
	    I[27] += C[0]*W[a]*(Cx*(Cz*zij + Pz)*(Dy*yij + Qy));
	    I[28] += C[0]*W[a]*(Ix*Qy*(Cz*zij + Pz));
	    I[29] += C[0]*W[a]*(Dx*f22*(Cz*zij + Pz));
	    I[30] += C[0]*W[a]*((Cz*zij + Pz)*(Qx*xij + f6));
	    I[31] += C[0]*W[a]*(Cz*Iy*(Qx*xij + f6));
	    I[32] += C[0]*W[a]*(Cy*Cz*(xij*(2*B00 + Dx*(xij + 2*Cx)) + f6));
	    I[33] += C[0]*W[a]*((2*B00*(xij + 2*Cx)*(Cx*Ix + 3*B10) + Dx*(3*pow(B10,2) + B10*(6*Cx*xij + pow(xij,2) + 6*pow(Cx,2)) + pow(Cx,2)*pow(Ix,2))));
	    I[34] += C[0]*W[a]*(f22*(Qx*xij + f6));
	    I[35] += C[0]*W[a]*(Pz*(xij*(2*B00 + Dx*(xij + 2*Cx)) + f6));
	    I[36] += C[0]*W[a]*(Py*(xij*(2*B00 + Dx*(xij + 2*Cx)) + f6));
	    I[37] += C[0]*W[a]*(Iz*Py*(Dx*xij + Qx));
	    I[38] += C[0]*W[a]*(Cz*f22*(Dx*xij + Qx));
	    I[39] += C[0]*W[a]*(Iy*Pz*(Dx*xij + Qx));
	    I[40] += C[0]*W[a]*(Iy*Qx*(Cz*zij + Pz));
	    I[41] += C[0]*W[a]*(f5*(Cz*zij + Pz));
	    I[42] += C[0]*W[a]*(Cz*Ix*f5);
	    I[43] += C[0]*W[a]*(Ix*Pz*(Dy*yij + Qy));
	    I[44] += C[0]*W[a]*(f38*(Dy*yij + Qy));
	    I[45] += C[0]*W[a]*(f24*(Dy*yij + Qy));
	    I[46] += C[0]*W[a]*(Dy*Iz*f24);
	    I[47] += C[0]*W[a]*(Ix*Iz*f9);
	    I[48] += C[0]*W[a]*(Dy*Ix*f38);
	    I[49] += C[0]*W[a]*(Cx*Dy*f20);
	    I[50] += C[0]*W[a]*(Dy*Pz*(xij*(xij + 2*Cx) + Px));
	    I[51] += C[0]*W[a]*(f1*(xij*(xij + 2*Cx) + Px));
	    I[52] += C[0]*W[a]*(f9*(xij*(xij + 2*Cx) + Px));
	    I[53] += C[0]*W[a]*(f5*(Px + Cx*xij));
	    I[54] += C[0]*W[a]*(Cx*Iz*f5);
	    I[55] += C[0]*W[a]*(Cx*Qy*f29);
	    I[56] += C[0]*W[a]*(Dy*Px*f29);
	    I[57] += C[0]*W[a]*(Cz*Dy*f4);
	    I[58] += C[0]*W[a]*(Cy*Dz*f4);
	    I[59] += C[0]*W[a]*(Dz*Px*f28);
	    I[60] += C[0]*W[a]*(Cx*Dz*f23);
	    I[61] += C[0]*W[a]*(Cx*Qz*f28);
	    I[62] += C[0]*W[a]*(f12*(Dz*zij + Qz));
	    I[63] += C[0]*W[a]*(f24*(Dz*zij + Qz));
	    I[64] += C[0]*W[a]*(Dz*Iy*f24);
	    I[65] += C[0]*W[a]*(Dz*Ix*f12);
	    I[66] += C[0]*W[a]*(Ix*Iy*f1);
	    I[67] += C[0]*W[a]*(Ix*Qz*f22);
	    I[68] += C[0]*W[a]*(Iz*Qx*f22);
	    I[69] += C[0]*W[a]*(Iy*Iz*f6);
	    I[70] += C[0]*W[a]*(Dx*Iy*f38);
	    I[71] += C[0]*W[a]*(f38*(Dx*xij + Qx));
	    I[72] += C[0]*W[a]*(f12*(Dx*xij + Qx));
	    I[73] += C[0]*W[a]*(Dx*Iz*f12);
	    I[74] += C[0]*W[a]*(Cy*Dx*f20);
	    I[75] += C[0]*W[a]*(Cy*Qx*f29);
	    I[76] += C[0]*W[a]*(Dx*Py*f29);
	    I[77] += C[0]*W[a]*(Cz*Dx*f23);
	    I[78] += C[0]*W[a]*(Cz*Qx*f28);
	    I[79] += C[0]*W[a]*(Dx*Pz*f28);
	    I[80] += C[0]*W[a]*(Dx*f21);
	    I[81] += C[0]*W[a]*(Dx*f10);
	    I[82] += C[0]*W[a]*(Dy*f10);
	    I[83] += C[0]*W[a]*(Dy*f13);
	    I[84] += C[0]*W[a]*(Dz*f13);
	    I[85] += C[0]*W[a]*(Dz*f21);
	    I[86] += C[0]*W[a]*(f29*f9);
	    I[87] += C[0]*W[a]*(f29*f6);
	    I[88] += C[0]*W[a]*(f28*f6);
	    I[89] += C[0]*W[a]*(f1*f28);
	    I[90] += C[0]*W[a]*(Qz*f4);
	    I[91] += C[0]*W[a]*(Qy*f4);
	    I[92] += C[0]*W[a]*(Qy*f20);
	    I[93] += C[0]*W[a]*(Qx*f20);
	    I[94] += C[0]*W[a]*(Qx*f23);
	    I[95] += C[0]*W[a]*(Qz*f23);
	    I[96] += C[0]*W[a]*(Cy*f25);
	    I[97] += C[0]*W[a]*(Cz*f25);
	    I[98] += C[0]*W[a]*(Cz*f14);
	    I[99] += C[0]*W[a]*(Cx*f14);
	    I[100] += C[0]*W[a]*(Cx*f11);
	    I[101] += C[0]*W[a]*(Cy*f11);
	    I[102] += C[0]*W[a]*(Iy*f39);
	    I[103] += C[0]*W[a]*(Iz*f39);
	    I[104] += C[0]*W[a]*(Iz*f37);
	    I[105] += C[0]*W[a]*(Ix*f37);
	    I[106] += C[0]*W[a]*(Ix*f2);
	    I[107] += C[0]*W[a]*(Iy*f2);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[108]) {
	double T[108];
	for (int i = 0; i < 108; ++i) {
	    T[i] = I[i];
	}
	I[27] = T[0];
	I[29] = T[1];
	I[64] = T[2];
	I[106] = T[3];
	I[100] = T[4];
	I[99] = T[5];
	I[77] = T[6];
	I[94] = T[7];
	I[93] = T[8];
	I[73] = T[9];
	I[97] = T[10];
	I[102] = T[11];
	I[105] = T[12];
	I[86] = T[13];
	I[107] = T[14];
	I[85] = T[15];
	I[84] = T[16];
	I[101] = T[17];
	I[87] = T[18];
	I[46] = T[19];
	I[43] = T[20];
	I[44] = T[21];
	I[42] = T[22];
	I[66] = T[23];
	I[63] = T[24];
	I[41] = T[25];
	I[58] = T[26];
	I[70] = T[27];
	I[65] = T[28];
	I[35] = T[29];
	I[28] = T[30];
	I[22] = T[31];
	I[5] = T[32];
	I[0] = T[33];
	I[21] = T[34];
	I[2] = T[35];
	I[1] = T[36];
	I[25] = T[37];
	I[23] = T[38];
	I[20] = T[39];
	I[34] = T[40];
	I[71] = T[41];
	I[59] = T[42];
	I[56] = T[43];
	I[68] = T[44];
	I[54] = T[45];
	I[60] = T[46];
	I[61] = T[47];
	I[62] = T[48];
	I[52] = T[49];
	I[38] = T[50];
	I[74] = T[51];
	I[37] = T[52];
	I[57] = T[53];
	I[69] = T[54];
	I[51] = T[55];
	I[48] = T[56];
	I[40] = T[57];
	I[75] = T[58];
	I[78] = T[59];
	I[81] = T[60];
	I[82] = T[61];
	I[103] = T[62];
	I[96] = T[63];
	I[90] = T[64];
	I[91] = T[65];
	I[92] = T[66];
	I[95] = T[67];
	I[33] = T[68];
	I[30] = T[69];
	I[32] = T[70];
	I[26] = T[71];
	I[19] = T[72];
	I[31] = T[73];
	I[17] = T[74];
	I[15] = T[75];
	I[13] = T[76];
	I[11] = T[77];
	I[10] = T[78];
	I[8] = T[79];
	I[7] = T[80];
	I[14] = T[81];
	I[50] = T[82];
	I[36] = T[83];
	I[72] = T[84];
	I[79] = T[85];
	I[49] = T[86];
	I[12] = T[87];
	I[6] = T[88];
	I[80] = T[89];
	I[76] = T[90];
	I[39] = T[91];
	I[53] = T[92];
	I[16] = T[93];
	I[9] = T[94];
	I[83] = T[95];
	I[3] = T[96];
	I[4] = T[97];
	I[47] = T[98];
	I[45] = T[99];
	I[88] = T[100];
	I[89] = T[101];
	I[18] = T[102];
	I[24] = T[103];
	I[67] = T[104];
	I[55] = T[105];
	I[98] = T[106];
	I[104] = T[107];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::S, rysq::SP, rysq::SP> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[4], double (&I)[96]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;


	double xkl = rkl[0];
	double ykl = rkl[1];
	double zkl = rkl[2];

// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Kx = (xkl + Dx);
	    double Ky = (ykl + Dy);
	    double Kz = (Dz + zkl);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f10 = 2*B00*Cx;
	    double f11 = (B01 + Dx*Kx);
	    double f19 = 2*B00*Cz;
	    double f2 = (Cy*Dy*Ky + B01*Cy + B00*(ykl + 2*Dy));
	    double f20 = (Dz*Kz + B01);
	    double f21 = Dx*pow(Cx,2);
	    double f22 = (B00 + Cx*Kx);
	    double f23 = (B00 + Cy*Ky);
	    double f25 = B10*Dz;
	    double f27 = 2*pow(B00,2);
	    double f29 = (B01 + Dy*Ky);
	    double f31 = B10*Dx;
	    double f33 = Dy*pow(Cy,2);
	    double f35 = (B01*Cz + B00*(2*Dz + zkl) + Cz*Dz*Kz);
	    double f4 = (B01*Cx + Cx*Dx*Kx + B00*(xkl + 2*Dx));
	    double f5 = 2*B00*Cy;
	    double f7 = B10*Dy;
	    double f8 = Dz*pow(Cz,2);

	    I[0] += C[3]*W[a]*((f27 + xkl*(f10 + f21 + f31) + B01*Px + Dx*(2*f10 + Dx*Px)));
	    I[1] += C[3]*W[a]*((f27 + Dz*(Dz*Pz + 2*f19) + B01*Pz + zkl*(f19 + f25 + f8)));
	    I[2] += C[3]*W[a]*((f27 + B01*Py + Dy*(2*f5 + Dy*Py) + ykl*(f33 + f5 + f7)));
	    I[3] += C[3]*W[a]*(Dx*(f33 + B10*Ky + ykl*pow(Cy,2) + f5));
	    I[4] += C[3]*W[a]*(Dz*(f33 + B10*Ky + ykl*pow(Cy,2) + f5));
	    I[5] += C[2]*W[a]*((f33 + B10*Ky + ykl*pow(Cy,2) + f5));
	    I[6] += C[3]*W[a]*(Kx*(f33 + f5 + f7));
	    I[7] += C[3]*W[a]*(Kz*(f33 + f5 + f7));
	    I[8] += C[1]*W[a]*((f33 + f5 + f7));
	    I[9] += C[3]*W[a]*(Dy*(zkl*pow(Cz,2) + f19 + B10*Kz + f8));
	    I[10] += C[3]*W[a]*(Dx*(zkl*pow(Cz,2) + f19 + B10*Kz + f8));
	    I[11] += C[2]*W[a]*((zkl*pow(Cz,2) + f19 + B10*Kz + f8));
	    I[12] += C[3]*W[a]*(Ky*(f19 + f25 + f8));
	    I[13] += C[3]*W[a]*(Kx*(f19 + f25 + f8));
	    I[14] += C[1]*W[a]*((f19 + f25 + f8));
	    I[15] += C[3]*W[a]*(Dz*(f10 + f21 + B10*Kx + xkl*pow(Cx,2)));
	    I[16] += C[3]*W[a]*(Dy*(f10 + f21 + B10*Kx + xkl*pow(Cx,2)));
	    I[17] += C[2]*W[a]*((f10 + f21 + B10*Kx + xkl*pow(Cx,2)));
	    I[18] += C[3]*W[a]*(Kz*(f10 + f21 + f31));
	    I[19] += C[3]*W[a]*(Ky*(f10 + f21 + f31));
	    I[20] += C[1]*W[a]*((f10 + f21 + f31));
	    I[21] += C[3]*W[a]*(Cy*Dx*(Cz*zkl + Qz));
	    I[22] += C[3]*W[a]*(Cx*Dy*(Cz*zkl + Qz));
	    I[23] += C[2]*W[a]*(Cx*(Cz*zkl + Qz));
	    I[24] += C[3]*W[a]*(Qx*(Cz*zkl + Qz));
	    I[25] += C[3]*W[a]*(Qy*(Cz*zkl + Qz));
	    I[26] += C[2]*W[a]*(Cy*(Cz*zkl + Qz));
	    I[27] += C[3]*W[a]*(Cy*Kx*Qz);
	    I[28] += C[3]*W[a]*(Px*f20);
	    I[29] += C[3]*W[a]*(Py*f20);
	    I[30] += C[3]*W[a]*(Cx*Cy*f20);
	    I[31] += C[3]*W[a]*(Cx*Dz*f23);
	    I[32] += C[3]*W[a]*(Cz*Dx*f23);
	    I[33] += C[2]*W[a]*(Cz*f23);
	    I[34] += C[3]*W[a]*(Qx*f23);
	    I[35] += C[2]*W[a]*(Cx*f23);
	    I[36] += C[3]*W[a]*(Cx*f2);
	    I[37] += C[3]*W[a]*(Cz*f2);
	    I[38] += C[3]*W[a]*(Qy*f22);
	    I[39] += C[3]*W[a]*(Qz*f22);
	    I[40] += C[3]*W[a]*(Cx*Ky*Qz);
	    I[41] += C[1]*W[a]*(Cx*Qz);
	    I[42] += C[3]*W[a]*(Qz*f23);
	    I[43] += C[1]*W[a]*(Cy*Qz);
	    I[44] += C[3]*W[a]*(Cy*f4);
	    I[45] += C[3]*W[a]*(Cz*f4);
	    I[46] += C[3]*W[a]*(Cz*Dy*f22);
	    I[47] += C[2]*W[a]*(Cz*f22);
	    I[48] += C[3]*W[a]*(Cz*Kx*Qy);
	    I[49] += C[1]*W[a]*(Cz*Qy);
	    I[50] += C[1]*W[a]*(Cx*Qy);
	    I[51] += C[3]*W[a]*(Cx*Kz*Qy);
	    I[52] += C[2]*W[a]*(Cx*Cy*Kz);
	    I[53] += C[3]*W[a]*(Dx*Kz*Py);
	    I[54] += C[2]*W[a]*(Kz*Py);
	    I[55] += C[3]*W[a]*(Cy*Kz*Qx);
	    I[56] += C[1]*W[a]*(Cy*Qx);
	    I[57] += C[2]*W[a]*(Cy*f22);
	    I[58] += C[3]*W[a]*(Cy*Dz*f22);
	    I[59] += C[1]*W[a]*(Cx*Cy*Dz);
	    I[60] += C[0]*W[a]*(Cx*Cy);
	    I[61] += C[3]*W[a]*(Cx*Cz*f29);
	    I[62] += C[3]*W[a]*(Px*f29);
	    I[63] += C[3]*W[a]*(Pz*f29);
	    I[64] += C[3]*W[a]*(Pz*f11);
	    I[65] += C[3]*W[a]*(Py*f11);
	    I[66] += C[3]*W[a]*(Cy*Cz*f11);
	    I[67] += C[2]*W[a]*(Cy*Cz*Kx);
	    I[68] += C[1]*W[a]*(Cy*Cz*Dx);
	    I[69] += C[0]*W[a]*(Cy*Cz);
	    I[70] += C[1]*W[a]*(Cz*Qx);
	    I[71] += C[3]*W[a]*(Cz*Ky*Qx);
	    I[72] += C[2]*W[a]*(Cx*Cz*Ky);
	    I[73] += C[0]*W[a]*(Cx*Cz);
	    I[74] += C[1]*W[a]*(Cx*Cz*Dy);
	    I[75] += C[1]*W[a]*(Dy*Pz);
	    I[76] += C[3]*W[a]*(Dy*Kx*Pz);
	    I[77] += C[2]*W[a]*(Kx*Pz);
	    I[78] += C[2]*W[a]*(Kx*Py);
	    I[79] += C[3]*W[a]*(Dz*Kx*Py);
	    I[80] += C[1]*W[a]*(Dz*Py);
	    I[81] += C[0]*W[a]*(Py);
	    I[82] += C[1]*W[a]*(Dx*Py);
	    I[83] += C[3]*W[a]*(Dx*Ky*Pz);
	    I[84] += C[1]*W[a]*(Dx*Pz);
	    I[85] += C[0]*W[a]*(Pz);
	    I[86] += C[2]*W[a]*(Ky*Pz);
	    I[87] += C[2]*W[a]*(Ky*Px);
	    I[88] += C[3]*W[a]*(Dz*Ky*Px);
	    I[89] += C[1]*W[a]*(Dz*Px);
	    I[90] += C[2]*W[a]*(Kz*Px);
	    I[91] += C[3]*W[a]*(Dy*Kz*Px);
	    I[92] += C[1]*W[a]*(Dy*Px);
	    I[93] += C[0]*W[a]*(Px);
	    I[94] += C[3]*W[a]*(Cx*f35);
	    I[95] += C[3]*W[a]*(Cy*f35);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[96]) {
	double T[96];
	for (int i = 0; i < 96; ++i) {
	    T[i] = I[i];
	}
	I[30] = T[0];
	I[92] = T[1];
	I[61] = T[2];
	I[55] = T[3];
	I[67] = T[4];
	I[49] = T[5];
	I[37] = T[6];
	I[85] = T[7];
	I[13] = T[8];
	I[86] = T[9];
	I[80] = T[10];
	I[74] = T[11];
	I[68] = T[12];
	I[44] = T[13];
	I[20] = T[14];
	I[42] = T[15];
	I[36] = T[16];
	I[24] = T[17];
	I[78] = T[18];
	I[54] = T[19];
	I[6] = T[20];
	I[83] = T[21];
	I[88] = T[22];
	I[76] = T[23];
	I[82] = T[24];
	I[89] = T[25];
	I[77] = T[26];
	I[47] = T[27];
	I[90] = T[28];
	I[91] = T[29];
	I[93] = T[30];
	I[69] = T[31];
	I[59] = T[32];
	I[53] = T[33];
	I[57] = T[34];
	I[51] = T[35];
	I[63] = T[36];
	I[65] = T[37];
	I[39] = T[38];
	I[46] = T[39];
	I[70] = T[40];
	I[22] = T[41];
	I[71] = T[42];
	I[23] = T[43];
	I[33] = T[44];
	I[34] = T[45];
	I[40] = T[46];
	I[28] = T[47];
	I[41] = T[48];
	I[17] = T[49];
	I[15] = T[50];
	I[87] = T[51];
	I[75] = T[52];
	I[79] = T[53];
	I[73] = T[54];
	I[81] = T[55];
	I[9] = T[56];
	I[27] = T[57];
	I[45] = T[58];
	I[21] = T[59];
	I[3] = T[60];
	I[64] = T[61];
	I[60] = T[62];
	I[62] = T[63];
	I[32] = T[64];
	I[31] = T[65];
	I[35] = T[66];
	I[29] = T[67];
	I[11] = T[68];
	I[5] = T[69];
	I[10] = T[70];
	I[58] = T[71];
	I[52] = T[72];
	I[4] = T[73];
	I[16] = T[74];
	I[14] = T[75];
	I[38] = T[76];
	I[26] = T[77];
	I[25] = T[78];
	I[43] = T[79];
	I[19] = T[80];
	I[1] = T[81];
	I[7] = T[82];
	I[56] = T[83];
	I[8] = T[84];
	I[2] = T[85];
	I[50] = T[86];
	I[48] = T[87];
	I[66] = T[88];
	I[18] = T[89];
	I[72] = T[90];
	I[84] = T[91];
	I[12] = T[92];
	I[0] = T[93];
	I[94] = T[94];
	I[95] = T[95];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::S, rysq::SP, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[2], double (&I)[24]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;



// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);

	    I[0] += C[1]*W[a]*((Dx*Px + 2*B00*Cx));
	    I[1] += C[1]*W[a]*((Dy*Py + 2*B00*Cy));
	    I[2] += C[1]*W[a]*((Dz*Pz + 2*B00*Cz));
	    I[3] += C[1]*W[a]*(Cy*Cz*Dx);
	    I[4] += C[1]*W[a]*(Cx*Cz*Dy);
	    I[5] += C[0]*W[a]*(Cx*Cz);
	    I[6] += C[0]*W[a]*(Cy*Cz);
	    I[7] += C[1]*W[a]*(Cy*Qx);
	    I[8] += C[1]*W[a]*(Cz*Qx);
	    I[9] += C[1]*W[a]*(Dx*Pz);
	    I[10] += C[1]*W[a]*(Dx*Py);
	    I[11] += C[0]*W[a]*(Py);
	    I[12] += C[1]*W[a]*(Dz*Py);
	    I[13] += C[1]*W[a]*(Dz*Px);
	    I[14] += C[1]*W[a]*(Cx*Cy*Dz);
	    I[15] += C[0]*W[a]*(Cx*Cy);
	    I[16] += C[1]*W[a]*(Cx*Qy);
	    I[17] += C[1]*W[a]*(Cz*Qy);
	    I[18] += C[0]*W[a]*(Pz);
	    I[19] += C[1]*W[a]*(Dy*Pz);
	    I[20] += C[1]*W[a]*(Dy*Px);
	    I[21] += C[0]*W[a]*(Px);
	    I[22] += C[1]*W[a]*(Cx*Qz);
	    I[23] += C[1]*W[a]*(Cy*Qz);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[24]) {
	double T[24];
	for (int i = 0; i < 24; ++i) {
	    T[i] = I[i];
	}
	I[6] = T[0];
	I[13] = T[1];
	I[20] = T[2];
	I[11] = T[3];
	I[16] = T[4];
	I[4] = T[5];
	I[5] = T[6];
	I[9] = T[7];
	I[10] = T[8];
	I[8] = T[9];
	I[7] = T[10];
	I[1] = T[11];
	I[19] = T[12];
	I[18] = T[13];
	I[21] = T[14];
	I[3] = T[15];
	I[15] = T[16];
	I[17] = T[17];
	I[2] = T[18];
	I[14] = T[19];
	I[12] = T[20];
	I[0] = T[21];
	I[22] = T[22];
	I[23] = T[23];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::SP, rysq::P, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[2], double (&I)[72]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f1 = (B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10));
	    double f10 = (Dy*Iy + B00);
	    double f11 = (Cy*Iy + B10);
	    double f12 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f14 = (Dx*Px + 2*B00*Cx);
	    double f15 = (Dx*Ix + B00);
	    double f22 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f3 = 3*B00*B10;
	    double f4 = (Dy*Py + 2*B00*Cy);
	    double f5 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));
	    double f6 = (Dz*Pz + 2*B00*Cz);
	    double f8 = (Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx));

	    I[0] += C[1]*W[a]*((Dx*Ix*pow(Cx,2) + B00*Cx*(3*Cx + 2*xij) + B10*Dx*(3*Cx + xij) + f3));
	    I[1] += C[1]*W[a]*((Dy*Iy*pow(Cy,2) + f3 + B00*Cy*(3*Cy + 2*yij) + B10*Dy*(3*Cy + yij)));
	    I[2] += C[1]*W[a]*((B00*Cz*(3*Cz + 2*zij) + Dz*Iz*pow(Cz,2) + B10*Dz*(3*Cz + zij) + f3));
	    I[3] += C[1]*W[a]*(Cx*(f6 + Qz*zij));
	    I[4] += C[1]*W[a]*(Cy*(f6 + Qz*zij));
	    I[5] += C[1]*W[a]*(Cx*Cy*(Dz*zij + Qz));
	    I[6] += C[1]*W[a]*(Cy*Dz*(Px + Cx*xij));
	    I[7] += C[1]*W[a]*(Qy*(Px + Cx*xij));
	    I[8] += C[1]*W[a]*(Qz*(Px + Cx*xij));
	    I[9] += C[1]*W[a]*(Cz*Dy*(Px + Cx*xij));
	    I[10] += C[1]*W[a]*(Cx*Dy*(Cz*zij + Pz));
	    I[11] += C[1]*W[a]*(Qy*(Cz*zij + Pz));
	    I[12] += C[1]*W[a]*(Qx*(Cz*zij + Pz));
	    I[13] += C[1]*W[a]*(Cy*Dx*(Cz*zij + Pz));
	    I[14] += C[1]*W[a]*(Dx*Iy*Pz);
	    I[15] += C[1]*W[a]*(Cx*Cz*f10);
	    I[16] += C[0]*W[a]*(Cx*Cz*Dy);
	    I[17] += C[1]*W[a]*(Dy*Iz*Px);
	    I[18] += C[1]*W[a]*(Px*(Dz*zij + Qz));
	    I[19] += C[1]*W[a]*(Py*(Dz*zij + Qz));
	    I[20] += C[1]*W[a]*(Dz*Ix*Py);
	    I[21] += C[0]*W[a]*(Cx*Cy*Dz);
	    I[22] += C[1]*W[a]*(Cx*Dz*f11);
	    I[23] += C[1]*W[a]*(Cz*Dx*f11);
	    I[24] += C[1]*W[a]*(Cy*Cz*f15);
	    I[25] += C[0]*W[a]*(Cy*Cz*Dx);
	    I[26] += C[1]*W[a]*(Dx*f5);
	    I[27] += C[1]*W[a]*(Py*f15);
	    I[28] += C[1]*W[a]*(Pz*f15);
	    I[29] += C[1]*W[a]*(Dy*Ix*Pz);
	    I[30] += C[0]*W[a]*(Dy*Pz);
	    I[31] += C[1]*W[a]*(Dy*f12);
	    I[32] += C[1]*W[a]*(Px*f10);
	    I[33] += C[1]*W[a]*(Pz*f10);
	    I[34] += C[0]*W[a]*(Dx*Pz);
	    I[35] += C[1]*W[a]*(Dx*Iz*Py);
	    I[36] += C[0]*W[a]*(Dx*Py);
	    I[37] += C[1]*W[a]*(Dx*f22);
	    I[38] += C[1]*W[a]*(Dy*f22);
	    I[39] += C[0]*W[a]*(Dy*Px);
	    I[40] += C[1]*W[a]*(Dz*Iy*Px);
	    I[41] += C[0]*W[a]*(Dz*Px);
	    I[42] += C[0]*W[a]*(Dz*Py);
	    I[43] += C[1]*W[a]*(Dz*f12);
	    I[44] += C[1]*W[a]*(Dz*f5);
	    I[45] += C[1]*W[a]*(Cx*f1);
	    I[46] += C[1]*W[a]*(Cz*f1);
	    I[47] += C[1]*W[a]*(Cz*Ix*Qy);
	    I[48] += C[0]*W[a]*(Cz*Qy);
	    I[49] += C[1]*W[a]*(Cx*Iz*Qy);
	    I[50] += C[0]*W[a]*(Cx*Qy);
	    I[51] += C[1]*W[a]*(Cx*Iy*Qz);
	    I[52] += C[0]*W[a]*(Cx*Qz);
	    I[53] += C[1]*W[a]*(Cy*Ix*Qz);
	    I[54] += C[0]*W[a]*(Cy*Qz);
	    I[55] += C[1]*W[a]*(Cy*f8);
	    I[56] += C[1]*W[a]*(Cz*f8);
	    I[57] += C[1]*W[a]*(Cz*Iy*Qx);
	    I[58] += C[0]*W[a]*(Cz*Qx);
	    I[59] += C[1]*W[a]*(Cy*Iz*Qx);
	    I[60] += C[0]*W[a]*(Cy*Qx);
	    I[61] += C[1]*W[a]*(Qx*f11);
	    I[62] += C[1]*W[a]*(Qz*f11);
	    I[63] += C[1]*W[a]*(Iy*f6);
	    I[64] += C[1]*W[a]*(Ix*f6);
	    I[65] += C[1]*W[a]*(Ix*f4);
	    I[66] += C[1]*W[a]*(Iz*f4);
	    I[67] += C[1]*W[a]*(Iz*f14);
	    I[68] += C[1]*W[a]*(Iy*f14);
	    I[69] += C[0]*W[a]*(f14);
	    I[70] += C[0]*W[a]*(f4);
	    I[71] += C[0]*W[a]*(f6);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[72]) {
	double T[72];
	for (int i = 0; i < 72; ++i) {
	    T[i] = I[i];
	}
	I[6] = T[0];
	I[37] = T[1];
	I[68] = T[2];
	I[70] = T[3];
	I[71] = T[4];
	I[69] = T[5];
	I[57] = T[6];
	I[33] = T[7];
	I[58] = T[8];
	I[34] = T[9];
	I[46] = T[10];
	I[47] = T[11];
	I[22] = T[12];
	I[23] = T[13];
	I[14] = T[14];
	I[40] = T[15];
	I[28] = T[16];
	I[42] = T[17];
	I[66] = T[18];
	I[67] = T[19];
	I[55] = T[20];
	I[51] = T[21];
	I[63] = T[22];
	I[17] = T[23];
	I[11] = T[24];
	I[5] = T[25];
	I[13] = T[26];
	I[7] = T[27];
	I[8] = T[28];
	I[32] = T[29];
	I[26] = T[30];
	I[30] = T[31];
	I[36] = T[32];
	I[38] = T[33];
	I[2] = T[34];
	I[19] = T[35];
	I[1] = T[36];
	I[20] = T[37];
	I[44] = T[38];
	I[24] = T[39];
	I[60] = T[40];
	I[48] = T[41];
	I[49] = T[42];
	I[54] = T[43];
	I[61] = T[44];
	I[39] = T[45];
	I[41] = T[46];
	I[35] = T[47];
	I[29] = T[48];
	I[45] = T[49];
	I[27] = T[50];
	I[64] = T[51];
	I[52] = T[52];
	I[59] = T[53];
	I[53] = T[54];
	I[9] = T[55];
	I[10] = T[56];
	I[16] = T[57];
	I[4] = T[58];
	I[21] = T[59];
	I[3] = T[60];
	I[15] = T[61];
	I[65] = T[62];
	I[62] = T[63];
	I[56] = T[64];
	I[31] = T[65];
	I[43] = T[66];
	I[18] = T[67];
	I[12] = T[68];
	I[0] = T[69];
	I[25] = T[70];
	I[50] = T[71];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::S, rysq::D, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[1], double (&I)[36]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;



// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double Rx = (B01 + pow(Dx,2));
	    double Ry = (B01 + pow(Dy,2));
	    double Rz = (pow(Dz,2) + B01);
	    double f0 = (Dz*Pz + 2*B00*Cz);
	    double f1 = (Dy*Py + 2*B00*Cy);
	    double f13 = 2*pow(B00,2);
	    double f15 = (2*B00*Dy + Cy*Ry);
	    double f3 = (2*B00*Dz + Cz*Rz);
	    double f7 = B01*B10;
	    double f8 = (2*B00*Dx + Cx*Rx);
	    double f9 = (Dx*Px + 2*B00*Cx);

	    I[0] += C[0]*W[a]*((4*B00*Cy*Dy + f13 + B01*pow(Cy,2) + f7 + Py*pow(Dy,2)));
	    I[1] += C[0]*W[a]*((f13 + B01*pow(Cz,2) + Pz*pow(Dz,2) + f7 + 4*B00*Cz*Dz));
	    I[2] += C[0]*W[a]*((B01*pow(Cx,2) + f13 + Px*pow(Dx,2) + 4*B00*Cx*Dx + f7));
	    I[3] += C[0]*W[a]*(Cx*Cy*Rz);
	    I[4] += C[0]*W[a]*(Cz*Dy*Qx);
	    I[5] += C[0]*W[a]*(Cz*Dx*Qy);
	    I[6] += C[0]*W[a]*(Dx*Dz*Py);
	    I[7] += C[0]*W[a]*(Cy*Dz*Qx);
	    I[8] += C[0]*W[a]*(Cy*Dx*Qz);
	    I[9] += C[0]*W[a]*(Dx*Dy*Pz);
	    I[10] += C[0]*W[a]*(Dy*Dz*Px);
	    I[11] += C[0]*W[a]*(Cx*Dz*Qy);
	    I[12] += C[0]*W[a]*(Cx*Dy*Qz);
	    I[13] += C[0]*W[a]*(Cx*Cz*Ry);
	    I[14] += C[0]*W[a]*(Cy*Cz*Rx);
	    I[15] += C[0]*W[a]*(Py*Rx);
	    I[16] += C[0]*W[a]*(Pz*Rx);
	    I[17] += C[0]*W[a]*(Pz*Ry);
	    I[18] += C[0]*W[a]*(Px*Ry);
	    I[19] += C[0]*W[a]*(Px*Rz);
	    I[20] += C[0]*W[a]*(Py*Rz);
	    I[21] += C[0]*W[a]*(Cy*f8);
	    I[22] += C[0]*W[a]*(Cz*f8);
	    I[23] += C[0]*W[a]*(Cz*f15);
	    I[24] += C[0]*W[a]*(Cx*f15);
	    I[25] += C[0]*W[a]*(Cx*f3);
	    I[26] += C[0]*W[a]*(Cy*f3);
	    I[27] += C[0]*W[a]*(Dx*f0);
	    I[28] += C[0]*W[a]*(Dx*f1);
	    I[29] += C[0]*W[a]*(Dz*f1);
	    I[30] += C[0]*W[a]*(Dz*f9);
	    I[31] += C[0]*W[a]*(Dy*f9);
	    I[32] += C[0]*W[a]*(Dy*f0);
	    I[33] += C[0]*W[a]*(Qx*Qz);
	    I[34] += C[0]*W[a]*(Qx*Qy);
	    I[35] += C[0]*W[a]*(Qy*Qz);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[36]) {
	double T[36];
	for (int i = 0; i < 36; ++i) {
	    T[i] = I[i];
	}
	I[7] = T[0];
	I[14] = T[1];
	I[0] = T[2];
	I[15] = T[3];
	I[22] = T[4];
	I[23] = T[5];
	I[25] = T[6];
	I[27] = T[7];
	I[29] = T[8];
	I[20] = T[9];
	I[30] = T[10];
	I[33] = T[11];
	I[34] = T[12];
	I[10] = T[13];
	I[5] = T[14];
	I[1] = T[15];
	I[2] = T[16];
	I[8] = T[17];
	I[6] = T[18];
	I[12] = T[19];
	I[13] = T[20];
	I[3] = T[21];
	I[4] = T[22];
	I[11] = T[23];
	I[9] = T[24];
	I[16] = T[25];
	I[17] = T[26];
	I[26] = T[27];
	I[19] = T[28];
	I[31] = T[29];
	I[24] = T[30];
	I[18] = T[31];
	I[32] = T[32];
	I[28] = T[33];
	I[21] = T[34];
	I[35] = T[35];
    }

};

template<>
struct impl<meta::braket<rysq::SP, rysq::P, rysq::S, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[2], double (&I)[12]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 

// 	for (int a = 0; a < N; ++a) {
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);


	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);

	    I[0] += C[1]*W[a]*((Cx*Ix + B10));
	    I[1] += C[1]*W[a]*((Cy*Iy + B10));
	    I[2] += C[1]*W[a]*((B10 + Cz*Iz));
	    I[3] += C[1]*W[a]*(Cz*Ix);
	    I[4] += C[1]*W[a]*(Cz*Iy);
	    I[5] += C[1]*W[a]*(Cx*Iy);
	    I[6] += C[1]*W[a]*(Cx*Iz);
	    I[7] += C[1]*W[a]*(Cy*Iz);
	    I[8] += C[1]*W[a]*(Cy*Ix);
	    I[9] += C[0]*W[a]*(Ix);
	    I[10] += C[0]*W[a]*(Iy);
	    I[11] += C[0]*W[a]*(Iz);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[12]) {
	double T[12];
	for (int i = 0; i < 12; ++i) {
	    T[i] = I[i];
	}
	I[1] = T[0];
	I[6] = T[1];
	I[11] = T[2];
	I[3] = T[3];
	I[7] = T[4];
	I[5] = T[5];
	I[9] = T[6];
	I[10] = T[7];
	I[2] = T[8];
	I[0] = T[9];
	I[4] = T[10];
	I[8] = T[11];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::SP, rysq::D, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[2], double (&I)[144]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double Rx = (B01 + pow(Dx,2));
	    double Ry = (B01 + pow(Dy,2));
	    double Rz = (pow(Dz,2) + B01);
	    double f1 = (Dz*Pz + 2*B00*Cz);
	    double f14 = (2*B00*Dx + Ix*Rx);
	    double f16 = (2*B00*Cy*yij + 3*B00*Py + Dy*(B10*(3*Cy + yij) + Iy*pow(Cy,2)));
	    double f2 = (4*B00*Cy*Dy + Py*Ry + 2*pow(B00,2));
	    double f21 = (Cy*Iy + B10);
	    double f22 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f23 = (2*B00*Dx + Cx*Rx);
	    double f24 = (B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10));
	    double f28 = (2*B00*Dy + Iy*Ry);
	    double f29 = (Dx*(B10*(3*Cx + xij) + Ix*pow(Cx,2)) + 2*B00*Cx*xij + 3*B00*Px);
	    double f3 = B01*B10;
	    double f30 = (Iz*Rz + 2*B00*Dz);
	    double f31 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));
	    double f32 = (Dy*Py + 2*B00*Cy);
	    double f33 = (2*B00*Dz + Cz*Rz);
	    double f36 = 2*pow(B00,2);
	    double f4 = (Dx*Px + 2*B00*Cx);
	    double f40 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f42 = (2*B00*Dy + Cy*Ry);
	    double f7 = (Dz*(Iz*pow(Cz,2) + B10*(3*Cz + zij)) + 3*B00*Pz + 2*B00*Cz*zij);
	    double f8 = (B00 + Dz*Iz);

	    I[0] += C[1]*W[a]*(Cz*(f36 + 2*B00*Dx*(xij + 2*Cx) + B01*Cx*Ix + pow(Dx,2)*(Cx*Ix + B10) + f3));
	    I[1] += C[1]*W[a]*(Cx*(yij*(2*B00*Dy + Cy*Ry) + f2));
	    I[2] += C[1]*W[a]*(Cz*(yij*(2*B00*Dy + Cy*Ry) + f2));
	    I[3] += C[1]*W[a]*((B01*Iy*pow(Cy,2) + 4*B00*Cy*Dy*yij + 6*B00*Dy*Py + yij*(f36 + f3 + Py*pow(Dy,2)) + Cy*(3*f36 + 3*f3 + pow(Dy,2)*(3*B10 + pow(Cy,2)))));
	    I[4] += C[1]*W[a]*(Cy*(f36 + pow(Dz,2)*(B10 + Cz*Iz) + 2*B00*Dz*(2*Cz + zij) + f3 + B01*Cz*Iz));
	    I[5] += C[1]*W[a]*(Cx*(f36 + pow(Dz,2)*(B10 + Cz*Iz) + 2*B00*Dz*(2*Cz + zij) + f3 + B01*Cz*Iz));
	    I[6] += C[1]*W[a]*((6*B00*Dz*Pz + B01*Iz*pow(Cz,2) + Cz*(3*f36 + 3*f3 + pow(Dz,2)*(3*B10 + pow(Cz,2))) + zij*(f36 + Pz*pow(Dz,2) + f3) + 4*B00*Cz*Dz*zij));
	    I[7] += C[1]*W[a]*(Iy*(f36 + B01*pow(Cz,2) + Pz*pow(Dz,2) + f3 + 4*B00*Cz*Dz));
	    I[8] += C[1]*W[a]*(Ix*(f36 + B01*pow(Cz,2) + Pz*pow(Dz,2) + f3 + 4*B00*Cz*Dz));
	    I[9] += C[0]*W[a]*((f36 + B01*pow(Cz,2) + Pz*pow(Dz,2) + f3 + 4*B00*Cz*Dz));
	    I[10] += C[1]*W[a]*(Iz*(B01*pow(Cx,2) + f36 + Px*pow(Dx,2) + 4*B00*Cx*Dx + f3));
	    I[11] += C[1]*W[a]*(Iy*(B01*pow(Cx,2) + f36 + Px*pow(Dx,2) + 4*B00*Cx*Dx + f3));
	    I[12] += C[0]*W[a]*((B01*pow(Cx,2) + f36 + Px*pow(Dx,2) + 4*B00*Cx*Dx + f3));
	    I[13] += C[1]*W[a]*((B01*Ix*pow(Cx,2) + 6*B00*Dx*Px + Cx*(3*f36 + 3*f3 + pow(Dx,2)*(3*B10 + pow(Cx,2))) + 4*B00*Cx*Dx*xij + xij*(f36 + Px*pow(Dx,2) + f3)));
	    I[14] += C[1]*W[a]*(Cy*(f36 + 2*B00*Dx*(xij + 2*Cx) + B01*Cx*Ix + pow(Dx,2)*(Cx*Ix + B10) + f3));
	    I[15] += C[1]*W[a]*(Cy*Dz*(Qx*xij + f4));
	    I[16] += C[1]*W[a]*(Qy*(Qx*xij + f4));
	    I[17] += C[1]*W[a]*(Qz*(Qx*xij + f4));
	    I[18] += C[1]*W[a]*(Cz*Dy*(Qx*xij + f4));
	    I[19] += C[1]*W[a]*(Cx*Dy*(f1 + Qz*zij));
	    I[20] += C[1]*W[a]*(Cy*Dx*(f1 + Qz*zij));
	    I[21] += C[1]*W[a]*(Qx*(f1 + Qz*zij));
	    I[22] += C[1]*W[a]*(Qy*(f1 + Qz*zij));
	    I[23] += C[1]*W[a]*(Cz*Qx*(Dy*yij + Qy));
	    I[24] += C[1]*W[a]*(Dx*Qy*(Cz*zij + Pz));
	    I[25] += C[1]*W[a]*(Dz*Py*(Dx*xij + Qx));
	    I[26] += C[1]*W[a]*(Cy*Qz*(Dx*xij + Qx));
	    I[27] += C[1]*W[a]*(Cy*Rz*(Px + Cx*xij));
	    I[28] += C[1]*W[a]*(Dz*Qy*(Px + Cx*xij));
	    I[29] += C[1]*W[a]*(Dz*Px*(Dy*yij + Qy));
	    I[30] += C[1]*W[a]*(Cx*Qz*(Dy*yij + Qy));
	    I[31] += C[1]*W[a]*(f4*(Dy*yij + Qy));
	    I[32] += C[1]*W[a]*(f1*(Dy*yij + Qy));
	    I[33] += C[1]*W[a]*(Dx*Pz*(Dy*yij + Qy));
	    I[34] += C[1]*W[a]*(Dy*Pz*(Dx*xij + Qx));
	    I[35] += C[1]*W[a]*(f1*(Dx*xij + Qx));
	    I[36] += C[1]*W[a]*(f32*(Dx*xij + Qx));
	    I[37] += C[1]*W[a]*(Cz*Qy*(Dx*xij + Qx));
	    I[38] += C[0]*W[a]*(Cz*Dx*Qy);
	    I[39] += C[1]*W[a]*(Cz*Dx*f24);
	    I[40] += C[0]*W[a]*(Dx*Dy*Pz);
	    I[41] += C[1]*W[a]*(Dx*Dy*f40);
	    I[42] += C[1]*W[a]*(Dx*Dz*f31);
	    I[43] += C[0]*W[a]*(Dx*Dz*Py);
	    I[44] += C[1]*W[a]*(Dx*Py*f8);
	    I[45] += C[1]*W[a]*(Dx*Qz*f21);
	    I[46] += C[0]*W[a]*(Cy*Dx*Qz);
	    I[47] += C[1]*W[a]*(Cy*Ix*f33);
	    I[48] += C[1]*W[a]*(Cx*Cy*f30);
	    I[49] += C[0]*W[a]*(Cx*Cy*Rz);
	    I[50] += C[1]*W[a]*(Cx*Rz*f21);
	    I[51] += C[1]*W[a]*(Dz*Qx*f21);
	    I[52] += C[0]*W[a]*(Cy*Dz*Qx);
	    I[53] += C[1]*W[a]*(Cy*Qx*f8);
	    I[54] += C[0]*W[a]*(Cz*Dy*Qx);
	    I[55] += C[1]*W[a]*(Dy*Qx*(Cz*zij + Pz));
	    I[56] += C[1]*W[a]*(Cx*Ry*(Cz*zij + Pz));
	    I[57] += C[1]*W[a]*(Cz*Ry*(Px + Cx*xij));
	    I[58] += C[1]*W[a]*(Dy*Qz*(Px + Cx*xij));
	    I[59] += C[0]*W[a]*(Cx*Dy*Qz);
	    I[60] += C[1]*W[a]*(Dy*Dz*f22);
	    I[61] += C[0]*W[a]*(Dy*Dz*Px);
	    I[62] += C[1]*W[a]*(Dy*Px*f8);
	    I[63] += C[1]*W[a]*(f33*(Px + Cx*xij));
	    I[64] += C[1]*W[a]*(f42*(Px + Cx*xij));
	    I[65] += C[1]*W[a]*(Cx*Iz*f42);
	    I[66] += C[1]*W[a]*(Cx*Dz*f24);
	    I[67] += C[0]*W[a]*(Cx*Dz*Qy);
	    I[68] += C[1]*W[a]*(Cx*Qy*f8);
	    I[69] += C[1]*W[a]*(Cx*Cz*f28);
	    I[70] += C[0]*W[a]*(Cx*Cz*Ry);
	    I[71] += C[1]*W[a]*(Cz*Rx*f21);
	    I[72] += C[1]*W[a]*(Cy*Rx*(Cz*zij + Pz));
	    I[73] += C[1]*W[a]*(f42*(Cz*zij + Pz));
	    I[74] += C[1]*W[a]*(f23*(Cz*zij + Pz));
	    I[75] += C[1]*W[a]*(Cz*Iy*f23);
	    I[76] += C[1]*W[a]*(Cy*Cz*f14);
	    I[77] += C[0]*W[a]*(Cy*Cz*Rx);
	    I[78] += C[1]*W[a]*(Cy*Iz*f23);
	    I[79] += C[0]*W[a]*(Cy*f23);
	    I[80] += C[1]*W[a]*(f21*f23);
	    I[81] += C[0]*W[a]*(Cz*f23);
	    I[82] += C[1]*W[a]*(Cz*Ix*f42);
	    I[83] += C[0]*W[a]*(Cz*f42);
	    I[84] += C[0]*W[a]*(Cx*f42);
	    I[85] += C[1]*W[a]*(Cx*Iy*f33);
	    I[86] += C[0]*W[a]*(Cx*f33);
	    I[87] += C[0]*W[a]*(Cy*f33);
	    I[88] += C[1]*W[a]*(f21*f33);
	    I[89] += C[1]*W[a]*(Rz*f31);
	    I[90] += C[1]*W[a]*(Rx*f31);
	    I[91] += C[1]*W[a]*(Rx*f40);
	    I[92] += C[1]*W[a]*(Ry*f40);
	    I[93] += C[1]*W[a]*(Ry*f22);
	    I[94] += C[1]*W[a]*(Rz*f22);
	    I[95] += C[1]*W[a]*(Ix*Py*Rz);
	    I[96] += C[0]*W[a]*(Py*Rz);
	    I[97] += C[1]*W[a]*(Iy*Px*Rz);
	    I[98] += C[0]*W[a]*(Px*Rz);
	    I[99] += C[1]*W[a]*(Px*f28);
	    I[100] += C[1]*W[a]*(Pz*f28);
	    I[101] += C[1]*W[a]*(Pz*f14);
	    I[102] += C[1]*W[a]*(Py*f14);
	    I[103] += C[1]*W[a]*(Iz*Py*Rx);
	    I[104] += C[0]*W[a]*(Py*Rx);
	    I[105] += C[1]*W[a]*(Iy*Pz*Rx);
	    I[106] += C[0]*W[a]*(Pz*Rx);
	    I[107] += C[1]*W[a]*(Ix*Pz*Ry);
	    I[108] += C[0]*W[a]*(Pz*Ry);
	    I[109] += C[1]*W[a]*(Iz*Px*Ry);
	    I[110] += C[0]*W[a]*(Px*Ry);
	    I[111] += C[1]*W[a]*(Px*f30);
	    I[112] += C[1]*W[a]*(Py*f30);
	    I[113] += C[1]*W[a]*(Dy*Ix*f1);
	    I[114] += C[0]*W[a]*(Dy*f1);
	    I[115] += C[1]*W[a]*(Dz*Ix*f32);
	    I[116] += C[0]*W[a]*(Dz*f32);
	    I[117] += C[1]*W[a]*(Dz*f16);
	    I[118] += C[1]*W[a]*(Ix*Qy*Qz);
	    I[119] += C[0]*W[a]*(Qy*Qz);
	    I[120] += C[1]*W[a]*(Iy*Qx*Qz);
	    I[121] += C[0]*W[a]*(Qx*Qz);
	    I[122] += C[1]*W[a]*(Iz*Qx*Qy);
	    I[123] += C[0]*W[a]*(Qx*Qy);
	    I[124] += C[1]*W[a]*(Qx*f24);
	    I[125] += C[1]*W[a]*(Qz*f24);
	    I[126] += C[1]*W[a]*(Dx*f7);
	    I[127] += C[1]*W[a]*(Dx*Iy*f1);
	    I[128] += C[0]*W[a]*(Dx*f1);
	    I[129] += C[1]*W[a]*(Dx*f16);
	    I[130] += C[1]*W[a]*(Dx*Iz*f32);
	    I[131] += C[0]*W[a]*(Dx*f32);
	    I[132] += C[1]*W[a]*(f32*f8);
	    I[133] += C[1]*W[a]*(f4*f8);
	    I[134] += C[1]*W[a]*(Dz*Iy*f4);
	    I[135] += C[0]*W[a]*(Dz*f4);
	    I[136] += C[1]*W[a]*(Dz*f29);
	    I[137] += C[1]*W[a]*(Dy*f29);
	    I[138] += C[1]*W[a]*(Dy*Iz*f4);
	    I[139] += C[0]*W[a]*(Dy*f4);
	    I[140] += C[1]*W[a]*(Dy*f7);
	    I[141] += C[1]*W[a]*(Iz*f2);
	    I[142] += C[1]*W[a]*(Ix*f2);
	    I[143] += C[0]*W[a]*(f2);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[144]) {
	double T[144];
	for (int i = 0; i < 144; ++i) {
	    T[i] = I[i];
	}
	I[10] = T[0];
	I[39] = T[1];
	I[41] = T[2];
	I[37] = T[3];
	I[71] = T[4];
	I[70] = T[5];
	I[68] = T[6];
	I[62] = T[7];
	I[56] = T[8];
	I[50] = T[9];
	I[18] = T[10];
	I[12] = T[11];
	I[0] = T[12];
	I[6] = T[13];
	I[9] = T[14];
	I[105] = T[15];
	I[81] = T[16];
	I[106] = T[17];
	I[82] = T[18];
	I[142] = T[19];
	I[119] = T[20];
	I[118] = T[21];
	I[143] = T[22];
	I[88] = T[23];
	I[95] = T[24];
	I[103] = T[25];
	I[107] = T[26];
	I[57] = T[27];
	I[129] = T[28];
	I[132] = T[29];
	I[136] = T[30];
	I[84] = T[31];
	I[134] = T[32];
	I[86] = T[33];
	I[80] = T[34];
	I[104] = T[35];
	I[79] = T[36];
	I[83] = T[37];
	I[77] = T[38];
	I[89] = T[39];
	I[74] = T[40];
	I[92] = T[41];
	I[109] = T[42];
	I[97] = T[43];
	I[115] = T[44];
	I[113] = T[45];
	I[101] = T[46];
	I[59] = T[47];
	I[69] = T[48];
	I[51] = T[49];
	I[63] = T[50];
	I[111] = T[51];
	I[99] = T[52];
	I[117] = T[53];
	I[76] = T[54];
	I[94] = T[55];
	I[46] = T[56];
	I[34] = T[57];
	I[130] = T[58];
	I[124] = T[59];
	I[126] = T[60];
	I[120] = T[61];
	I[138] = T[62];
	I[58] = T[63];
	I[33] = T[64];
	I[45] = T[65];
	I[135] = T[66];
	I[123] = T[67];
	I[141] = T[68];
	I[40] = T[69];
	I[28] = T[70];
	I[17] = T[71];
	I[23] = T[72];
	I[47] = T[73];
	I[22] = T[74];
	I[16] = T[75];
	I[11] = T[76];
	I[5] = T[77];
	I[21] = T[78];
	I[3] = T[79];
	I[15] = T[80];
	I[4] = T[81];
	I[35] = T[82];
	I[29] = T[83];
	I[27] = T[84];
	I[64] = T[85];
	I[52] = T[86];
	I[53] = T[87];
	I[65] = T[88];
	I[61] = T[89];
	I[13] = T[90];
	I[20] = T[91];
	I[44] = T[92];
	I[30] = T[93];
	I[54] = T[94];
	I[55] = T[95];
	I[49] = T[96];
	I[60] = T[97];
	I[48] = T[98];
	I[36] = T[99];
	I[38] = T[100];
	I[8] = T[101];
	I[7] = T[102];
	I[19] = T[103];
	I[1] = T[104];
	I[14] = T[105];
	I[2] = T[106];
	I[32] = T[107];
	I[26] = T[108];
	I[42] = T[109];
	I[24] = T[110];
	I[66] = T[111];
	I[67] = T[112];
	I[128] = T[113];
	I[122] = T[114];
	I[127] = T[115];
	I[121] = T[116];
	I[133] = T[117];
	I[131] = T[118];
	I[125] = T[119];
	I[112] = T[120];
	I[100] = T[121];
	I[93] = T[122];
	I[75] = T[123];
	I[87] = T[124];
	I[137] = T[125];
	I[116] = T[126];
	I[110] = T[127];
	I[98] = T[128];
	I[85] = T[129];
	I[91] = T[130];
	I[73] = T[131];
	I[139] = T[132];
	I[114] = T[133];
	I[108] = T[134];
	I[96] = T[135];
	I[102] = T[136];
	I[78] = T[137];
	I[90] = T[138];
	I[72] = T[139];
	I[140] = T[140];
	I[43] = T[141];
	I[31] = T[142];
	I[25] = T[143];
    }

};

template<>
struct impl<meta::braket<rysq::P, rysq::P, rysq::S, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[1], double (&I)[9]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 

// 	for (int a = 0; a < N; ++a) {
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);


	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);

	    I[0] += C[0]*W[a]*((Cx*Ix + B10));
	    I[1] += C[0]*W[a]*((B10 + Cz*Iz));
	    I[2] += C[0]*W[a]*((Cy*Iy + B10));
	    I[3] += C[0]*W[a]*(Cy*Ix);
	    I[4] += C[0]*W[a]*(Cz*Ix);
	    I[5] += C[0]*W[a]*(Cz*Iy);
	    I[6] += C[0]*W[a]*(Cx*Iy);
	    I[7] += C[0]*W[a]*(Cx*Iz);
	    I[8] += C[0]*W[a]*(Cy*Iz);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[9]) {
	double T[9];
	for (int i = 0; i < 9; ++i) {
	    T[i] = I[i];
	}
	I[0] = T[0];
	I[8] = T[1];
	I[4] = T[2];
	I[1] = T[3];
	I[2] = T[4];
	I[5] = T[5];
	I[3] = T[6];
	I[6] = T[7];
	I[7] = T[8];
    }

};

template<>
struct impl<meta::braket<rysq::SP, rysq::S, rysq::S, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<1> &t2, const Vector<1> &W,
			    const double (&C)[2], double (&I)[4]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 1;



// 
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 

// 	for (int a = 0; a < N; ++a) {
// 
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {


	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);



	    I[0] += C[1]*W[a]*(Cx);
	    I[1] += C[1]*W[a]*(Cy);
	    I[2] += C[1]*W[a]*(Cz);
	    I[3] += C[0]*W[a]*(1);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[4]) {
	double T[4];
	for (int i = 0; i < 4; ++i) {
	    T[i] = I[i];
	}
	I[1] = T[0];
	I[2] = T[1];
	I[3] = T[2];
	I[0] = T[3];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::S, rysq::SP, rysq::SP> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[4], double (&I)[160]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;


	double xkl = rkl[0];
	double ykl = rkl[1];
	double zkl = rkl[2];

// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Kx = (xkl + Dx);
	    double Ky = (ykl + Dy);
	    double Kz = (Dz + zkl);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f0 = (3*B10 + pow(Cx,2));
	    double f1 = 3*B00*pow(Cy,2);
	    double f11 = 3*B00*B10;
	    double f12 = (Dy*Py + 2*B00*Cy);
	    double f14 = (2*B00*Cx*(xkl + 2*Dx) + 2*pow(B00,2) + Px*(B01 + Dx*Kx));
	    double f15 = 3*B10*Cy*Dy;
	    double f16 = (B01 + Dx*Kx);
	    double f19 = (Kx*Px + 2*B00*Cx);
	    double f2 = (B01*Cx + Cx*Dx*Kx + B00*(xkl + 2*Dx));
	    double f20 = 3*B00*pow(Cz,2);
	    double f21 = (2*pow(B00,2) + Pz*(Dz*Kz + B01) + 2*B00*Cz*(2*Dz + zkl));
	    double f22 = 3*B00*pow(Cx,2);
	    double f25 = (Dx*Px + 2*B00*Cx);
	    double f26 = 3*B10*Cx*Dx;
	    double f27 = (Dz*Kz + B01);
	    double f28 = (B01 + Dy*Ky);
	    double f29 = (B00 + Cx*Kx);
	    double f3 = (Dz*Pz + 2*B00*Cz);
	    double f36 = Dx*pow(Cx,3);
	    double f37 = (3*B10 + pow(Cz,2));
	    double f38 = (3*B10 + pow(Cy,2));
	    double f4 = 3*B10*Cz*Dz;
	    double f42 = Dy*pow(Cy,3);
	    double f45 = (B01*Cz + B00*(2*Dz + zkl) + Cz*Dz*Kz);
	    double f5 = (Py*(B01 + Dy*Ky) + 2*pow(B00,2) + 2*B00*Cy*(ykl + 2*Dy));
	    double f6 = (Cy*Dy*Ky + B01*Cy + B00*(ykl + 2*Dy));
	    double f7 = Dz*pow(Cz,3);

	    I[0] += C[3]*W[a]*((6*Cy*pow(B00,2) + B01*(3*B10*Cy + pow(Cy,3)) + Dy*(2*f11 + Cy*Dy*f38 + 2*f1) + ykl*(f11 + f15 + f42 + f1)));
	    I[1] += C[3]*W[a]*(Dx*(f11 + f15 + f42 + ykl*pow(Cy,3) + f1 + 3*B10*Cy*ykl));
	    I[2] += C[3]*W[a]*(Dz*(f11 + f15 + f42 + ykl*pow(Cy,3) + f1 + 3*B10*Cy*ykl));
	    I[3] += C[2]*W[a]*((f11 + f15 + f42 + ykl*pow(Cy,3) + f1 + 3*B10*Cy*ykl));
	    I[4] += C[3]*W[a]*((Dz*(2*f20 + Cz*Dz*f37 + 2*f11) + zkl*(f11 + f20 + f4 + f7) + B01*(3*B10*Cz + pow(Cz,3)) + 6*Cz*pow(B00,2)));
	    I[5] += C[3]*W[a]*(Dy*(f11 + f20 + 3*B10*Cz*zkl + f4 + f7 + zkl*pow(Cz,3)));
	    I[6] += C[3]*W[a]*(Dx*(f11 + f20 + 3*B10*Cz*zkl + f4 + f7 + zkl*pow(Cz,3)));
	    I[7] += C[2]*W[a]*((f11 + f20 + 3*B10*Cz*zkl + f4 + f7 + zkl*pow(Cz,3)));
	    I[8] += C[3]*W[a]*((6*Cx*pow(B00,2) + B01*(3*B10*Cx + pow(Cx,3)) + xkl*(f11 + f22 + f26 + f36) + Dx*(Cx*Dx*f0 + 2*f22 + 2*f11)));
	    I[9] += C[3]*W[a]*(Dz*(f11 + f22 + f26 + f36 + 3*B10*Cx*xkl + xkl*pow(Cx,3)));
	    I[10] += C[3]*W[a]*(Dy*(f11 + f22 + f26 + f36 + 3*B10*Cx*xkl + xkl*pow(Cx,3)));
	    I[11] += C[2]*W[a]*((f11 + f22 + f26 + f36 + 3*B10*Cx*xkl + xkl*pow(Cx,3)));
	    I[12] += C[3]*W[a]*(Ky*(f11 + f22 + f26 + f36));
	    I[13] += C[3]*W[a]*(Kz*(f11 + f22 + f26 + f36));
	    I[14] += C[1]*W[a]*((f11 + f22 + f26 + f36));
	    I[15] += C[3]*W[a]*(Kz*(f11 + f15 + f42 + f1));
	    I[16] += C[3]*W[a]*(Kx*(f11 + f15 + f42 + f1));
	    I[17] += C[1]*W[a]*((f11 + f15 + f42 + f1));
	    I[18] += C[3]*W[a]*(Ky*(f11 + f20 + f4 + f7));
	    I[19] += C[3]*W[a]*(Kx*(f11 + f20 + f4 + f7));
	    I[20] += C[1]*W[a]*((f11 + f20 + f4 + f7));
	    I[21] += C[3]*W[a]*(Dx*Py*(Cz*zkl + Qz));
	    I[22] += C[3]*W[a]*(Cx*Qy*(Cz*zkl + Qz));
	    I[23] += C[3]*W[a]*(Qy*(Pz*zkl + f3));
	    I[24] += C[3]*W[a]*(Qx*(Pz*zkl + f3));
	    I[25] += C[3]*W[a]*(Cx*Dy*(Pz*zkl + f3));
	    I[26] += C[2]*W[a]*(Cx*(Pz*zkl + f3));
	    I[27] += C[3]*W[a]*(Cy*Dx*(Pz*zkl + f3));
	    I[28] += C[2]*W[a]*(Cy*(Pz*zkl + f3));
	    I[29] += C[3]*W[a]*(Cy*Qx*(Cz*zkl + Qz));
	    I[30] += C[2]*W[a]*(Cx*Cy*(Cz*zkl + Qz));
	    I[31] += C[3]*W[a]*(Dy*Px*(Cz*zkl + Qz));
	    I[32] += C[2]*W[a]*(Px*(Cz*zkl + Qz));
	    I[33] += C[2]*W[a]*(Py*(Cz*zkl + Qz));
	    I[34] += C[3]*W[a]*(f25*(Cz*zkl + Qz));
	    I[35] += C[3]*W[a]*(f12*(Cz*zkl + Qz));
	    I[36] += C[3]*W[a]*(Cy*f27*f38);
	    I[37] += C[3]*W[a]*(Cy*Px*f27);
	    I[38] += C[3]*W[a]*(Cx*Cz*f6);
	    I[39] += C[3]*W[a]*(Cx*Pz*f28);
	    I[40] += C[3]*W[a]*(Cx*f0*f28);
	    I[41] += C[3]*W[a]*(Cx*f0*f27);
	    I[42] += C[3]*W[a]*(Cx*Py*f27);
	    I[43] += C[3]*W[a]*(Cx*Qz*(Cy*ykl + Qy));
	    I[44] += C[3]*W[a]*(Qz*(f12 + Py*ykl));
	    I[45] += C[3]*W[a]*(Qx*(f12 + Py*ykl));
	    I[46] += C[3]*W[a]*(Cx*Dz*(f12 + Py*ykl));
	    I[47] += C[2]*W[a]*(Cx*(f12 + Py*ykl));
	    I[48] += C[3]*W[a]*(Cz*Dx*(f12 + Py*ykl));
	    I[49] += C[2]*W[a]*(Cz*(f12 + Py*ykl));
	    I[50] += C[3]*W[a]*(Cz*Qx*(Cy*ykl + Qy));
	    I[51] += C[2]*W[a]*(Cx*Cz*(Cy*ykl + Qy));
	    I[52] += C[3]*W[a]*(f3*(Cy*ykl + Qy));
	    I[53] += C[3]*W[a]*(f25*(Cy*ykl + Qy));
	    I[54] += C[3]*W[a]*(Dz*Px*(Cy*ykl + Qy));
	    I[55] += C[2]*W[a]*(Px*(Cy*ykl + Qy));
	    I[56] += C[3]*W[a]*(Dx*Pz*(Cy*ykl + Qy));
	    I[57] += C[2]*W[a]*(Pz*(Cy*ykl + Qy));
	    I[58] += C[3]*W[a]*(Kx*Pz*Qy);
	    I[59] += C[3]*W[a]*(Cz*Qy*f29);
	    I[60] += C[3]*W[a]*(Cy*Cz*f2);
	    I[61] += C[1]*W[a]*(Cy*Cz*Qx);
	    I[62] += C[2]*W[a]*(Cy*Cz*f29);
	    I[63] += C[3]*W[a]*(Cy*Qz*f29);
	    I[64] += C[1]*W[a]*(Cx*Cy*Qz);
	    I[65] += C[3]*W[a]*(Cx*Cy*f45);
	    I[66] += C[0]*W[a]*(Cx*Cy*Cz);
	    I[67] += C[1]*W[a]*(Cx*Cz*Qy);
	    I[68] += C[3]*W[a]*(Qy*f19);
	    I[69] += C[3]*W[a]*(Qz*f19);
	    I[70] += C[3]*W[a]*(Px*f6);
	    I[71] += C[3]*W[a]*(f12*f29);
	    I[72] += C[3]*W[a]*(f29*f3);
	    I[73] += C[3]*W[a]*(Cx*Ky*f3);
	    I[74] += C[1]*W[a]*(Cx*f3);
	    I[75] += C[3]*W[a]*(Cy*Kx*f3);
	    I[76] += C[1]*W[a]*(Cy*f3);
	    I[77] += C[3]*W[a]*(Cy*f14);
	    I[78] += C[3]*W[a]*(Cz*f14);
	    I[79] += C[3]*W[a]*(Cz*Dy*f19);
	    I[80] += C[2]*W[a]*(Cz*f19);
	    I[81] += C[3]*W[a]*(Cz*Kx*f12);
	    I[82] += C[1]*W[a]*(Cz*f12);
	    I[83] += C[3]*W[a]*(Cx*Kz*f12);
	    I[84] += C[1]*W[a]*(Cx*f12);
	    I[85] += C[3]*W[a]*(Cx*f5);
	    I[86] += C[3]*W[a]*(Cz*f5);
	    I[87] += C[3]*W[a]*(Pz*f2);
	    I[88] += C[3]*W[a]*(Py*f2);
	    I[89] += C[3]*W[a]*(Kx*Py*Qz);
	    I[90] += C[1]*W[a]*(Py*Qz);
	    I[91] += C[2]*W[a]*(Py*f29);
	    I[92] += C[3]*W[a]*(Dz*Py*f29);
	    I[93] += C[1]*W[a]*(Cx*Dz*Py);
	    I[94] += C[0]*W[a]*(Cx*Py);
	    I[95] += C[2]*W[a]*(Cx*Kz*Py);
	    I[96] += C[3]*W[a]*(Kz*Py*Qx);
	    I[97] += C[1]*W[a]*(Py*Qx);
	    I[98] += C[3]*W[a]*(Ky*Pz*Qx);
	    I[99] += C[1]*W[a]*(Pz*Qx);
	    I[100] += C[3]*W[a]*(Pz*f6);
	    I[101] += C[1]*W[a]*(Pz*Qy);
	    I[102] += C[3]*W[a]*(Kz*Px*Qy);
	    I[103] += C[1]*W[a]*(Px*Qy);
	    I[104] += C[1]*W[a]*(Px*Qz);
	    I[105] += C[3]*W[a]*(Ky*Px*Qz);
	    I[106] += C[2]*W[a]*(Cz*Ky*Px);
	    I[107] += C[3]*W[a]*(Cz*Px*f28);
	    I[108] += C[3]*W[a]*(Cz*f28*f37);
	    I[109] += C[3]*W[a]*(Cz*f16*f37);
	    I[110] += C[3]*W[a]*(Cz*Py*f16);
	    I[111] += C[2]*W[a]*(Cz*Kx*Py);
	    I[112] += C[2]*W[a]*(Cz*Kx*f37);
	    I[113] += C[3]*W[a]*(Cz*Dy*Kx*f37);
	    I[114] += C[1]*W[a]*(Cz*Dy*f37);
	    I[115] += C[1]*W[a]*(Cz*Dy*Px);
	    I[116] += C[0]*W[a]*(Cz*Px);
	    I[117] += C[3]*W[a]*(Px*f45);
	    I[118] += C[3]*W[a]*(Py*f45);
	    I[119] += C[0]*W[a]*(Cz*Py);
	    I[120] += C[1]*W[a]*(Cz*Dx*Py);
	    I[121] += C[3]*W[a]*(Cz*Dx*Ky*f37);
	    I[122] += C[1]*W[a]*(Cz*Dx*f37);
	    I[123] += C[0]*W[a]*(Cz*f37);
	    I[124] += C[2]*W[a]*(Cz*Ky*f37);
	    I[125] += C[3]*W[a]*(Cz*Ky*f25);
	    I[126] += C[1]*W[a]*(Cz*f25);
	    I[127] += C[3]*W[a]*(Cy*Kz*f25);
	    I[128] += C[1]*W[a]*(Cy*f25);
	    I[129] += C[2]*W[a]*(Cy*f19);
	    I[130] += C[3]*W[a]*(Cy*Dz*f19);
	    I[131] += C[1]*W[a]*(Cy*Dz*Px);
	    I[132] += C[0]*W[a]*(Cy*Px);
	    I[133] += C[2]*W[a]*(Cy*Kz*Px);
	    I[134] += C[3]*W[a]*(Cy*Dx*Kz*f38);
	    I[135] += C[2]*W[a]*(Cy*Kz*f38);
	    I[136] += C[3]*W[a]*(Cy*f16*f38);
	    I[137] += C[3]*W[a]*(Cy*Pz*f16);
	    I[138] += C[2]*W[a]*(Cy*Kx*Pz);
	    I[139] += C[2]*W[a]*(Cy*Kx*f38);
	    I[140] += C[3]*W[a]*(Cy*Dz*Kx*f38);
	    I[141] += C[1]*W[a]*(Cy*Dz*f38);
	    I[142] += C[0]*W[a]*(Cy*f38);
	    I[143] += C[1]*W[a]*(Cy*Dx*f38);
	    I[144] += C[1]*W[a]*(Cy*Dx*Pz);
	    I[145] += C[0]*W[a]*(Cy*Pz);
	    I[146] += C[2]*W[a]*(Pz*f29);
	    I[147] += C[3]*W[a]*(Dy*Pz*f29);
	    I[148] += C[1]*W[a]*(Cx*Dy*Pz);
	    I[149] += C[0]*W[a]*(Cx*Pz);
	    I[150] += C[2]*W[a]*(Cx*Ky*Pz);
	    I[151] += C[2]*W[a]*(Cx*Ky*f0);
	    I[152] += C[3]*W[a]*(Cx*Dz*Ky*f0);
	    I[153] += C[1]*W[a]*(Cx*Dz*f0);
	    I[154] += C[2]*W[a]*(Cx*Kz*f0);
	    I[155] += C[3]*W[a]*(Cx*Dy*Kz*f0);
	    I[156] += C[1]*W[a]*(Cx*Dy*f0);
	    I[157] += C[0]*W[a]*(Cx*f0);
	    I[158] += C[3]*W[a]*(Cx*f21);
	    I[159] += C[3]*W[a]*(Cy*f21);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[160]) {
	double T[160];
	for (int i = 0; i < 160; ++i) {
	    T[i] = I[i];
	}
	I[101] = T[0];
	I[91] = T[1];
	I[111] = T[2];
	I[81] = T[3];
	I[152] = T[4];
	I[142] = T[5];
	I[132] = T[6];
	I[122] = T[7];
	I[50] = T[8];
	I[70] = T[9];
	I[60] = T[10];
	I[40] = T[11];
	I[90] = T[12];
	I[130] = T[13];
	I[10] = T[14];
	I[141] = T[15];
	I[61] = T[16];
	I[21] = T[17];
	I[112] = T[18];
	I[72] = T[19];
	I[32] = T[20];
	I[136] = T[21];
	I[149] = T[22];
	I[148] = T[23];
	I[137] = T[24];
	I[147] = T[25];
	I[127] = T[26];
	I[138] = T[27];
	I[128] = T[28];
	I[139] = T[29];
	I[129] = T[30];
	I[144] = T[31];
	I[124] = T[32];
	I[126] = T[33];
	I[134] = T[34];
	I[146] = T[35];
	I[151] = T[36];
	I[153] = T[37];
	I[109] = T[38];
	I[107] = T[39];
	I[100] = T[40];
	I[150] = T[41];
	I[155] = T[42];
	I[119] = T[43];
	I[116] = T[44];
	I[95] = T[45];
	I[115] = T[46];
	I[85] = T[47];
	I[96] = T[48];
	I[86] = T[49];
	I[99] = T[50];
	I[89] = T[51];
	I[118] = T[52];
	I[93] = T[53];
	I[113] = T[54];
	I[83] = T[55];
	I[98] = T[56];
	I[88] = T[57];
	I[68] = T[58];
	I[69] = T[59];
	I[59] = T[60];
	I[19] = T[61];
	I[49] = T[62];
	I[79] = T[63];
	I[39] = T[64];
	I[159] = T[65];
	I[9] = T[66];
	I[29] = T[67];
	I[63] = T[68];
	I[74] = T[69];
	I[103] = T[70];
	I[65] = T[71];
	I[77] = T[72];
	I[117] = T[73];
	I[37] = T[74];
	I[78] = T[75];
	I[38] = T[76];
	I[53] = T[77];
	I[54] = T[78];
	I[64] = T[79];
	I[44] = T[80];
	I[66] = T[81];
	I[26] = T[82];
	I[145] = T[83];
	I[25] = T[84];
	I[105] = T[85];
	I[106] = T[86];
	I[57] = T[87];
	I[55] = T[88];
	I[76] = T[89];
	I[36] = T[90];
	I[45] = T[91];
	I[75] = T[92];
	I[35] = T[93];
	I[5] = T[94];
	I[125] = T[95];
	I[135] = T[96];
	I[15] = T[97];
	I[97] = T[98];
	I[17] = T[99];
	I[108] = T[100];
	I[28] = T[101];
	I[143] = T[102];
	I[23] = T[103];
	I[34] = T[104];
	I[114] = T[105];
	I[84] = T[106];
	I[104] = T[107];
	I[102] = T[108];
	I[52] = T[109];
	I[56] = T[110];
	I[46] = T[111];
	I[42] = T[112];
	I[62] = T[113];
	I[22] = T[114];
	I[24] = T[115];
	I[4] = T[116];
	I[154] = T[117];
	I[156] = T[118];
	I[6] = T[119];
	I[16] = T[120];
	I[92] = T[121];
	I[12] = T[122];
	I[2] = T[123];
	I[82] = T[124];
	I[94] = T[125];
	I[14] = T[126];
	I[133] = T[127];
	I[13] = T[128];
	I[43] = T[129];
	I[73] = T[130];
	I[33] = T[131];
	I[3] = T[132];
	I[123] = T[133];
	I[131] = T[134];
	I[121] = T[135];
	I[51] = T[136];
	I[58] = T[137];
	I[48] = T[138];
	I[41] = T[139];
	I[71] = T[140];
	I[31] = T[141];
	I[1] = T[142];
	I[11] = T[143];
	I[18] = T[144];
	I[8] = T[145];
	I[47] = T[146];
	I[67] = T[147];
	I[27] = T[148];
	I[7] = T[149];
	I[87] = T[150];
	I[80] = T[151];
	I[110] = T[152];
	I[30] = T[153];
	I[120] = T[154];
	I[140] = T[155];
	I[20] = T[156];
	I[0] = T[157];
	I[157] = T[158];
	I[158] = T[159];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::S, rysq::SP, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[2], double (&I)[40]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;



// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f0 = (3*B10 + pow(Cx,2));
	    double f1 = (Dz*Pz + 2*B00*Cz);
	    double f10 = (3*B10 + pow(Cz,2));
	    double f11 = (3*B10 + pow(Cy,2));
	    double f3 = 3*B00*B10;
	    double f4 = (Dy*Py + 2*B00*Cy);
	    double f6 = (Dx*Px + 2*B00*Cx);

	    I[0] += C[1]*W[a]*((Dx*pow(Cx,3) + 3*B00*pow(Cx,2) + 3*B10*Cx*Dx + f3));
	    I[1] += C[1]*W[a]*((3*B00*pow(Cy,2) + f3 + Dy*pow(Cy,3) + 3*B10*Cy*Dy));
	    I[2] += C[1]*W[a]*((3*B10*Cz*Dz + f3 + 3*B00*pow(Cz,2) + Dz*pow(Cz,3)));
	    I[3] += C[1]*W[a]*(Cy*Cz*Qx);
	    I[4] += C[1]*W[a]*(Cz*Dy*f10);
	    I[5] += C[1]*W[a]*(Cz*Dy*Px);
	    I[6] += C[1]*W[a]*(Cx*Cz*Qy);
	    I[7] += C[1]*W[a]*(Cx*Dz*f0);
	    I[8] += C[1]*W[a]*(Cx*Dz*Py);
	    I[9] += C[1]*W[a]*(Cx*Cy*Qz);
	    I[10] += C[0]*W[a]*(Cx*Cy*Cz);
	    I[11] += C[1]*W[a]*(Cz*Dx*Py);
	    I[12] += C[1]*W[a]*(Cz*Dx*f10);
	    I[13] += C[0]*W[a]*(Cz*f10);
	    I[14] += C[1]*W[a]*(Cy*f6);
	    I[15] += C[1]*W[a]*(Cz*f6);
	    I[16] += C[0]*W[a]*(Cz*Px);
	    I[17] += C[1]*W[a]*(Px*Qy);
	    I[18] += C[1]*W[a]*(Cx*f4);
	    I[19] += C[1]*W[a]*(Cz*f4);
	    I[20] += C[0]*W[a]*(Cz*Py);
	    I[21] += C[1]*W[a]*(Py*Qx);
	    I[22] += C[1]*W[a]*(Pz*Qx);
	    I[23] += C[1]*W[a]*(Pz*Qy);
	    I[24] += C[0]*W[a]*(Cy*Pz);
	    I[25] += C[1]*W[a]*(Cy*Dx*Pz);
	    I[26] += C[1]*W[a]*(Cy*Dx*f11);
	    I[27] += C[0]*W[a]*(Cy*f11);
	    I[28] += C[1]*W[a]*(Cy*Dz*f11);
	    I[29] += C[1]*W[a]*(Cy*Dz*Px);
	    I[30] += C[0]*W[a]*(Cy*Px);
	    I[31] += C[1]*W[a]*(Px*Qz);
	    I[32] += C[1]*W[a]*(Py*Qz);
	    I[33] += C[0]*W[a]*(Cx*Py);
	    I[34] += C[0]*W[a]*(Cx*Pz);
	    I[35] += C[1]*W[a]*(Cx*Dy*Pz);
	    I[36] += C[1]*W[a]*(Cx*Dy*f0);
	    I[37] += C[0]*W[a]*(Cx*f0);
	    I[38] += C[1]*W[a]*(Cx*f1);
	    I[39] += C[1]*W[a]*(Cy*f1);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[40]) {
	double T[40];
	for (int i = 0; i < 40; ++i) {
	    T[i] = I[i];
	}
	I[10] = T[0];
	I[21] = T[1];
	I[32] = T[2];
	I[19] = T[3];
	I[22] = T[4];
	I[24] = T[5];
	I[29] = T[6];
	I[30] = T[7];
	I[35] = T[8];
	I[39] = T[9];
	I[9] = T[10];
	I[16] = T[11];
	I[12] = T[12];
	I[2] = T[13];
	I[13] = T[14];
	I[14] = T[15];
	I[4] = T[16];
	I[23] = T[17];
	I[25] = T[18];
	I[26] = T[19];
	I[6] = T[20];
	I[15] = T[21];
	I[17] = T[22];
	I[28] = T[23];
	I[8] = T[24];
	I[18] = T[25];
	I[11] = T[26];
	I[1] = T[27];
	I[31] = T[28];
	I[33] = T[29];
	I[3] = T[30];
	I[34] = T[31];
	I[36] = T[32];
	I[5] = T[33];
	I[7] = T[34];
	I[27] = T[35];
	I[20] = T[36];
	I[0] = T[37];
	I[37] = T[38];
	I[38] = T[39];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::S, rysq::P, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[1], double (&I)[30]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;



// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f0 = (3*B10 + pow(Cx,2));
	    double f1 = (Dz*Pz + 2*B00*Cz);
	    double f10 = (3*B10 + pow(Cz,2));
	    double f11 = (3*B10 + pow(Cy,2));
	    double f3 = 3*B00*B10;
	    double f4 = (Dy*Py + 2*B00*Cy);
	    double f6 = (Dx*Px + 2*B00*Cx);

	    I[0] += C[0]*W[a]*((Dx*pow(Cx,3) + 3*B00*pow(Cx,2) + 3*B10*Cx*Dx + f3));
	    I[1] += C[0]*W[a]*((3*B00*pow(Cy,2) + f3 + Dy*pow(Cy,3) + 3*B10*Cy*Dy));
	    I[2] += C[0]*W[a]*((3*B10*Cz*Dz + f3 + 3*B00*pow(Cz,2) + Dz*pow(Cz,3)));
	    I[3] += C[0]*W[a]*(Cx*Cz*Qy);
	    I[4] += C[0]*W[a]*(Cy*Dz*Px);
	    I[5] += C[0]*W[a]*(Cy*Dz*f11);
	    I[6] += C[0]*W[a]*(Cy*Dx*f11);
	    I[7] += C[0]*W[a]*(Cy*Dx*Pz);
	    I[8] += C[0]*W[a]*(Cx*Dy*Pz);
	    I[9] += C[0]*W[a]*(Cx*Dy*f0);
	    I[10] += C[0]*W[a]*(Cx*Dz*f0);
	    I[11] += C[0]*W[a]*(Cx*Dz*Py);
	    I[12] += C[0]*W[a]*(Cx*Cy*Qz);
	    I[13] += C[0]*W[a]*(Cy*Cz*Qx);
	    I[14] += C[0]*W[a]*(Cz*Dy*Px);
	    I[15] += C[0]*W[a]*(Cz*Dy*f10);
	    I[16] += C[0]*W[a]*(Cz*Dx*f10);
	    I[17] += C[0]*W[a]*(Cz*Dx*Py);
	    I[18] += C[0]*W[a]*(Py*Qx);
	    I[19] += C[0]*W[a]*(Pz*Qx);
	    I[20] += C[0]*W[a]*(Pz*Qy);
	    I[21] += C[0]*W[a]*(Px*Qy);
	    I[22] += C[0]*W[a]*(Px*Qz);
	    I[23] += C[0]*W[a]*(Py*Qz);
	    I[24] += C[0]*W[a]*(Cy*f6);
	    I[25] += C[0]*W[a]*(Cz*f6);
	    I[26] += C[0]*W[a]*(Cz*f4);
	    I[27] += C[0]*W[a]*(Cx*f4);
	    I[28] += C[0]*W[a]*(Cx*f1);
	    I[29] += C[0]*W[a]*(Cy*f1);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[30]) {
	double T[30];
	for (int i = 0; i < 30; ++i) {
	    T[i] = I[i];
	}
	I[0] = T[0];
	I[11] = T[1];
	I[22] = T[2];
	I[19] = T[3];
	I[23] = T[4];
	I[21] = T[5];
	I[1] = T[6];
	I[8] = T[7];
	I[17] = T[8];
	I[10] = T[9];
	I[20] = T[10];
	I[25] = T[11];
	I[29] = T[12];
	I[9] = T[13];
	I[14] = T[14];
	I[12] = T[15];
	I[2] = T[16];
	I[6] = T[17];
	I[5] = T[18];
	I[7] = T[19];
	I[18] = T[20];
	I[13] = T[21];
	I[24] = T[22];
	I[26] = T[23];
	I[3] = T[24];
	I[4] = T[25];
	I[16] = T[26];
	I[15] = T[27];
	I[27] = T[28];
	I[28] = T[29];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::SP, rysq::SP, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[4], double (&I)[160]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f0 = (3*B10 + pow(Cx,2));
	    double f12 = (B00 + Dz*Iz);
	    double f2 = (Dz*Pz + 2*B00*Cz);
	    double f20 = (3*B10 + pow(Cz,2));
	    double f21 = (3*pow(B10,2) + Iy*pow(Cy,3) + 3*B10*Cy*(yij + 2*Cy));
	    double f23 = (Cy*Iy + B10);
	    double f24 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f25 = (B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10));
	    double f28 = (2*B00*Cy*yij + 3*B00*Py + Dy*(B10*(3*Cy + yij) + Iy*pow(Cy,2)));
	    double f29 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));
	    double f3 = (3*B00*Pz + Cz*Dz*(3*B10 + pow(Cz,2)));
	    double f30 = (Cy*Dy*(3*B10 + pow(Cy,2)) + 3*B00*Py);
	    double f31 = 3*B00*B10;
	    double f32 = (Dy*Py + 2*B00*Cy);
	    double f33 = 3*pow(B10,2);
	    double f36 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f4 = (3*pow(B10,2) + 3*B10*Cx*(xij + 2*Cx) + Ix*pow(Cx,3));
	    double f5 = (Dx*Px + 2*B00*Cx);
	    double f7 = (3*pow(B10,2) + 3*B10*Cz*(2*Cz + zij) + Iz*pow(Cz,3));
	    double f8 = (Cx*Dx*(3*B10 + pow(Cx,2)) + 3*B00*Px);
	    double f9 = (3*B10 + pow(Cy,2));

	    I[0] += C[3]*W[a]*((B00*pow(Cy,2)*(4*Cy + 3*yij) + f31*yij + Dy*f33 + Dy*pow(Cy,4) + Dy*yij*pow(Cy,3) + 4*Cy*f31 + 3*B10*Cy*Dy*(yij + 2*Cy)));
	    I[1] += C[3]*W[a]*(Cz*(xij*(Dx*Px + 2*B00*Cx) + f8));
	    I[2] += C[3]*W[a]*(Cy*(xij*(Dx*Px + 2*B00*Cx) + f8));
	    I[3] += C[3]*W[a]*(Cy*Cz*(Qx*xij + f5));
	    I[4] += C[3]*W[a]*(Py*(Qx*xij + f5));
	    I[5] += C[3]*W[a]*(Pz*(Qx*xij + f5));
	    I[6] += C[3]*W[a]*((Dx*f33 + Dx*xij*pow(Cx,3) + f31*xij + Dx*pow(Cx,4) + 4*Cx*f31 + B00*pow(Cx,2)*(4*Cx + 3*xij) + 3*B10*Cx*Dx*(xij + 2*Cx)));
	    I[7] += C[3]*W[a]*((Dz*f33 + B00*pow(Cz,2)*(4*Cz + 3*zij) + Dz*pow(Cz,4) + 3*B10*Cz*Dz*(2*Cz + zij) + 4*Cz*f31 + Dz*zij*pow(Cz,3) + f31*zij));
	    I[8] += C[3]*W[a]*(Cy*(f3 + zij*(Dz*Pz + 2*B00*Cz)));
	    I[9] += C[3]*W[a]*(Cx*(f3 + zij*(Dz*Pz + 2*B00*Cz)));
	    I[10] += C[3]*W[a]*(Cx*Cy*(f2 + Qz*zij));
	    I[11] += C[3]*W[a]*(Px*(f2 + Qz*zij));
	    I[12] += C[3]*W[a]*(Py*(f2 + Qz*zij));
	    I[13] += C[3]*W[a]*(Cz*f20*(Dy*yij + Qy));
	    I[14] += C[3]*W[a]*(Cz*Px*(Dy*yij + Qy));
	    I[15] += C[3]*W[a]*(Cx*f0*(Dy*yij + Qy));
	    I[16] += C[3]*W[a]*(Cx*Pz*(Dy*yij + Qy));
	    I[17] += C[3]*W[a]*(Cy*Qz*(Px + Cx*xij));
	    I[18] += C[1]*W[a]*(Cy*Cz*(Px + Cx*xij));
	    I[19] += C[3]*W[a]*(Cz*Qy*(Px + Cx*xij));
	    I[20] += C[3]*W[a]*(Cx*Qy*(Cz*zij + Pz));
	    I[21] += C[3]*W[a]*(Dx*Py*(Cz*zij + Pz));
	    I[22] += C[3]*W[a]*(Cz*Py*(Dx*xij + Qx));
	    I[23] += C[3]*W[a]*(Cz*f20*(Dx*xij + Qx));
	    I[24] += C[3]*W[a]*(Cy*f9*(Dx*xij + Qx));
	    I[25] += C[3]*W[a]*(Cy*Pz*(Dx*xij + Qx));
	    I[26] += C[3]*W[a]*(Cy*Qx*(Cz*zij + Pz));
	    I[27] += C[1]*W[a]*(Cx*Cy*(Cz*zij + Pz));
	    I[28] += C[3]*W[a]*(Dy*Px*(Cz*zij + Pz));
	    I[29] += C[1]*W[a]*(Px*(Cz*zij + Pz));
	    I[30] += C[1]*W[a]*(Py*(Cz*zij + Pz));
	    I[31] += C[3]*W[a]*(f5*(Cz*zij + Pz));
	    I[32] += C[3]*W[a]*(f32*(Cz*zij + Pz));
	    I[33] += C[3]*W[a]*(Cy*Dz*f24);
	    I[34] += C[3]*W[a]*(Ix*Py*Qz);
	    I[35] += C[3]*W[a]*(Cy*Ix*f2);
	    I[36] += C[2]*W[a]*(Cx*Cy*Qz);
	    I[37] += C[3]*W[a]*(Cx*Cz*f25);
	    I[38] += C[2]*W[a]*(Cx*Cz*Qy);
	    I[39] += C[3]*W[a]*(Cz*Qx*f23);
	    I[40] += C[2]*W[a]*(Cy*Cz*Qx);
	    I[41] += C[0]*W[a]*(Cx*Cy*Cz);
	    I[42] += C[1]*W[a]*(Cx*Cz*f23);
	    I[43] += C[3]*W[a]*(Cx*Qz*f23);
	    I[44] += C[3]*W[a]*(Cx*f0*f12);
	    I[45] += C[3]*W[a]*(Cy*f12*f9);
	    I[46] += C[3]*W[a]*(Cy*Px*f12);
	    I[47] += C[3]*W[a]*(f2*(Px + Cx*xij));
	    I[48] += C[3]*W[a]*(f32*(Px + Cx*xij));
	    I[49] += C[3]*W[a]*(Dy*Pz*(Px + Cx*xij));
	    I[50] += C[1]*W[a]*(Pz*(Px + Cx*xij));
	    I[51] += C[1]*W[a]*(Py*(Px + Cx*xij));
	    I[52] += C[3]*W[a]*(Dz*Py*(Px + Cx*xij));
	    I[53] += C[2]*W[a]*(Cx*Dz*Py);
	    I[54] += C[3]*W[a]*(Cx*Py*f12);
	    I[55] += C[3]*W[a]*(Cx*f28);
	    I[56] += C[3]*W[a]*(Cz*f28);
	    I[57] += C[3]*W[a]*(Pz*f25);
	    I[58] += C[3]*W[a]*(Px*f25);
	    I[59] += C[3]*W[a]*(Iy*Px*Qz);
	    I[60] += C[2]*W[a]*(Px*Qz);
	    I[61] += C[2]*W[a]*(Py*Qz);
	    I[62] += C[3]*W[a]*(Qx*f36);
	    I[63] += C[3]*W[a]*(Qx*f29);
	    I[64] += C[3]*W[a]*(Iy*Pz*Qx);
	    I[65] += C[2]*W[a]*(Pz*Qx);
	    I[66] += C[3]*W[a]*(Dx*Pz*f23);
	    I[67] += C[1]*W[a]*(Pz*f23);
	    I[68] += C[3]*W[a]*(f23*f5);
	    I[69] += C[3]*W[a]*(Cy*Iz*f5);
	    I[70] += C[2]*W[a]*(Cy*f5);
	    I[71] += C[3]*W[a]*(Cy*Dx*f36);
	    I[72] += C[1]*W[a]*(Cy*f36);
	    I[73] += C[3]*W[a]*(Qy*f36);
	    I[74] += C[3]*W[a]*(Qy*f24);
	    I[75] += C[3]*W[a]*(Iz*Px*Qy);
	    I[76] += C[2]*W[a]*(Px*Qy);
	    I[77] += C[2]*W[a]*(Cz*Dy*Px);
	    I[78] += C[3]*W[a]*(Cz*Dy*Ix*f20);
	    I[79] += C[2]*W[a]*(Cz*Dy*f20);
	    I[80] += C[3]*W[a]*(Cz*Dy*f24);
	    I[81] += C[1]*W[a]*(Cz*f24);
	    I[82] += C[2]*W[a]*(Cz*f5);
	    I[83] += C[3]*W[a]*(Cz*Iy*f5);
	    I[84] += C[1]*W[a]*(Cz*Iy*Px);
	    I[85] += C[0]*W[a]*(Cz*Px);
	    I[86] += C[3]*W[a]*(Cz*Ix*f32);
	    I[87] += C[2]*W[a]*(Cz*f32);
	    I[88] += C[3]*W[a]*(Cx*Iz*f32);
	    I[89] += C[2]*W[a]*(Cx*f32);
	    I[90] += C[3]*W[a]*(Cx*Iy*f2);
	    I[91] += C[2]*W[a]*(Cx*f2);
	    I[92] += C[2]*W[a]*(Cy*f2);
	    I[93] += C[1]*W[a]*(Cy*f24);
	    I[94] += C[3]*W[a]*(Qz*f24);
	    I[95] += C[3]*W[a]*(Ix*f3);
	    I[96] += C[3]*W[a]*(Ix*f30);
	    I[97] += C[3]*W[a]*(Iz*f30);
	    I[98] += C[3]*W[a]*(Iz*f8);
	    I[99] += C[3]*W[a]*(Iy*f8);
	    I[100] += C[3]*W[a]*(Iy*f3);
	    I[101] += C[3]*W[a]*(Cx*Dz*Iy*f0);
	    I[102] += C[2]*W[a]*(Cx*Dz*f0);
	    I[103] += C[3]*W[a]*(Cx*Dz*f29);
	    I[104] += C[1]*W[a]*(Cx*f29);
	    I[105] += C[3]*W[a]*(Qz*f29);
	    I[106] += C[1]*W[a]*(Cz*f29);
	    I[107] += C[3]*W[a]*(Cz*Dx*f29);
	    I[108] += C[2]*W[a]*(Cz*Dx*Py);
	    I[109] += C[2]*W[a]*(Cz*Dx*f20);
	    I[110] += C[3]*W[a]*(Cz*Dx*Iy*f20);
	    I[111] += C[1]*W[a]*(Cz*Iy*f20);
	    I[112] += C[0]*W[a]*(Cz*f20);
	    I[113] += C[1]*W[a]*(Cz*Ix*f20);
	    I[114] += C[1]*W[a]*(Cz*Ix*Py);
	    I[115] += C[0]*W[a]*(Cz*Py);
	    I[116] += C[2]*W[a]*(Py*Qx);
	    I[117] += C[3]*W[a]*(Iz*Py*Qx);
	    I[118] += C[1]*W[a]*(Cx*Iz*Py);
	    I[119] += C[0]*W[a]*(Cx*Py);
	    I[120] += C[1]*W[a]*(Cx*f36);
	    I[121] += C[3]*W[a]*(Cx*Dy*f36);
	    I[122] += C[2]*W[a]*(Cx*Dy*Pz);
	    I[123] += C[2]*W[a]*(Cx*Dy*f0);
	    I[124] += C[3]*W[a]*(Cx*Dy*Iz*f0);
	    I[125] += C[1]*W[a]*(Cx*Iz*f0);
	    I[126] += C[0]*W[a]*(Cx*f0);
	    I[127] += C[1]*W[a]*(Cx*Iy*f0);
	    I[128] += C[1]*W[a]*(Cx*Iy*Pz);
	    I[129] += C[0]*W[a]*(Cx*Pz);
	    I[130] += C[2]*W[a]*(Pz*Qy);
	    I[131] += C[3]*W[a]*(Ix*Pz*Qy);
	    I[132] += C[1]*W[a]*(Cy*Ix*Pz);
	    I[133] += C[0]*W[a]*(Cy*Pz);
	    I[134] += C[2]*W[a]*(Cy*Dx*Pz);
	    I[135] += C[3]*W[a]*(Cy*Dx*Iz*f9);
	    I[136] += C[2]*W[a]*(Cy*Dx*f9);
	    I[137] += C[2]*W[a]*(Cy*Dz*f9);
	    I[138] += C[3]*W[a]*(Cy*Dz*Ix*f9);
	    I[139] += C[1]*W[a]*(Cy*Ix*f9);
	    I[140] += C[0]*W[a]*(Cy*f9);
	    I[141] += C[1]*W[a]*(Cy*Iz*f9);
	    I[142] += C[1]*W[a]*(Cy*Iz*Px);
	    I[143] += C[0]*W[a]*(Cy*Px);
	    I[144] += C[2]*W[a]*(Cy*Dz*Px);
	    I[145] += C[3]*W[a]*(Dz*Px*f23);
	    I[146] += C[1]*W[a]*(Px*f23);
	    I[147] += C[3]*W[a]*(f2*f23);
	    I[148] += C[3]*W[a]*(Dz*f4);
	    I[149] += C[3]*W[a]*(Dz*f21);
	    I[150] += C[3]*W[a]*(Dx*f21);
	    I[151] += C[3]*W[a]*(Dx*f7);
	    I[152] += C[3]*W[a]*(Dy*f7);
	    I[153] += C[3]*W[a]*(Dy*f4);
	    I[154] += C[1]*W[a]*(f4);
	    I[155] += C[1]*W[a]*(f21);
	    I[156] += C[1]*W[a]*(f7);
	    I[157] += C[2]*W[a]*(f8);
	    I[158] += C[2]*W[a]*(f30);
	    I[159] += C[2]*W[a]*(f3);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[160]) {
	double T[160];
	for (int i = 0; i < 160; ++i) {
	    T[i] = I[i];
	}
	I[101] = T[0];
	I[54] = T[1];
	I[53] = T[2];
	I[59] = T[3];
	I[55] = T[4];
	I[57] = T[5];
	I[50] = T[6];
	I[152] = T[7];
	I[158] = T[8];
	I[157] = T[9];
	I[159] = T[10];
	I[154] = T[11];
	I[156] = T[12];
	I[102] = T[13];
	I[104] = T[14];
	I[100] = T[15];
	I[107] = T[16];
	I[139] = T[17];
	I[19] = T[18];
	I[99] = T[19];
	I[119] = T[20];
	I[76] = T[21];
	I[56] = T[22];
	I[52] = T[23];
	I[51] = T[24];
	I[58] = T[25];
	I[79] = T[26];
	I[39] = T[27];
	I[114] = T[28];
	I[34] = T[29];
	I[36] = T[30];
	I[74] = T[31];
	I[116] = T[32];
	I[133] = T[33];
	I[136] = T[34];
	I[138] = T[35];
	I[129] = T[36];
	I[109] = T[37];
	I[89] = T[38];
	I[69] = T[39];
	I[49] = T[40];
	I[9] = T[41];
	I[29] = T[42];
	I[149] = T[43];
	I[150] = T[44];
	I[151] = T[45];
	I[153] = T[46];
	I[137] = T[47];
	I[95] = T[48];
	I[97] = T[49];
	I[17] = T[50];
	I[15] = T[51];
	I[135] = T[52];
	I[125] = T[53];
	I[155] = T[54];
	I[105] = T[55];
	I[106] = T[56];
	I[108] = T[57];
	I[103] = T[58];
	I[144] = T[59];
	I[124] = T[60];
	I[126] = T[61];
	I[77] = T[62];
	I[65] = T[63];
	I[67] = T[64];
	I[47] = T[65];
	I[68] = T[66];
	I[28] = T[67];
	I[63] = T[68];
	I[73] = T[69];
	I[43] = T[70];
	I[78] = T[71];
	I[38] = T[72];
	I[118] = T[73];
	I[93] = T[74];
	I[113] = T[75];
	I[83] = T[76];
	I[84] = T[77];
	I[92] = T[78];
	I[82] = T[79];
	I[94] = T[80];
	I[14] = T[81];
	I[44] = T[82];
	I[64] = T[83];
	I[24] = T[84];
	I[4] = T[85];
	I[96] = T[86];
	I[86] = T[87];
	I[115] = T[88];
	I[85] = T[89];
	I[147] = T[90];
	I[127] = T[91];
	I[128] = T[92];
	I[13] = T[93];
	I[134] = T[94];
	I[132] = T[95];
	I[91] = T[96];
	I[111] = T[97];
	I[70] = T[98];
	I[60] = T[99];
	I[142] = T[100];
	I[140] = T[101];
	I[120] = T[102];
	I[145] = T[103];
	I[25] = T[104];
	I[146] = T[105];
	I[26] = T[106];
	I[66] = T[107];
	I[46] = T[108];
	I[42] = T[109];
	I[62] = T[110];
	I[22] = T[111];
	I[2] = T[112];
	I[12] = T[113];
	I[16] = T[114];
	I[6] = T[115];
	I[45] = T[116];
	I[75] = T[117];
	I[35] = T[118];
	I[5] = T[119];
	I[37] = T[120];
	I[117] = T[121];
	I[87] = T[122];
	I[80] = T[123];
	I[110] = T[124];
	I[30] = T[125];
	I[0] = T[126];
	I[20] = T[127];
	I[27] = T[128];
	I[7] = T[129];
	I[88] = T[130];
	I[98] = T[131];
	I[18] = T[132];
	I[8] = T[133];
	I[48] = T[134];
	I[71] = T[135];
	I[41] = T[136];
	I[121] = T[137];
	I[131] = T[138];
	I[11] = T[139];
	I[1] = T[140];
	I[31] = T[141];
	I[33] = T[142];
	I[3] = T[143];
	I[123] = T[144];
	I[143] = T[145];
	I[23] = T[146];
	I[148] = T[147];
	I[130] = T[148];
	I[141] = T[149];
	I[61] = T[150];
	I[72] = T[151];
	I[112] = T[152];
	I[90] = T[153];
	I[10] = T[154];
	I[21] = T[155];
	I[32] = T[156];
	I[40] = T[157];
	I[81] = T[158];
	I[122] = T[159];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::S, rysq::SP, rysq::P> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[2], double (&I)[72]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;


	double xkl = rkl[0];
	double ykl = rkl[1];
	double zkl = rkl[2];

// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Kx = (xkl + Dx);
	    double Ky = (ykl + Dy);
	    double Kz = (Dz + zkl);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f0 = (Dy*Py + 2*B00*Cy);
	    double f1 = (Dz*Pz + 2*B00*Cz);
	    double f10 = B01*B10;
	    double f13 = (Dx*Px + 2*B00*Cx);
	    double f15 = (Dz*Kz + B01);
	    double f16 = (B00 + Cx*Kx);
	    double f17 = (B00 + Cy*Ky);
	    double f19 = 2*pow(B00,2);
	    double f2 = (Cy*Dy*Ky + B01*Cy + B00*(ykl + 2*Dy));
	    double f21 = (B01 + Dy*Ky);
	    double f4 = (B01*Cx + Cx*Dx*Kx + B00*(xkl + 2*Dx));
	    double f6 = (B01*Cz + B00*(2*Dz + zkl) + Cz*Dz*Kz);
	    double f7 = (B01 + Dx*Kx);
	    double f9 = (Kx*Px + 2*B00*Cx);

	    I[0] += C[1]*W[a]*((B01*pow(Cx,2) + f10 + f19 + 2*B00*Cx*(xkl + 2*Dx) + Dx*Kx*Px));
	    I[1] += C[1]*W[a]*((f10 + f19 + B01*pow(Cz,2) + Dz*Kz*Pz + 2*B00*Cz*(2*Dz + zkl)));
	    I[2] += C[1]*W[a]*((f10 + f19 + Dy*Ky*Py + B01*pow(Cy,2) + 2*B00*Cy*(ykl + 2*Dy)));
	    I[3] += C[1]*W[a]*(Dx*(f0 + Py*ykl));
	    I[4] += C[1]*W[a]*(Dz*(f0 + Py*ykl));
	    I[5] += C[0]*W[a]*((f0 + Py*ykl));
	    I[6] += C[1]*W[a]*(Dy*(Pz*zkl + f1));
	    I[7] += C[1]*W[a]*(Dx*(Pz*zkl + f1));
	    I[8] += C[0]*W[a]*((Pz*zkl + f1));
	    I[9] += C[1]*W[a]*(Cx*Cz*f21);
	    I[10] += C[1]*W[a]*(Cx*Ky*Qz);
	    I[11] += C[1]*W[a]*(Cy*Kx*Qz);
	    I[12] += C[1]*W[a]*(Cy*Kz*Qx);
	    I[13] += C[1]*W[a]*(Dy*Kz*Px);
	    I[14] += C[1]*W[a]*(Cx*Kz*Qy);
	    I[15] += C[0]*W[a]*(Cx*Cy*Kz);
	    I[16] += C[1]*W[a]*(Cx*Cy*f15);
	    I[17] += C[1]*W[a]*(Cy*Cz*f7);
	    I[18] += C[0]*W[a]*(Cy*Cz*Kx);
	    I[19] += C[1]*W[a]*(Cz*Kx*Qy);
	    I[20] += C[1]*W[a]*(Qy*(Cz*zkl + Qz));
	    I[21] += C[1]*W[a]*(Cy*Dx*(Cz*zkl + Qz));
	    I[22] += C[0]*W[a]*(Cy*(Cz*zkl + Qz));
	    I[23] += C[1]*W[a]*(Cx*Dy*(Cz*zkl + Qz));
	    I[24] += C[0]*W[a]*(Cx*(Cz*zkl + Qz));
	    I[25] += C[1]*W[a]*(Qx*(Cz*zkl + Qz));
	    I[26] += C[1]*W[a]*(Cz*Ky*Qx);
	    I[27] += C[0]*W[a]*(Cx*Cz*Ky);
	    I[28] += C[1]*W[a]*(Ky*f13);
	    I[29] += C[1]*W[a]*(Qx*f17);
	    I[30] += C[1]*W[a]*(Kx*f0);
	    I[31] += C[1]*W[a]*(Kx*f1);
	    I[32] += C[1]*W[a]*(Ky*f1);
	    I[33] += C[1]*W[a]*(Qy*f16);
	    I[34] += C[1]*W[a]*(Qz*f16);
	    I[35] += C[1]*W[a]*(Qz*f17);
	    I[36] += C[1]*W[a]*(Dx*Kz*Py);
	    I[37] += C[0]*W[a]*(Kz*Py);
	    I[38] += C[1]*W[a]*(Kz*f13);
	    I[39] += C[1]*W[a]*(Kz*f0);
	    I[40] += C[0]*W[a]*(Kz*Px);
	    I[41] += C[1]*W[a]*(Px*f21);
	    I[42] += C[1]*W[a]*(Pz*f21);
	    I[43] += C[1]*W[a]*(Pz*f7);
	    I[44] += C[1]*W[a]*(Py*f7);
	    I[45] += C[1]*W[a]*(Dz*Kx*Py);
	    I[46] += C[0]*W[a]*(Kx*Py);
	    I[47] += C[1]*W[a]*(Dy*Kx*Pz);
	    I[48] += C[0]*W[a]*(Kx*Pz);
	    I[49] += C[1]*W[a]*(Dx*Ky*Pz);
	    I[50] += C[0]*W[a]*(Ky*Pz);
	    I[51] += C[1]*W[a]*(Dz*Ky*Px);
	    I[52] += C[0]*W[a]*(Ky*Px);
	    I[53] += C[1]*W[a]*(Px*f15);
	    I[54] += C[1]*W[a]*(Py*f15);
	    I[55] += C[1]*W[a]*(Cx*f2);
	    I[56] += C[1]*W[a]*(Cz*f2);
	    I[57] += C[1]*W[a]*(Cz*f4);
	    I[58] += C[1]*W[a]*(Cy*f4);
	    I[59] += C[1]*W[a]*(Cy*Dz*f16);
	    I[60] += C[0]*W[a]*(Cy*f16);
	    I[61] += C[1]*W[a]*(Cz*Dy*f16);
	    I[62] += C[0]*W[a]*(Cz*f16);
	    I[63] += C[1]*W[a]*(Cz*Dx*f17);
	    I[64] += C[0]*W[a]*(Cz*f17);
	    I[65] += C[1]*W[a]*(Cx*Dz*f17);
	    I[66] += C[0]*W[a]*(Cx*f17);
	    I[67] += C[1]*W[a]*(Cx*f6);
	    I[68] += C[1]*W[a]*(Cy*f6);
	    I[69] += C[1]*W[a]*(Dz*f9);
	    I[70] += C[1]*W[a]*(Dy*f9);
	    I[71] += C[0]*W[a]*(f9);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[72]) {
	double T[72];
	for (int i = 0; i < 72; ++i) {
	    T[i] = I[i];
	}
	I[6] = T[0];
	I[68] = T[1];
	I[37] = T[2];
	I[31] = T[3];
	I[43] = T[4];
	I[25] = T[5];
	I[62] = T[6];
	I[56] = T[7];
	I[50] = T[8];
	I[40] = T[9];
	I[46] = T[10];
	I[23] = T[11];
	I[57] = T[12];
	I[60] = T[13];
	I[63] = T[14];
	I[51] = T[15];
	I[69] = T[16];
	I[11] = T[17];
	I[5] = T[18];
	I[17] = T[19];
	I[65] = T[20];
	I[59] = T[21];
	I[53] = T[22];
	I[64] = T[23];
	I[52] = T[24];
	I[58] = T[25];
	I[34] = T[26];
	I[28] = T[27];
	I[30] = T[28];
	I[33] = T[29];
	I[13] = T[30];
	I[20] = T[31];
	I[44] = T[32];
	I[15] = T[33];
	I[22] = T[34];
	I[47] = T[35];
	I[55] = T[36];
	I[49] = T[37];
	I[54] = T[38];
	I[61] = T[39];
	I[48] = T[40];
	I[36] = T[41];
	I[38] = T[42];
	I[8] = T[43];
	I[7] = T[44];
	I[19] = T[45];
	I[1] = T[46];
	I[14] = T[47];
	I[2] = T[48];
	I[32] = T[49];
	I[26] = T[50];
	I[42] = T[51];
	I[24] = T[52];
	I[66] = T[53];
	I[67] = T[54];
	I[39] = T[55];
	I[41] = T[56];
	I[10] = T[57];
	I[9] = T[58];
	I[21] = T[59];
	I[3] = T[60];
	I[16] = T[61];
	I[4] = T[62];
	I[35] = T[63];
	I[29] = T[64];
	I[45] = T[65];
	I[27] = T[66];
	I[70] = T[67];
	I[71] = T[68];
	I[18] = T[69];
	I[12] = T[70];
	I[0] = T[71];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::S, rysq::S, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[1], double (&I)[6]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;



// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 

// 	for (int a = 0; a < N; ++a) {
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);


	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));

	    I[0] += C[0]*W[a]*(Py);
	    I[1] += C[0]*W[a]*(Pz);
	    I[2] += C[0]*W[a]*(Px);
	    I[3] += C[0]*W[a]*(Cx*Cy);
	    I[4] += C[0]*W[a]*(Cx*Cz);
	    I[5] += C[0]*W[a]*(Cy*Cz);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[6]) {
	double T[6];
	for (int i = 0; i < 6; ++i) {
	    T[i] = I[i];
	}
	I[1] = T[0];
	I[2] = T[1];
	I[0] = T[2];
	I[3] = T[3];
	I[4] = T[4];
	I[5] = T[5];
    }

};

template<>
struct impl<meta::braket<rysq::S, rysq::S, rysq::S, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<0> &t2, const Vector<1> &W,
			    const double (&C)[1], double (&I)[1]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 1;



// 
// 
// 

// 

// 

// 	for (int a = 0; a < N; ++a) {
// 
// 
// 

// 

// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {





	    I[0] += C[0]*W[a]*(1);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[1]) {
	double T[1];
	for (int i = 0; i < 1; ++i) {
	    T[i] = I[i];
	}
	I[0] = T[0];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::S, rysq::P, rysq::P> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[1], double (&I)[90]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;


	double xkl = rkl[0];
	double ykl = rkl[1];
	double zkl = rkl[2];

// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Kx = (xkl + Dx);
	    double Ky = (ykl + Dy);
	    double Kz = (Dz + zkl);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f0 = (3*B10 + pow(Cx,2));
	    double f1 = (Cx*Kx*(3*B10 + pow(Cx,2)) + 3*B00*Px);
	    double f10 = (Cz*Kz*(3*B10 + pow(Cz,2)) + 3*B00*Pz);
	    double f11 = (B01*Cz + B00*(2*Dz + zkl) + Cz*Dz*Kz);
	    double f12 = (B01 + Dx*Kx);
	    double f14 = (3*B00*Pz + Cz*Dz*(3*B10 + pow(Cz,2)));
	    double f15 = (Kx*Px + 2*B00*Cx);
	    double f16 = (Dy*Py + 2*B00*Cy);
	    double f17 = (2*pow(B00,2) + Pz*(Dz*Kz + B01) + 2*B00*Cz*(2*Dz + zkl));
	    double f2 = (Dz*Pz + 2*B00*Cz);
	    double f20 = (Dx*Px + 2*B00*Cx);
	    double f22 = (Dz*Kz + B01);
	    double f23 = (B00 + Cx*Kx);
	    double f29 = (Cx*Dx*(3*B10 + pow(Cx,2)) + 3*B00*Px);
	    double f3 = (Py*(B01 + Dy*Ky) + 2*pow(B00,2) + 2*B00*Cy*(ykl + 2*Dy));
	    double f30 = (B01 + Dy*Ky);
	    double f31 = (3*B10 + pow(Cz,2));
	    double f32 = (3*B10 + pow(Cy,2));
	    double f4 = (Cy*Dy*Ky + B01*Cy + B00*(ykl + 2*Dy));
	    double f5 = (3*B00*Py + Cy*Ky*(3*B10 + pow(Cy,2)));
	    double f6 = (Cy*Dy*(3*B10 + pow(Cy,2)) + 3*B00*Py);
	    double f7 = (B01*Cx + Cx*Dx*Kx + B00*(xkl + 2*Dx));
	    double f9 = (2*B00*Cx*(xkl + 2*Dx) + 2*pow(B00,2) + Px*(B01 + Dx*Kx));

	    I[0] += C[0]*W[a]*((3*B00*Px*(xkl + 2*Dx) + Cx*f0*(B01 + Dx*Kx) + 6*Cx*pow(B00,2)));
	    I[1] += C[0]*W[a]*((6*Cy*pow(B00,2) + 3*B00*Py*(ykl + 2*Dy) + Cy*f32*(B01 + Dy*Ky)));
	    I[2] += C[0]*W[a]*((Cz*f31*(Dz*Kz + B01) + 3*B00*Pz*(2*Dz + zkl) + 6*Cz*pow(B00,2)));
	    I[3] += C[0]*W[a]*(Dx*Py*(Cz*zkl + Qz));
	    I[4] += C[0]*W[a]*(Cy*Qx*(Cz*zkl + Qz));
	    I[5] += C[0]*W[a]*(Dy*Px*(Cz*zkl + Qz));
	    I[6] += C[0]*W[a]*(Cx*Dy*(Pz*zkl + f2));
	    I[7] += C[0]*W[a]*(Cy*Dx*(Pz*zkl + f2));
	    I[8] += C[0]*W[a]*(Qx*(Pz*zkl + f2));
	    I[9] += C[0]*W[a]*(Qy*(Pz*zkl + f2));
	    I[10] += C[0]*W[a]*(Cx*Qy*(Cz*zkl + Qz));
	    I[11] += C[0]*W[a]*(Cx*Qz*(Cy*ykl + Qy));
	    I[12] += C[0]*W[a]*(Dz*Px*(Cy*ykl + Qy));
	    I[13] += C[0]*W[a]*(Cx*Dz*(f16 + Py*ykl));
	    I[14] += C[0]*W[a]*(Cz*Dx*(f16 + Py*ykl));
	    I[15] += C[0]*W[a]*(Qz*(f16 + Py*ykl));
	    I[16] += C[0]*W[a]*(Qx*(f16 + Py*ykl));
	    I[17] += C[0]*W[a]*(Cz*Qx*(Cy*ykl + Qy));
	    I[18] += C[0]*W[a]*(Dx*Pz*(Cy*ykl + Qy));
	    I[19] += C[0]*W[a]*(Cy*Dx*Kz*f32);
	    I[20] += C[0]*W[a]*(Cx*Dz*Ky*f0);
	    I[21] += C[0]*W[a]*(Cx*Dy*Kz*f0);
	    I[22] += C[0]*W[a]*(Cx*Cz*f4);
	    I[23] += C[0]*W[a]*(Cy*Qz*f23);
	    I[24] += C[0]*W[a]*(Ky*Px*Qz);
	    I[25] += C[0]*W[a]*(Cx*Ky*f2);
	    I[26] += C[0]*W[a]*(Ky*Pz*Qx);
	    I[27] += C[0]*W[a]*(Kz*Py*Qx);
	    I[28] += C[0]*W[a]*(Kz*Px*Qy);
	    I[29] += C[0]*W[a]*(Cx*Kz*f16);
	    I[30] += C[0]*W[a]*(Cy*Kz*f20);
	    I[31] += C[0]*W[a]*(f20*(Cy*ykl + Qy));
	    I[32] += C[0]*W[a]*(f2*(Cy*ykl + Qy));
	    I[33] += C[0]*W[a]*(Cy*Kx*f2);
	    I[34] += C[0]*W[a]*(Kx*Py*Qz);
	    I[35] += C[0]*W[a]*(Dz*Py*f23);
	    I[36] += C[0]*W[a]*(Cz*Qy*f23);
	    I[37] += C[0]*W[a]*(Cy*Cz*f7);
	    I[38] += C[0]*W[a]*(Cy*Dz*f15);
	    I[39] += C[0]*W[a]*(Cy*Dz*Kx*f32);
	    I[40] += C[0]*W[a]*(Cy*f22*f32);
	    I[41] += C[0]*W[a]*(Cy*Px*f22);
	    I[42] += C[0]*W[a]*(Cz*Px*f30);
	    I[43] += C[0]*W[a]*(Cz*f30*f31);
	    I[44] += C[0]*W[a]*(Cz*Dx*Ky*f31);
	    I[45] += C[0]*W[a]*(Cz*Ky*f20);
	    I[46] += C[0]*W[a]*(f20*(Cz*zkl + Qz));
	    I[47] += C[0]*W[a]*(f16*(Cz*zkl + Qz));
	    I[48] += C[0]*W[a]*(Cz*Kx*f16);
	    I[49] += C[0]*W[a]*(Kx*Pz*Qy);
	    I[50] += C[0]*W[a]*(Dy*Pz*f23);
	    I[51] += C[0]*W[a]*(Cx*Pz*f30);
	    I[52] += C[0]*W[a]*(Cx*f0*f30);
	    I[53] += C[0]*W[a]*(Cx*f0*f22);
	    I[54] += C[0]*W[a]*(Cx*Py*f22);
	    I[55] += C[0]*W[a]*(Cx*Cy*f11);
	    I[56] += C[0]*W[a]*(Cy*Pz*f12);
	    I[57] += C[0]*W[a]*(Cy*f12*f32);
	    I[58] += C[0]*W[a]*(Cz*Py*f12);
	    I[59] += C[0]*W[a]*(Cz*f12*f31);
	    I[60] += C[0]*W[a]*(Cz*Dy*Kx*f31);
	    I[61] += C[0]*W[a]*(Cz*Dy*f15);
	    I[62] += C[0]*W[a]*(Qy*f15);
	    I[63] += C[0]*W[a]*(Qz*f15);
	    I[64] += C[0]*W[a]*(f16*f23);
	    I[65] += C[0]*W[a]*(f2*f23);
	    I[66] += C[0]*W[a]*(Kx*f6);
	    I[67] += C[0]*W[a]*(Kx*f14);
	    I[68] += C[0]*W[a]*(Ky*f14);
	    I[69] += C[0]*W[a]*(Ky*f29);
	    I[70] += C[0]*W[a]*(Kz*f29);
	    I[71] += C[0]*W[a]*(Kz*f6);
	    I[72] += C[0]*W[a]*(Dx*f10);
	    I[73] += C[0]*W[a]*(Dx*f5);
	    I[74] += C[0]*W[a]*(Dz*f5);
	    I[75] += C[0]*W[a]*(Dz*f1);
	    I[76] += C[0]*W[a]*(Dy*f1);
	    I[77] += C[0]*W[a]*(Dy*f10);
	    I[78] += C[0]*W[a]*(Py*f7);
	    I[79] += C[0]*W[a]*(Pz*f7);
	    I[80] += C[0]*W[a]*(Pz*f4);
	    I[81] += C[0]*W[a]*(Px*f4);
	    I[82] += C[0]*W[a]*(Px*f11);
	    I[83] += C[0]*W[a]*(Py*f11);
	    I[84] += C[0]*W[a]*(Cy*f9);
	    I[85] += C[0]*W[a]*(Cz*f9);
	    I[86] += C[0]*W[a]*(Cz*f3);
	    I[87] += C[0]*W[a]*(Cx*f3);
	    I[88] += C[0]*W[a]*(Cx*f17);
	    I[89] += C[0]*W[a]*(Cy*f17);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[90]) {
	double T[90];
	for (int i = 0; i < 90; ++i) {
	    T[i] = I[i];
	}
	I[0] = T[0];
	I[41] = T[1];
	I[82] = T[2];
	I[66] = T[3];
	I[69] = T[4];
	I[74] = T[5];
	I[77] = T[6];
	I[68] = T[7];
	I[67] = T[8];
	I[78] = T[9];
	I[79] = T[10];
	I[59] = T[11];
	I[53] = T[12];
	I[55] = T[13];
	I[36] = T[14];
	I[56] = T[15];
	I[35] = T[16];
	I[39] = T[17];
	I[38] = T[18];
	I[61] = T[19];
	I[50] = T[20];
	I[70] = T[21];
	I[49] = T[22];
	I[29] = T[23];
	I[54] = T[24];
	I[57] = T[25];
	I[37] = T[26];
	I[65] = T[27];
	I[73] = T[28];
	I[75] = T[29];
	I[63] = T[30];
	I[33] = T[31];
	I[58] = T[32];
	I[28] = T[33];
	I[26] = T[34];
	I[25] = T[35];
	I[19] = T[36];
	I[9] = T[37];
	I[23] = T[38];
	I[21] = T[39];
	I[81] = T[40];
	I[83] = T[41];
	I[44] = T[42];
	I[42] = T[43];
	I[32] = T[44];
	I[34] = T[45];
	I[64] = T[46];
	I[76] = T[47];
	I[16] = T[48];
	I[18] = T[49];
	I[17] = T[50];
	I[47] = T[51];
	I[40] = T[52];
	I[80] = T[53];
	I[85] = T[54];
	I[89] = T[55];
	I[8] = T[56];
	I[1] = T[57];
	I[6] = T[58];
	I[2] = T[59];
	I[12] = T[60];
	I[14] = T[61];
	I[13] = T[62];
	I[24] = T[63];
	I[15] = T[64];
	I[27] = T[65];
	I[11] = T[66];
	I[22] = T[67];
	I[52] = T[68];
	I[30] = T[69];
	I[60] = T[70];
	I[71] = T[71];
	I[62] = T[72];
	I[31] = T[73];
	I[51] = T[74];
	I[20] = T[75];
	I[10] = T[76];
	I[72] = T[77];
	I[5] = T[78];
	I[7] = T[79];
	I[48] = T[80];
	I[43] = T[81];
	I[84] = T[82];
	I[86] = T[83];
	I[3] = T[84];
	I[4] = T[85];
	I[46] = T[86];
	I[45] = T[87];
	I[87] = T[88];
	I[88] = T[89];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::S, rysq::P, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[1], double (&I)[18]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;



// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);

	    I[0] += C[0]*W[a]*((Dx*Px + 2*B00*Cx));
	    I[1] += C[0]*W[a]*((Dy*Py + 2*B00*Cy));
	    I[2] += C[0]*W[a]*((Dz*Pz + 2*B00*Cz));
	    I[3] += C[0]*W[a]*(Cx*Cz*Dy);
	    I[4] += C[0]*W[a]*(Cx*Cy*Dz);
	    I[5] += C[0]*W[a]*(Cy*Cz*Dx);
	    I[6] += C[0]*W[a]*(Dx*Py);
	    I[7] += C[0]*W[a]*(Dx*Pz);
	    I[8] += C[0]*W[a]*(Dy*Pz);
	    I[9] += C[0]*W[a]*(Dy*Px);
	    I[10] += C[0]*W[a]*(Dz*Px);
	    I[11] += C[0]*W[a]*(Dz*Py);
	    I[12] += C[0]*W[a]*(Cy*Qx);
	    I[13] += C[0]*W[a]*(Cz*Qx);
	    I[14] += C[0]*W[a]*(Cz*Qy);
	    I[15] += C[0]*W[a]*(Cx*Qy);
	    I[16] += C[0]*W[a]*(Cx*Qz);
	    I[17] += C[0]*W[a]*(Cy*Qz);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[18]) {
	double T[18];
	for (int i = 0; i < 18; ++i) {
	    T[i] = I[i];
	}
	I[0] = T[0];
	I[7] = T[1];
	I[14] = T[2];
	I[10] = T[3];
	I[15] = T[4];
	I[5] = T[5];
	I[1] = T[6];
	I[2] = T[7];
	I[8] = T[8];
	I[6] = T[9];
	I[12] = T[10];
	I[13] = T[11];
	I[3] = T[12];
	I[4] = T[13];
	I[11] = T[14];
	I[9] = T[15];
	I[16] = T[16];
	I[17] = T[17];
    }

};

template<>
struct impl<meta::braket<rysq::P, rysq::S, rysq::S, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<1> &t2, const Vector<1> &W,
			    const double (&C)[1], double (&I)[3]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 1;



// 
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 

// 	for (int a = 0; a < N; ++a) {
// 
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {


	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);



	    I[0] += C[0]*W[a]*(Cx);
	    I[1] += C[0]*W[a]*(Cy);
	    I[2] += C[0]*W[a]*(Cz);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[3]) {
	double T[3];
	for (int i = 0; i < 3; ++i) {
	    T[i] = I[i];
	}
	I[0] = T[0];
	I[1] = T[1];
	I[2] = T[2];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::SP, rysq::S, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[2], double (&I)[24]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 

// 	for (int a = 0; a < N; ++a) {
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);


	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));

	    I[0] += C[1]*W[a]*(Iy*Px);
	    I[1] += C[1]*W[a]*(Iy*Pz);
	    I[2] += C[1]*W[a]*(Cx*Cz*Iy);
	    I[3] += C[1]*W[a]*(Cy*Cz*Ix);
	    I[4] += C[1]*W[a]*(Cy*(B10 + Cz*Iz));
	    I[5] += C[1]*W[a]*(Cz*(Cy*Iy + B10));
	    I[6] += C[1]*W[a]*(Cx*(Cy*Iy + B10));
	    I[7] += C[1]*W[a]*((B10*(3*Cy + yij) + Iy*pow(Cy,2)));
	    I[8] += C[1]*W[a]*(Iz*Py);
	    I[9] += C[1]*W[a]*(Cy*(Cx*Ix + B10));
	    I[10] += C[1]*W[a]*((B10*(3*Cx + xij) + Ix*pow(Cx,2)));
	    I[11] += C[1]*W[a]*(Cz*(Cx*Ix + B10));
	    I[12] += C[1]*W[a]*(Cx*(B10 + Cz*Iz));
	    I[13] += C[1]*W[a]*((Iz*pow(Cz,2) + B10*(3*Cz + zij)));
	    I[14] += C[1]*W[a]*(Ix*Pz);
	    I[15] += C[1]*W[a]*(Ix*Py);
	    I[16] += C[0]*W[a]*(Py);
	    I[17] += C[0]*W[a]*(Pz);
	    I[18] += C[0]*W[a]*(Px);
	    I[19] += C[1]*W[a]*(Iz*Px);
	    I[20] += C[1]*W[a]*(Cx*Cy*Iz);
	    I[21] += C[0]*W[a]*(Cx*Cy);
	    I[22] += C[0]*W[a]*(Cx*Cz);
	    I[23] += C[0]*W[a]*(Cy*Cz);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[24]) {
	double T[24];
	for (int i = 0; i < 24; ++i) {
	    T[i] = I[i];
	}
	I[12] = T[0];
	I[14] = T[1];
	I[16] = T[2];
	I[11] = T[3];
	I[23] = T[4];
	I[17] = T[5];
	I[15] = T[6];
	I[13] = T[7];
	I[19] = T[8];
	I[9] = T[9];
	I[6] = T[10];
	I[10] = T[11];
	I[22] = T[12];
	I[20] = T[13];
	I[8] = T[14];
	I[7] = T[15];
	I[1] = T[16];
	I[2] = T[17];
	I[0] = T[18];
	I[18] = T[19];
	I[21] = T[20];
	I[3] = T[21];
	I[4] = T[22];
	I[5] = T[23];
    }

};

template<>
struct impl<meta::braket<rysq::SP, rysq::S, rysq::SP, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[4], double (&I)[16]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;



// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);

	    I[0] += C[3]*W[a]*(Qx);
	    I[1] += C[3]*W[a]*(Qy);
	    I[2] += C[3]*W[a]*(Qz);
	    I[3] += C[3]*W[a]*(Cx*Dz);
	    I[4] += C[3]*W[a]*(Cy*Dz);
	    I[5] += C[3]*W[a]*(Cy*Dx);
	    I[6] += C[3]*W[a]*(Cz*Dx);
	    I[7] += C[3]*W[a]*(Cz*Dy);
	    I[8] += C[3]*W[a]*(Cx*Dy);
	    I[9] += C[1]*W[a]*(Cx);
	    I[10] += C[1]*W[a]*(Cy);
	    I[11] += C[1]*W[a]*(Cz);
	    I[12] += C[2]*W[a]*(Dx);
	    I[13] += C[2]*W[a]*(Dy);
	    I[14] += C[2]*W[a]*(Dz);
	    I[15] += C[0]*W[a]*(1);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[16]) {
	double T[16];
	for (int i = 0; i < 16; ++i) {
	    T[i] = I[i];
	}
	I[5] = T[0];
	I[10] = T[1];
	I[15] = T[2];
	I[13] = T[3];
	I[14] = T[4];
	I[6] = T[5];
	I[7] = T[6];
	I[11] = T[7];
	I[9] = T[8];
	I[1] = T[9];
	I[2] = T[10];
	I[3] = T[11];
	I[4] = T[12];
	I[8] = T[13];
	I[12] = T[14];
	I[0] = T[15];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::SP, rysq::P, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[2], double (&I)[120]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f0 = (3*B10 + pow(Cx,2));
	    double f10 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));
	    double f12 = (3*pow(B10,2) + 3*B10*Cz*(2*Cz + zij) + Iz*pow(Cz,3));
	    double f15 = (Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx));
	    double f17 = (Dy*Iy + B00);
	    double f18 = (Cy*Iy + B10);
	    double f19 = (3*pow(B10,2) + Iy*pow(Cy,3) + 3*B10*Cy*(yij + 2*Cy));
	    double f2 = (Dz*Pz + 2*B00*Cz);
	    double f20 = (3*pow(B10,2) + 3*B10*Cx*(xij + 2*Cx) + Ix*pow(Cx,3));
	    double f22 = (Dx*Px + 2*B00*Cx);
	    double f23 = (Dx*Ix + B00);
	    double f24 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f3 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f33 = (3*B10 + pow(Cz,2));
	    double f34 = (3*B10 + pow(Cy,2));
	    double f4 = (B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10));
	    double f5 = (Dx*(B10*(3*Cx + xij) + Ix*pow(Cx,2)) + 2*B00*Cx*xij + 3*B00*Px);
	    double f6 = (Cy*Dy*(3*B10 + pow(Cy,2)) + 3*B00*Py);
	    double f8 = 3*B00*B10;
	    double f9 = (Dy*Py + 2*B00*Cy);

	    I[0] += C[1]*W[a]*(Cx*(f6 + yij*(Dy*Py + 2*B00*Cy)));
	    I[1] += C[1]*W[a]*(Cz*(f6 + yij*(Dy*Py + 2*B00*Cy)));
	    I[2] += C[1]*W[a]*((B00*pow(Cy,2)*(4*Cy + 3*yij) + 3*Dy*pow(B10,2) + 4*Cy*f8 + Dy*pow(Cy,4) + Dy*yij*pow(Cy,3) + f8*yij + 3*B10*Cy*Dy*(yij + 2*Cy)));
	    I[3] += C[1]*W[a]*(Cx*Cy*(f2 + Qz*zij));
	    I[4] += C[1]*W[a]*(Px*(f2 + Qz*zij));
	    I[5] += C[1]*W[a]*(Py*(f2 + Qz*zij));
	    I[6] += C[1]*W[a]*(Cy*(B00*Cz*(3*Cz + 2*zij) + Dz*Iz*pow(Cz,2) + B10*Dz*(3*Cz + zij) + f8));
	    I[7] += C[1]*W[a]*(Cx*(B00*Cz*(3*Cz + 2*zij) + Dz*Iz*pow(Cz,2) + B10*Dz*(3*Cz + zij) + f8));
	    I[8] += C[1]*W[a]*((f8*zij + 3*Dz*pow(B10,2) + B00*pow(Cz,2)*(4*Cz + 3*zij) + Dz*pow(Cz,4) + 3*B10*Cz*Dz*(2*Cz + zij) + 4*Cz*f8 + Dz*zij*pow(Cz,3)));
	    I[9] += C[1]*W[a]*(Iy*(3*B10*Cz*Dz + f8 + 3*B00*pow(Cz,2) + Dz*pow(Cz,3)));
	    I[10] += C[1]*W[a]*(Ix*(3*B10*Cz*Dz + f8 + 3*B00*pow(Cz,2) + Dz*pow(Cz,3)));
	    I[11] += C[0]*W[a]*((3*B10*Cz*Dz + f8 + 3*B00*pow(Cz,2) + Dz*pow(Cz,3)));
	    I[12] += C[1]*W[a]*(Iz*(Dx*pow(Cx,3) + 3*B00*pow(Cx,2) + 3*B10*Cx*Dx + f8));
	    I[13] += C[1]*W[a]*(Iy*(Dx*pow(Cx,3) + 3*B00*pow(Cx,2) + 3*B10*Cx*Dx + f8));
	    I[14] += C[0]*W[a]*((Dx*pow(Cx,3) + 3*B00*pow(Cx,2) + 3*B10*Cx*Dx + f8));
	    I[15] += C[1]*W[a]*((Dx*xij*pow(Cx,3) + f8*xij + 3*Dx*pow(B10,2) + Dx*pow(Cx,4) + 4*Cx*f8 + B00*pow(Cx,2)*(4*Cx + 3*xij) + 3*B10*Cx*Dx*(xij + 2*Cx)));
	    I[16] += C[1]*W[a]*(Cy*Qz*(Px + Cx*xij));
	    I[17] += C[1]*W[a]*(Cy*f34*(Dz*zij + Qz));
	    I[18] += C[1]*W[a]*(Cy*Px*(Dz*zij + Qz));
	    I[19] += C[1]*W[a]*(Cx*f0*(Dz*zij + Qz));
	    I[20] += C[1]*W[a]*(Cx*Py*(Dz*zij + Qz));
	    I[21] += C[1]*W[a]*(Dz*Py*(Px + Cx*xij));
	    I[22] += C[1]*W[a]*(Dy*Pz*(Px + Cx*xij));
	    I[23] += C[1]*W[a]*(f9*(Px + Cx*xij));
	    I[24] += C[1]*W[a]*(f2*(Px + Cx*xij));
	    I[25] += C[1]*W[a]*(Cz*Qy*(Px + Cx*xij));
	    I[26] += C[1]*W[a]*(Cx*Qy*(Cz*zij + Pz));
	    I[27] += C[1]*W[a]*(Cy*Qx*(Cz*zij + Pz));
	    I[28] += C[1]*W[a]*(Dx*Py*(Cz*zij + Pz));
	    I[29] += C[1]*W[a]*(f9*(Cz*zij + Pz));
	    I[30] += C[1]*W[a]*(f22*(Cz*zij + Pz));
	    I[31] += C[1]*W[a]*(Dy*Px*(Cz*zij + Pz));
	    I[32] += C[1]*W[a]*(Cz*Dy*f24);
	    I[33] += C[1]*W[a]*(Cx*Dy*f3);
	    I[34] += C[0]*W[a]*(Cx*Dz*Py);
	    I[35] += C[1]*W[a]*(Dz*Px*f18);
	    I[36] += C[1]*W[a]*(Cz*Qx*f18);
	    I[37] += C[1]*W[a]*(Cx*Cz*f4);
	    I[38] += C[0]*W[a]*(Cx*Cz*Qy);
	    I[39] += C[1]*W[a]*(Cx*f0*f17);
	    I[40] += C[1]*W[a]*(Cx*Pz*f17);
	    I[41] += C[0]*W[a]*(Cx*Dy*Pz);
	    I[42] += C[1]*W[a]*(Cx*Dy*Iz*f0);
	    I[43] += C[0]*W[a]*(Cx*Dy*f0);
	    I[44] += C[1]*W[a]*(Cx*Dz*Iy*f0);
	    I[45] += C[0]*W[a]*(Cx*Dz*f0);
	    I[46] += C[1]*W[a]*(Cx*Dz*f10);
	    I[47] += C[0]*W[a]*(Cx*Cy*Qz);
	    I[48] += C[1]*W[a]*(Cx*Qz*f18);
	    I[49] += C[1]*W[a]*(Dx*Pz*f18);
	    I[50] += C[1]*W[a]*(Cz*Dx*f10);
	    I[51] += C[1]*W[a]*(Cy*Dx*f3);
	    I[52] += C[1]*W[a]*(Cy*Dz*f24);
	    I[53] += C[0]*W[a]*(Cy*Dz*Px);
	    I[54] += C[1]*W[a]*(Cy*Dz*Ix*f34);
	    I[55] += C[0]*W[a]*(Cy*Dz*f34);
	    I[56] += C[1]*W[a]*(Cy*f23*f34);
	    I[57] += C[1]*W[a]*(Cy*Pz*f23);
	    I[58] += C[0]*W[a]*(Cy*Dx*Pz);
	    I[59] += C[1]*W[a]*(Cy*Dx*Iz*f34);
	    I[60] += C[0]*W[a]*(Cy*Dx*f34);
	    I[61] += C[1]*W[a]*(Cy*Cz*f15);
	    I[62] += C[0]*W[a]*(Cy*Cz*Qx);
	    I[63] += C[1]*W[a]*(Cz*Px*f17);
	    I[64] += C[1]*W[a]*(Cz*f17*f33);
	    I[65] += C[1]*W[a]*(Cz*f23*f33);
	    I[66] += C[1]*W[a]*(Cz*Py*f23);
	    I[67] += C[0]*W[a]*(Cz*Dx*Py);
	    I[68] += C[1]*W[a]*(Cz*Dx*Iy*f33);
	    I[69] += C[0]*W[a]*(Cz*Dx*f33);
	    I[70] += C[1]*W[a]*(Cz*Dy*Ix*f33);
	    I[71] += C[0]*W[a]*(Cz*Dy*f33);
	    I[72] += C[0]*W[a]*(Cz*Dy*Px);
	    I[73] += C[1]*W[a]*(Px*f4);
	    I[74] += C[1]*W[a]*(Pz*f4);
	    I[75] += C[1]*W[a]*(Qz*f24);
	    I[76] += C[1]*W[a]*(Dz*f20);
	    I[77] += C[1]*W[a]*(Dy*f20);
	    I[78] += C[1]*W[a]*(Dy*f12);
	    I[79] += C[1]*W[a]*(Dx*f12);
	    I[80] += C[1]*W[a]*(Dx*f19);
	    I[81] += C[1]*W[a]*(Dz*f19);
	    I[82] += C[1]*W[a]*(Qx*f3);
	    I[83] += C[1]*W[a]*(Qy*f3);
	    I[84] += C[1]*W[a]*(Qy*f24);
	    I[85] += C[1]*W[a]*(Ix*Pz*Qy);
	    I[86] += C[0]*W[a]*(Pz*Qy);
	    I[87] += C[1]*W[a]*(Iz*Px*Qy);
	    I[88] += C[0]*W[a]*(Px*Qy);
	    I[89] += C[1]*W[a]*(Iy*Px*Qz);
	    I[90] += C[0]*W[a]*(Px*Qz);
	    I[91] += C[1]*W[a]*(Ix*Py*Qz);
	    I[92] += C[0]*W[a]*(Py*Qz);
	    I[93] += C[1]*W[a]*(Py*f15);
	    I[94] += C[1]*W[a]*(Pz*f15);
	    I[95] += C[1]*W[a]*(Iy*Pz*Qx);
	    I[96] += C[0]*W[a]*(Pz*Qx);
	    I[97] += C[1]*W[a]*(Iz*Py*Qx);
	    I[98] += C[0]*W[a]*(Py*Qx);
	    I[99] += C[1]*W[a]*(Qx*f10);
	    I[100] += C[1]*W[a]*(Qz*f10);
	    I[101] += C[1]*W[a]*(Cz*Ix*f9);
	    I[102] += C[0]*W[a]*(Cz*f9);
	    I[103] += C[1]*W[a]*(Cx*Iz*f9);
	    I[104] += C[0]*W[a]*(Cx*f9);
	    I[105] += C[1]*W[a]*(Cx*Iy*f2);
	    I[106] += C[0]*W[a]*(Cx*f2);
	    I[107] += C[1]*W[a]*(Cy*Ix*f2);
	    I[108] += C[0]*W[a]*(Cy*f2);
	    I[109] += C[1]*W[a]*(Cy*f5);
	    I[110] += C[1]*W[a]*(Cz*f5);
	    I[111] += C[1]*W[a]*(Cz*Iy*f22);
	    I[112] += C[0]*W[a]*(Cz*f22);
	    I[113] += C[1]*W[a]*(Cy*Iz*f22);
	    I[114] += C[0]*W[a]*(Cy*f22);
	    I[115] += C[1]*W[a]*(f18*f22);
	    I[116] += C[1]*W[a]*(f18*f2);
	    I[117] += C[1]*W[a]*(Iz*f6);
	    I[118] += C[1]*W[a]*(Ix*f6);
	    I[119] += C[0]*W[a]*(f6);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[120]) {
	double T[120];
	for (int i = 0; i < 120; ++i) {
	    T[i] = I[i];
	}
	I[65] = T[0];
	I[66] = T[1];
	I[61] = T[2];
	I[119] = T[3];
	I[114] = T[4];
	I[116] = T[5];
	I[118] = T[6];
	I[117] = T[7];
	I[112] = T[8];
	I[102] = T[9];
	I[92] = T[10];
	I[82] = T[11];
	I[30] = T[12];
	I[20] = T[13];
	I[0] = T[14];
	I[10] = T[15];
	I[99] = T[16];
	I[111] = T[17];
	I[113] = T[18];
	I[110] = T[19];
	I[115] = T[20];
	I[95] = T[21];
	I[57] = T[22];
	I[55] = T[23];
	I[97] = T[24];
	I[59] = T[25];
	I[79] = T[26];
	I[39] = T[27];
	I[36] = T[28];
	I[76] = T[29];
	I[34] = T[30];
	I[74] = T[31];
	I[54] = T[32];
	I[77] = T[33];
	I[85] = T[34];
	I[103] = T[35];
	I[29] = T[36];
	I[69] = T[37];
	I[49] = T[38];
	I[60] = T[39];
	I[67] = T[40];
	I[47] = T[41];
	I[70] = T[42];
	I[40] = T[43];
	I[100] = T[44];
	I[80] = T[45];
	I[105] = T[46];
	I[89] = T[47];
	I[109] = T[48];
	I[28] = T[49];
	I[26] = T[50];
	I[38] = T[51];
	I[93] = T[52];
	I[83] = T[53];
	I[91] = T[54];
	I[81] = T[55];
	I[11] = T[56];
	I[18] = T[57];
	I[8] = T[58];
	I[31] = T[59];
	I[1] = T[60];
	I[19] = T[61];
	I[9] = T[62];
	I[64] = T[63];
	I[62] = T[64];
	I[12] = T[65];
	I[16] = T[66];
	I[6] = T[67];
	I[22] = T[68];
	I[2] = T[69];
	I[52] = T[70];
	I[42] = T[71];
	I[44] = T[72];
	I[63] = T[73];
	I[68] = T[74];
	I[94] = T[75];
	I[90] = T[76];
	I[50] = T[77];
	I[72] = T[78];
	I[32] = T[79];
	I[21] = T[80];
	I[101] = T[81];
	I[37] = T[82];
	I[78] = T[83];
	I[53] = T[84];
	I[58] = T[85];
	I[48] = T[86];
	I[73] = T[87];
	I[43] = T[88];
	I[104] = T[89];
	I[84] = T[90];
	I[96] = T[91];
	I[86] = T[92];
	I[15] = T[93];
	I[17] = T[94];
	I[27] = T[95];
	I[7] = T[96];
	I[35] = T[97];
	I[5] = T[98];
	I[25] = T[99];
	I[106] = T[100];
	I[56] = T[101];
	I[46] = T[102];
	I[75] = T[103];
	I[45] = T[104];
	I[107] = T[105];
	I[87] = T[106];
	I[98] = T[107];
	I[88] = T[108];
	I[13] = T[109];
	I[14] = T[110];
	I[24] = T[111];
	I[4] = T[112];
	I[33] = T[113];
	I[3] = T[114];
	I[23] = T[115];
	I[108] = T[116];
	I[71] = T[117];
	I[51] = T[118];
	I[41] = T[119];
    }

};

template<>
struct impl<meta::braket<rysq::SP, rysq::S, rysq::P, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[2], double (&I)[12]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;



// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);

	    I[0] += C[1]*W[a]*(Qx);
	    I[1] += C[1]*W[a]*(Qy);
	    I[2] += C[1]*W[a]*(Qz);
	    I[3] += C[1]*W[a]*(Cz*Dx);
	    I[4] += C[1]*W[a]*(Cz*Dy);
	    I[5] += C[1]*W[a]*(Cx*Dy);
	    I[6] += C[1]*W[a]*(Cx*Dz);
	    I[7] += C[1]*W[a]*(Cy*Dz);
	    I[8] += C[1]*W[a]*(Cy*Dx);
	    I[9] += C[0]*W[a]*(Dx);
	    I[10] += C[0]*W[a]*(Dy);
	    I[11] += C[0]*W[a]*(Dz);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[12]) {
	double T[12];
	for (int i = 0; i < 12; ++i) {
	    T[i] = I[i];
	}
	I[1] = T[0];
	I[6] = T[1];
	I[11] = T[2];
	I[3] = T[3];
	I[7] = T[4];
	I[5] = T[5];
	I[9] = T[6];
	I[10] = T[7];
	I[2] = T[8];
	I[0] = T[9];
	I[4] = T[10];
	I[8] = T[11];
    }

};

template<>
struct impl<meta::braket<rysq::SP, rysq::SP, rysq::SP, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[8], double (&I)[64]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qy = (Cy*Dy + B00);
	    double f2 = Cz*Dz;
	    double f4 = Cx*Dx;
	    double f8 = Cy*Dy;

	    I[0] += C[7]*W[a]*((B00*(yij + 2*Cy) + f8*yij + Dy*Py));
	    I[1] += C[7]*W[a]*((f4*xij + Dx*Px + B00*(xij + 2*Cx)));
	    I[2] += C[7]*W[a]*(Cy*(B00 + f4 + Dx*xij));
	    I[3] += C[7]*W[a]*(Cz*(B00 + f4 + Dx*xij));
	    I[4] += C[6]*W[a]*((B00 + f4 + Dx*xij));
	    I[5] += C[7]*W[a]*((Dz*Pz + B00*(2*Cz + zij) + f2*zij));
	    I[6] += C[7]*W[a]*(Cy*(B00 + Dz*zij + f2));
	    I[7] += C[7]*W[a]*(Cx*(B00 + Dz*zij + f2));
	    I[8] += C[6]*W[a]*((B00 + Dz*zij + f2));
	    I[9] += C[7]*W[a]*(Dz*(Cy*Iy + B10));
	    I[10] += C[7]*W[a]*(Dx*(Cy*Iy + B10));
	    I[11] += C[3]*W[a]*((Cy*Iy + B10));
	    I[12] += C[7]*W[a]*(Dy*(B10 + Cz*Iz));
	    I[13] += C[7]*W[a]*(Dx*(B10 + Cz*Iz));
	    I[14] += C[3]*W[a]*((B10 + Cz*Iz));
	    I[15] += C[7]*W[a]*(Cz*(Dy*yij + Qy));
	    I[16] += C[7]*W[a]*(Cx*(Dy*yij + Qy));
	    I[17] += C[6]*W[a]*((Dy*yij + Qy));
	    I[18] += C[7]*W[a]*(Dy*(Cx*Ix + B10));
	    I[19] += C[3]*W[a]*((Cx*Ix + B10));
	    I[20] += C[7]*W[a]*(Dz*(Cx*Ix + B10));
	    I[21] += C[7]*W[a]*(Cx*Dz*Iy);
	    I[22] += C[5]*W[a]*(Cx*Dz);
	    I[23] += C[5]*W[a]*(Cy*Dz);
	    I[24] += C[6]*W[a]*(Dy*Iz);
	    I[25] += C[7]*W[a]*(Iz*(B00 + f4));
	    I[26] += C[5]*W[a]*((B00 + f4));
	    I[27] += C[7]*W[a]*(Iy*(B00 + f4));
	    I[28] += C[7]*W[a]*(Iy*(B00 + f2));
	    I[29] += C[7]*W[a]*(Ix*(B00 + f2));
	    I[30] += C[5]*W[a]*((B00 + f2));
	    I[31] += C[6]*W[a]*(Dx*Iy);
	    I[32] += C[6]*W[a]*(Dx*Iz);
	    I[33] += C[7]*W[a]*(Iz*Qy);
	    I[34] += C[7]*W[a]*(Ix*Qy);
	    I[35] += C[6]*W[a]*(Dy*Ix);
	    I[36] += C[5]*W[a]*(Cx*Dy);
	    I[37] += C[7]*W[a]*(Cx*Dy*Iz);
	    I[38] += C[3]*W[a]*(Cx*Iz);
	    I[39] += C[3]*W[a]*(Cy*Iz);
	    I[40] += C[7]*W[a]*(Cy*Dx*Iz);
	    I[41] += C[5]*W[a]*(Cy*Dx);
	    I[42] += C[5]*W[a]*(Cz*Dx);
	    I[43] += C[7]*W[a]*(Cz*Dx*Iy);
	    I[44] += C[3]*W[a]*(Cz*Iy);
	    I[45] += C[5]*W[a]*(Cz*Dy);
	    I[46] += C[7]*W[a]*(Cz*Dy*Ix);
	    I[47] += C[3]*W[a]*(Cz*Ix);
	    I[48] += C[3]*W[a]*(Cy*Ix);
	    I[49] += C[7]*W[a]*(Cy*Dz*Ix);
	    I[50] += C[6]*W[a]*(Dz*Ix);
	    I[51] += C[6]*W[a]*(Dz*Iy);
	    I[52] += C[3]*W[a]*(Cx*Iy);
	    I[53] += C[1]*W[a]*(Cx);
	    I[54] += C[1]*W[a]*(Cy);
	    I[55] += C[1]*W[a]*(Cz);
	    I[56] += C[2]*W[a]*(Ix);
	    I[57] += C[2]*W[a]*(Iy);
	    I[58] += C[2]*W[a]*(Iz);
	    I[59] += C[4]*W[a]*(Dx);
	    I[60] += C[4]*W[a]*(Dy);
	    I[61] += C[5]*W[a]*(Qy);
	    I[62] += C[4]*W[a]*(Dz);
	    I[63] += C[0]*W[a]*(1);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[64]) {
	double T[64];
	for (int i = 0; i < 64; ++i) {
	    T[i] = I[i];
	}
	I[42] = T[0];
	I[21] = T[1];
	I[22] = T[2];
	I[23] = T[3];
	I[20] = T[4];
	I[63] = T[5];
	I[62] = T[6];
	I[61] = T[7];
	I[60] = T[8];
	I[58] = T[9];
	I[26] = T[10];
	I[10] = T[11];
	I[47] = T[12];
	I[31] = T[13];
	I[15] = T[14];
	I[43] = T[15];
	I[41] = T[16];
	I[40] = T[17];
	I[37] = T[18];
	I[5] = T[19];
	I[53] = T[20];
	I[57] = T[21];
	I[49] = T[22];
	I[50] = T[23];
	I[44] = T[24];
	I[29] = T[25];
	I[17] = T[26];
	I[25] = T[27];
	I[59] = T[28];
	I[55] = T[29];
	I[51] = T[30];
	I[24] = T[31];
	I[28] = T[32];
	I[46] = T[33];
	I[38] = T[34];
	I[36] = T[35];
	I[33] = T[36];
	I[45] = T[37];
	I[13] = T[38];
	I[14] = T[39];
	I[30] = T[40];
	I[18] = T[41];
	I[19] = T[42];
	I[27] = T[43];
	I[11] = T[44];
	I[35] = T[45];
	I[39] = T[46];
	I[7] = T[47];
	I[6] = T[48];
	I[54] = T[49];
	I[52] = T[50];
	I[56] = T[51];
	I[9] = T[52];
	I[1] = T[53];
	I[2] = T[54];
	I[3] = T[55];
	I[4] = T[56];
	I[8] = T[57];
	I[12] = T[58];
	I[16] = T[59];
	I[32] = T[60];
	I[34] = T[61];
	I[48] = T[62];
	I[0] = T[63];
    }

};

template<>
struct impl<meta::braket<rysq::SP, rysq::S, rysq::P, rysq::P> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[2], double (&I)[36]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;


	double xkl = rkl[0];
	double ykl = rkl[1];
	double zkl = rkl[2];

// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Kx = (xkl + Dx);
	    double Ky = (ykl + Dy);
	    double Kz = (Dz + zkl);
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f7 = (B00 + Cy*Ky);
	    double f9 = (B00 + Cx*Kx);

	    I[0] += C[1]*W[a]*((Cy*Dy*Ky + B01*Cy + B00*(ykl + 2*Dy)));
	    I[1] += C[1]*W[a]*((B01*Cx + Cx*Dx*Kx + B00*(xkl + 2*Dx)));
	    I[2] += C[1]*W[a]*(Dx*(Cz*zkl + Qz));
	    I[3] += C[1]*W[a]*(Dy*(Cz*zkl + Qz));
	    I[4] += C[1]*W[a]*(Cx*Dy*Kz);
	    I[5] += C[1]*W[a]*(Cx*(B01 + Dy*Ky));
	    I[6] += C[0]*W[a]*((B01 + Dy*Ky));
	    I[7] += C[1]*W[a]*(Cz*(B01 + Dy*Ky));
	    I[8] += C[1]*W[a]*((B01*Cz + B00*(2*Dz + zkl) + Cz*Dz*Kz));
	    I[9] += C[1]*W[a]*(Cy*(Dz*Kz + B01));
	    I[10] += C[1]*W[a]*(Cx*(Dz*Kz + B01));
	    I[11] += C[0]*W[a]*((Dz*Kz + B01));
	    I[12] += C[1]*W[a]*(Cy*(B01 + Dx*Kx));
	    I[13] += C[0]*W[a]*((B01 + Dx*Kx));
	    I[14] += C[1]*W[a]*(Cz*(B01 + Dx*Kx));
	    I[15] += C[1]*W[a]*(Cz*Dx*Ky);
	    I[16] += C[1]*W[a]*(Ky*Qx);
	    I[17] += C[1]*W[a]*(Dx*f7);
	    I[18] += C[1]*W[a]*(Kx*Qy);
	    I[19] += C[1]*W[a]*(Cy*Dz*Kx);
	    I[20] += C[0]*W[a]*(Dz*Kx);
	    I[21] += C[1]*W[a]*(Cx*Dz*Ky);
	    I[22] += C[0]*W[a]*(Dz*Ky);
	    I[23] += C[1]*W[a]*(Dz*f7);
	    I[24] += C[1]*W[a]*(Dz*f9);
	    I[25] += C[1]*W[a]*(Dy*f9);
	    I[26] += C[1]*W[a]*(Cz*Dy*Kx);
	    I[27] += C[0]*W[a]*(Dy*Kx);
	    I[28] += C[1]*W[a]*(Kx*Qz);
	    I[29] += C[1]*W[a]*(Ky*Qz);
	    I[30] += C[0]*W[a]*(Dx*Ky);
	    I[31] += C[1]*W[a]*(Cy*Dx*Kz);
	    I[32] += C[0]*W[a]*(Dx*Kz);
	    I[33] += C[1]*W[a]*(Kz*Qx);
	    I[34] += C[0]*W[a]*(Dy*Kz);
	    I[35] += C[1]*W[a]*(Kz*Qy);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[36]) {
	double T[36];
	for (int i = 0; i < 36; ++i) {
	    T[i] = I[i];
	}
	I[18] = T[0];
	I[1] = T[1];
	I[27] = T[2];
	I[31] = T[3];
	I[29] = T[4];
	I[17] = T[5];
	I[16] = T[6];
	I[19] = T[7];
	I[35] = T[8];
	I[34] = T[9];
	I[33] = T[10];
	I[32] = T[11];
	I[2] = T[12];
	I[0] = T[13];
	I[3] = T[14];
	I[15] = T[15];
	I[13] = T[16];
	I[14] = T[17];
	I[6] = T[18];
	I[10] = T[19];
	I[8] = T[20];
	I[21] = T[21];
	I[20] = T[22];
	I[22] = T[23];
	I[9] = T[24];
	I[5] = T[25];
	I[7] = T[26];
	I[4] = T[27];
	I[11] = T[28];
	I[23] = T[29];
	I[12] = T[30];
	I[26] = T[31];
	I[24] = T[32];
	I[25] = T[33];
	I[28] = T[34];
	I[30] = T[35];
    }

};

template<>
struct impl<meta::braket<rysq::P, rysq::S, rysq::P, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[1], double (&I)[9]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;



// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);

	    I[0] += C[0]*W[a]*(Qy);
	    I[1] += C[0]*W[a]*(Qz);
	    I[2] += C[0]*W[a]*(Qx);
	    I[3] += C[0]*W[a]*(Cy*Dx);
	    I[4] += C[0]*W[a]*(Cz*Dx);
	    I[5] += C[0]*W[a]*(Cz*Dy);
	    I[6] += C[0]*W[a]*(Cx*Dy);
	    I[7] += C[0]*W[a]*(Cx*Dz);
	    I[8] += C[0]*W[a]*(Cy*Dz);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[9]) {
	double T[9];
	for (int i = 0; i < 9; ++i) {
	    T[i] = I[i];
	}
	I[4] = T[0];
	I[8] = T[1];
	I[0] = T[2];
	I[1] = T[3];
	I[2] = T[4];
	I[5] = T[5];
	I[3] = T[6];
	I[6] = T[7];
	I[7] = T[8];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::SP, rysq::S, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[2], double (&I)[40]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 

// 	for (int a = 0; a < N; ++a) {
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);


	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double f0 = (3*B10 + pow(Cx,2));
	    double f11 = (3*B10 + pow(Cz,2));
	    double f12 = (3*B10 + pow(Cy,2));
	    double f13 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f2 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));
	    double f5 = (Cy*Iy + B10);
	    double f6 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f8 = 3*pow(B10,2);

	    I[0] += C[1]*W[a]*(Cx*Cy*(Cz*zij + Pz));
	    I[1] += C[1]*W[a]*(Cy*Cz*(Px + Cx*xij));
	    I[2] += C[1]*W[a]*(Py*(Px + Cx*xij));
	    I[3] += C[1]*W[a]*(Pz*(Px + Cx*xij));
	    I[4] += C[1]*W[a]*((3*B10*Cx*(xij + 2*Cx) + Ix*pow(Cx,3) + f8));
	    I[5] += C[1]*W[a]*((Iy*pow(Cy,3) + f8 + 3*B10*Cy*(yij + 2*Cy)));
	    I[6] += C[1]*W[a]*((3*B10*Cz*(2*Cz + zij) + f8 + Iz*pow(Cz,3)));
	    I[7] += C[1]*W[a]*(Px*(Cz*zij + Pz));
	    I[8] += C[1]*W[a]*(Py*(Cz*zij + Pz));
	    I[9] += C[1]*W[a]*(Cz*Iy*f11);
	    I[10] += C[1]*W[a]*(Cz*Iy*Px);
	    I[11] += C[0]*W[a]*(Cx*Cy*Cz);
	    I[12] += C[1]*W[a]*(Cx*Cz*f5);
	    I[13] += C[1]*W[a]*(Cz*Ix*Py);
	    I[14] += C[1]*W[a]*(Cz*Ix*f11);
	    I[15] += C[0]*W[a]*(Cz*f11);
	    I[16] += C[0]*W[a]*(Cz*Py);
	    I[17] += C[1]*W[a]*(Cy*f6);
	    I[18] += C[1]*W[a]*(Cz*f6);
	    I[19] += C[1]*W[a]*(Px*f5);
	    I[20] += C[1]*W[a]*(Cx*Iz*f0);
	    I[21] += C[1]*W[a]*(Cx*Iz*Py);
	    I[22] += C[0]*W[a]*(Cx*Py);
	    I[23] += C[1]*W[a]*(Cx*f2);
	    I[24] += C[1]*W[a]*(Cz*f2);
	    I[25] += C[0]*W[a]*(Cz*Px);
	    I[26] += C[0]*W[a]*(Cy*Px);
	    I[27] += C[1]*W[a]*(Cy*Iz*Px);
	    I[28] += C[1]*W[a]*(Cy*Iz*f12);
	    I[29] += C[0]*W[a]*(Cy*f12);
	    I[30] += C[1]*W[a]*(Cy*Ix*f12);
	    I[31] += C[1]*W[a]*(Cy*Ix*Pz);
	    I[32] += C[0]*W[a]*(Cy*Pz);
	    I[33] += C[1]*W[a]*(Pz*f5);
	    I[34] += C[0]*W[a]*(Cx*Pz);
	    I[35] += C[1]*W[a]*(Cx*Iy*Pz);
	    I[36] += C[1]*W[a]*(Cx*Iy*f0);
	    I[37] += C[0]*W[a]*(Cx*f0);
	    I[38] += C[1]*W[a]*(Cx*f13);
	    I[39] += C[1]*W[a]*(Cy*f13);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[40]) {
	double T[40];
	for (int i = 0; i < 40; ++i) {
	    T[i] = I[i];
	}
	I[39] = T[0];
	I[19] = T[1];
	I[15] = T[2];
	I[17] = T[3];
	I[10] = T[4];
	I[21] = T[5];
	I[32] = T[6];
	I[34] = T[7];
	I[36] = T[8];
	I[22] = T[9];
	I[24] = T[10];
	I[9] = T[11];
	I[29] = T[12];
	I[16] = T[13];
	I[12] = T[14];
	I[2] = T[15];
	I[6] = T[16];
	I[13] = T[17];
	I[14] = T[18];
	I[23] = T[19];
	I[30] = T[20];
	I[35] = T[21];
	I[5] = T[22];
	I[25] = T[23];
	I[26] = T[24];
	I[4] = T[25];
	I[3] = T[26];
	I[33] = T[27];
	I[31] = T[28];
	I[1] = T[29];
	I[11] = T[30];
	I[18] = T[31];
	I[8] = T[32];
	I[28] = T[33];
	I[7] = T[34];
	I[27] = T[35];
	I[20] = T[36];
	I[0] = T[37];
	I[37] = T[38];
	I[38] = T[39];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::SP, rysq::SP, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[4], double (&I)[96]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f10 = (Dy*Iy + B00);
	    double f15 = (Dx*Ix + B00);

	    I[0] += C[3]*W[a]*(Py*(Dz*zij + Qz));
	    I[1] += C[3]*W[a]*(Cx*Cy*(Dz*zij + Qz));
	    I[2] += C[3]*W[a]*(Qy*(Cx*Ix + B10));
	    I[3] += C[3]*W[a]*(Qz*(Cx*Ix + B10));
	    I[4] += C[3]*W[a]*(Px*(Dz*zij + Qz));
	    I[5] += C[3]*W[a]*(Cx*(Dz*(B10 + Cz*Iz) + B00*(2*Cz + zij)));
	    I[6] += C[3]*W[a]*(Cy*(Dz*(B10 + Cz*Iz) + B00*(2*Cz + zij)));
	    I[7] += C[3]*W[a]*(Cz*(B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10)));
	    I[8] += C[3]*W[a]*(Cx*(B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10)));
	    I[9] += C[3]*W[a]*(Iz*(Dy*Py + 2*B00*Cy));
	    I[10] += C[3]*W[a]*(Ix*(Dy*Py + 2*B00*Cy));
	    I[11] += C[2]*W[a]*((Dy*Py + 2*B00*Cy));
	    I[12] += C[3]*W[a]*((2*B00*Cy*yij + 3*B00*Py + Dy*(B10*(3*Cy + yij) + Iy*pow(Cy,2))));
	    I[13] += C[3]*W[a]*(Qx*(Cy*Iy + B10));
	    I[14] += C[3]*W[a]*(Qz*(Cy*Iy + B10));
	    I[15] += C[2]*W[a]*(Cy*Cz*Dx);
	    I[16] += C[3]*W[a]*(Cy*Cz*f15);
	    I[17] += C[3]*W[a]*(Cy*Iz*Qx);
	    I[18] += C[3]*W[a]*(Iz*(Dx*Px + 2*B00*Cx));
	    I[19] += C[3]*W[a]*(Iy*(Dx*Px + 2*B00*Cx));
	    I[20] += C[3]*W[a]*(Dz*Iy*Px);
	    I[21] += C[3]*W[a]*(Iy*(Dz*Pz + 2*B00*Cz));
	    I[22] += C[3]*W[a]*(Ix*(Dz*Pz + 2*B00*Cz));
	    I[23] += C[2]*W[a]*((Dz*Pz + 2*B00*Cz));
	    I[24] += C[3]*W[a]*((Dz*(Iz*pow(Cz,2) + B10*(3*Cz + zij)) + 3*B00*Pz + 2*B00*Cz*zij));
	    I[25] += C[3]*W[a]*(Qy*(B10 + Cz*Iz));
	    I[26] += C[3]*W[a]*(Qx*(B10 + Cz*Iz));
	    I[27] += C[3]*W[a]*(Dx*(Iz*pow(Cz,2) + B10*(3*Cz + zij)));
	    I[28] += C[3]*W[a]*(Cy*Dx*(B10 + Cz*Iz));
	    I[29] += C[1]*W[a]*(Cy*(B10 + Cz*Iz));
	    I[30] += C[3]*W[a]*(Cx*Dy*(B10 + Cz*Iz));
	    I[31] += C[2]*W[a]*(Cx*Cz*Dy);
	    I[32] += C[3]*W[a]*(Cz*Ix*Qy);
	    I[33] += C[3]*W[a]*(Cz*Iy*Qx);
	    I[34] += C[1]*W[a]*(Cx*Cz*Iy);
	    I[35] += C[1]*W[a]*(Iy*Pz);
	    I[36] += C[3]*W[a]*(Dx*Iy*Pz);
	    I[37] += C[2]*W[a]*(Dx*Pz);
	    I[38] += C[3]*W[a]*(Cz*(Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx)));
	    I[39] += C[3]*W[a]*(Cy*(Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx)));
	    I[40] += C[3]*W[a]*((Dx*(B10*(3*Cx + xij) + Ix*pow(Cx,2)) + 2*B00*Cx*xij + 3*B00*Px));
	    I[41] += C[2]*W[a]*((Dx*Px + 2*B00*Cx));
	    I[42] += C[3]*W[a]*(Px*f10);
	    I[43] += C[3]*W[a]*(Pz*f10);
	    I[44] += C[3]*W[a]*(Cx*Cz*f10);
	    I[45] += C[1]*W[a]*(Cx*(B10 + Cz*Iz));
	    I[46] += C[1]*W[a]*((Iz*pow(Cz,2) + B10*(3*Cz + zij)));
	    I[47] += C[3]*W[a]*(Dy*(Iz*pow(Cz,2) + B10*(3*Cz + zij)));
	    I[48] += C[2]*W[a]*(Dy*Pz);
	    I[49] += C[3]*W[a]*(Dy*Ix*Pz);
	    I[50] += C[1]*W[a]*(Ix*Pz);
	    I[51] += C[3]*W[a]*(Pz*f15);
	    I[52] += C[3]*W[a]*(Py*f15);
	    I[53] += C[1]*W[a]*(Iz*Py);
	    I[54] += C[3]*W[a]*(Dx*Iz*Py);
	    I[55] += C[2]*W[a]*(Dx*Py);
	    I[56] += C[3]*W[a]*(Dx*(B10*(3*Cy + yij) + Iy*pow(Cy,2)));
	    I[57] += C[3]*W[a]*(Cz*Dx*(Cy*Iy + B10));
	    I[58] += C[1]*W[a]*(Cz*(Cy*Iy + B10));
	    I[59] += C[3]*W[a]*(Cx*Dz*(Cy*Iy + B10));
	    I[60] += C[3]*W[a]*(Cy*Dz*(Cx*Ix + B10));
	    I[61] += C[3]*W[a]*(Dz*(B10*(3*Cx + xij) + Ix*pow(Cx,2)));
	    I[62] += C[2]*W[a]*(Dz*Px);
	    I[63] += C[2]*W[a]*(Cx*Cy*Dz);
	    I[64] += C[3]*W[a]*(Cy*Ix*Qz);
	    I[65] += C[1]*W[a]*(Cy*Cz*Ix);
	    I[66] += C[0]*W[a]*(Cy*Cz);
	    I[67] += C[2]*W[a]*(Cy*Qx);
	    I[68] += C[2]*W[a]*(Cz*Qx);
	    I[69] += C[0]*W[a]*(Cx*Cz);
	    I[70] += C[3]*W[a]*(Cx*Iz*Qy);
	    I[71] += C[2]*W[a]*(Cx*Qy);
	    I[72] += C[2]*W[a]*(Cz*Qy);
	    I[73] += C[0]*W[a]*(Pz);
	    I[74] += C[1]*W[a]*(Cz*(Cx*Ix + B10));
	    I[75] += C[3]*W[a]*(Cz*Dy*(Cx*Ix + B10));
	    I[76] += C[3]*W[a]*(Dy*(B10*(3*Cx + xij) + Ix*pow(Cx,2)));
	    I[77] += C[2]*W[a]*(Dy*Px);
	    I[78] += C[3]*W[a]*(Dy*Iz*Px);
	    I[79] += C[1]*W[a]*(Iz*Px);
	    I[80] += C[1]*W[a]*(Cx*Cy*Iz);
	    I[81] += C[0]*W[a]*(Cx*Cy);
	    I[82] += C[1]*W[a]*(Cx*(Cy*Iy + B10));
	    I[83] += C[1]*W[a]*((B10*(3*Cy + yij) + Iy*pow(Cy,2)));
	    I[84] += C[3]*W[a]*(Dz*(B10*(3*Cy + yij) + Iy*pow(Cy,2)));
	    I[85] += C[2]*W[a]*(Dz*Py);
	    I[86] += C[3]*W[a]*(Dz*Ix*Py);
	    I[87] += C[1]*W[a]*(Ix*Py);
	    I[88] += C[0]*W[a]*(Py);
	    I[89] += C[1]*W[a]*(Cy*(Cx*Ix + B10));
	    I[90] += C[1]*W[a]*((B10*(3*Cx + xij) + Ix*pow(Cx,2)));
	    I[91] += C[0]*W[a]*(Px);
	    I[92] += C[1]*W[a]*(Iy*Px);
	    I[93] += C[3]*W[a]*(Cx*Iy*Qz);
	    I[94] += C[2]*W[a]*(Cx*Qz);
	    I[95] += C[2]*W[a]*(Cy*Qz);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[96]) {
	double T[96];
	for (int i = 0; i < 96; ++i) {
	    T[i] = I[i];
	}
	I[91] = T[0];
	I[93] = T[1];
	I[57] = T[2];
	I[82] = T[3];
	I[90] = T[4];
	I[94] = T[5];
	I[95] = T[6];
	I[65] = T[7];
	I[63] = T[8];
	I[67] = T[9];
	I[55] = T[10];
	I[49] = T[11];
	I[61] = T[12];
	I[39] = T[13];
	I[89] = T[14];
	I[29] = T[15];
	I[35] = T[16];
	I[45] = T[17];
	I[42] = T[18];
	I[36] = T[19];
	I[84] = T[20];
	I[86] = T[21];
	I[80] = T[22];
	I[74] = T[23];
	I[92] = T[24];
	I[71] = T[25];
	I[46] = T[26];
	I[44] = T[27];
	I[47] = T[28];
	I[23] = T[29];
	I[70] = T[30];
	I[52] = T[31];
	I[59] = T[32];
	I[40] = T[33];
	I[16] = T[34];
	I[14] = T[35];
	I[38] = T[36];
	I[26] = T[37];
	I[34] = T[38];
	I[33] = T[39];
	I[30] = T[40];
	I[24] = T[41];
	I[60] = T[42];
	I[62] = T[43];
	I[64] = T[44];
	I[22] = T[45];
	I[20] = T[46];
	I[68] = T[47];
	I[50] = T[48];
	I[56] = T[49];
	I[8] = T[50];
	I[32] = T[51];
	I[31] = T[52];
	I[19] = T[53];
	I[43] = T[54];
	I[25] = T[55];
	I[37] = T[56];
	I[41] = T[57];
	I[17] = T[58];
	I[87] = T[59];
	I[81] = T[60];
	I[78] = T[61];
	I[72] = T[62];
	I[75] = T[63];
	I[83] = T[64];
	I[11] = T[65];
	I[5] = T[66];
	I[27] = T[67];
	I[28] = T[68];
	I[4] = T[69];
	I[69] = T[70];
	I[51] = T[71];
	I[53] = T[72];
	I[2] = T[73];
	I[10] = T[74];
	I[58] = T[75];
	I[54] = T[76];
	I[48] = T[77];
	I[66] = T[78];
	I[18] = T[79];
	I[21] = T[80];
	I[3] = T[81];
	I[15] = T[82];
	I[13] = T[83];
	I[85] = T[84];
	I[73] = T[85];
	I[79] = T[86];
	I[7] = T[87];
	I[1] = T[88];
	I[9] = T[89];
	I[6] = T[90];
	I[0] = T[91];
	I[12] = T[92];
	I[88] = T[93];
	I[76] = T[94];
	I[77] = T[95];
    }

};

template<>
struct impl<meta::braket<rysq::SP, rysq::SP, rysq::SP, rysq::P> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[8], double (&I)[192]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];

	double xkl = rkl[0];
	double ykl = rkl[1];
	double zkl = rkl[2];

// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Kx = (xkl + Dx);
	    double Ky = (ykl + Dy);
	    double Kz = (Dz + zkl);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qy = (Cy*Dy + B00);
	    double f1 = B00*zkl;
	    double f11 = 2*B00*Dx;
	    double f14 = Cz*pow(Dz,2);
	    double f19 = Cz*zkl;
	    double f2 = Cx*Dx;
	    double f25 = B01*Cy;
	    double f27 = Cx*xkl;
	    double f3 = B01*B10;
	    double f31 = (Cx*Ix + B10);
	    double f32 = 2*B00*Dy;
	    double f35 = (B01 + Dx*Kx);
	    double f36 = (Cy*Iy + B10);
	    double f38 = Cy*ykl;
	    double f41 = 2*B00*Dz;
	    double f42 = (B01 + Dy*Ky);
	    double f44 = B01*Cx;
	    double f47 = Cx*pow(Dx,2);
	    double f49 = Cz*Dz;
	    double f52 = Cy*pow(Dy,2);
	    double f53 = 2*pow(B00,2);
	    double f58 = B01*Cz;
	    double f59 = (B10 + Cz*Iz);
	    double f6 = B00*ykl;
	    double f7 = Cy*Dy;

	    I[0] += C[7]*W[a]*((f53 + Px*pow(Dx,2) + f44*(Cx + xij) + f3 + Dx*Px*xkl + B00*(2*f27 + 4*f2 + xij*(xkl + 2*Dx)) + xij*(f47 + f2*xkl)));
	    I[1] += C[7]*W[a]*(Cy*(xkl*(B00 + f2) + f11 + f47 + B01*Ix + Dx*Kx*xij));
	    I[2] += C[7]*W[a]*(Cz*(xkl*(B00 + f2) + f11 + f47 + B01*Ix + Dx*Kx*xij));
	    I[3] += C[6]*W[a]*((xkl*(B00 + f2) + f11 + f47 + B01*Ix + Dx*Kx*xij));
	    I[4] += C[7]*W[a]*(Iy*(xkl*(B00 + f2) + f11 + f44 + f47));
	    I[5] += C[7]*W[a]*(Iz*(xkl*(B00 + f2) + f11 + f44 + f47));
	    I[6] += C[5]*W[a]*((xkl*(B00 + f2) + f11 + f44 + f47));
	    I[7] += C[7]*W[a]*(Dz*(xij*(f27 + f2) + Kx*Px + B00*(xij + 2*Cx)));
	    I[8] += C[7]*W[a]*(Dy*(xij*(f27 + f2) + Kx*Px + B00*(xij + 2*Cx)));
	    I[9] += C[3]*W[a]*((xij*(f27 + f2) + Kx*Px + B00*(xij + 2*Cx)));
	    I[10] += C[7]*W[a]*(Ky*(Dx*Px + B00*(xij + 2*Cx) + f2*xij));
	    I[11] += C[7]*W[a]*(Kz*(Dx*Px + B00*(xij + 2*Cx) + f2*xij));
	    I[12] += C[7]*W[a]*((yij*(f52 + f25 + f7*ykl + f6) + f53 + Dy*Py*ykl + B01*pow(Cy,2) + f3 + Py*pow(Dy,2) + 2*Cy*f6 + 2*B00*(Dy*yij + 2*f7)));
	    I[13] += C[7]*W[a]*(Cz*(f52 + f25 + f32 + f7*ykl + yij*(B01 + Dy*Ky) + f6));
	    I[14] += C[7]*W[a]*(Cx*(f52 + f25 + f32 + f7*ykl + yij*(B01 + Dy*Ky) + f6));
	    I[15] += C[6]*W[a]*((f52 + f25 + f32 + f7*ykl + yij*(B01 + Dy*Ky) + f6));
	    I[16] += C[7]*W[a]*(Iz*(f52 + f25 + f32 + f7*ykl + f6));
	    I[17] += C[7]*W[a]*(Ix*(f52 + f25 + f32 + f7*ykl + f6));
	    I[18] += C[5]*W[a]*((f52 + f25 + f32 + f7*ykl + f6));
	    I[19] += C[7]*W[a]*(Dz*(B00*(yij + 2*Cy) + yij*(f38 + f7) + Ky*Py));
	    I[20] += C[7]*W[a]*(Dx*(B00*(yij + 2*Cy) + yij*(f38 + f7) + Ky*Py));
	    I[21] += C[3]*W[a]*((B00*(yij + 2*Cy) + yij*(f38 + f7) + Ky*Py));
	    I[22] += C[7]*W[a]*(Kx*(f7*yij + B00*(yij + 2*Cy) + Dy*Py));
	    I[23] += C[7]*W[a]*(Kz*(f7*yij + B00*(yij + 2*Cy) + Dy*Py));
	    I[24] += C[7]*W[a]*((f53 + Pz*pow(Dz,2) + Dz*Pz*zkl + f58*(Cz + zij) + zij*(f14 + f1 + Dz*f19) + f3 + 2*Cz*f1 + 2*B00*(2*f49 + Dz*zij)));
	    I[25] += C[7]*W[a]*(Cy*(f14 + Dz*(f19 + 2*B00 + Kz*zij) + f1 + B01*Iz));
	    I[26] += C[7]*W[a]*(Cx*(f14 + Dz*(f19 + 2*B00 + Kz*zij) + f1 + B01*Iz));
	    I[27] += C[6]*W[a]*((f14 + Dz*(f19 + 2*B00 + Kz*zij) + f1 + B01*Iz));
	    I[28] += C[7]*W[a]*(Iy*(f14 + f41 + f58 + f1 + Dz*f19));
	    I[29] += C[7]*W[a]*(Ix*(f14 + f41 + f58 + f1 + Dz*f19));
	    I[30] += C[5]*W[a]*((f14 + f41 + f58 + f1 + Dz*f19));
	    I[31] += C[7]*W[a]*(Dy*(Kz*Pz + zij*(f19 + f49) + B00*(2*Cz + zij)));
	    I[32] += C[7]*W[a]*(Dx*(Kz*Pz + zij*(f19 + f49) + B00*(2*Cz + zij)));
	    I[33] += C[3]*W[a]*((Kz*Pz + zij*(f19 + f49) + B00*(2*Cz + zij)));
	    I[34] += C[7]*W[a]*(Ky*(f49*zij + Dz*Pz + B00*(2*Cz + zij)));
	    I[35] += C[7]*W[a]*(Kx*(f49*zij + Dz*Pz + B00*(2*Cz + zij)));
	    I[36] += C[7]*W[a]*((f38 + Qy)*(f49 + B00 + Dz*zij));
	    I[37] += C[7]*W[a]*(Qy*(f19 + f49 + B00 + Kz*zij));
	    I[38] += C[7]*W[a]*((B00 + f2)*(f19 + f49 + B00 + Kz*zij));
	    I[39] += C[6]*W[a]*(Dx*(f19 + f49 + B00 + Kz*zij));
	    I[40] += C[7]*W[a]*(Cy*Dx*(f19 + f49 + B00 + Kz*zij));
	    I[41] += C[3]*W[a]*(Cy*(f19 + f49 + B00 + Kz*zij));
	    I[42] += C[6]*W[a]*(Dy*(f19 + f49 + B00 + Kz*zij));
	    I[43] += C[7]*W[a]*(Cx*Dy*(f19 + f49 + B00 + Kz*zij));
	    I[44] += C[3]*W[a]*(Cx*(f19 + f49 + B00 + Kz*zij));
	    I[45] += C[2]*W[a]*((f19 + f49 + B00 + Kz*zij));
	    I[46] += C[7]*W[a]*((f27 + B00 + f2)*(f49 + B00 + Dz*zij));
	    I[47] += C[7]*W[a]*((Dy*yij + Qy)*(f27 + B00 + f2));
	    I[48] += C[7]*W[a]*((Dy*yij + Qy)*(f19 + f49 + B00));
	    I[49] += C[7]*W[a]*((f49 + B00)*(f38 + Ky*yij + Qy));
	    I[50] += C[7]*W[a]*((B00 + f2)*(f38 + Ky*yij + Qy));
	    I[51] += C[6]*W[a]*(Dx*(f38 + Ky*yij + Qy));
	    I[52] += C[7]*W[a]*(Cz*Dx*(f38 + Ky*yij + Qy));
	    I[53] += C[3]*W[a]*(Cz*(f38 + Ky*yij + Qy));
	    I[54] += C[6]*W[a]*(Dz*(f38 + Ky*yij + Qy));
	    I[55] += C[7]*W[a]*(Cx*Dz*(f38 + Ky*yij + Qy));
	    I[56] += C[3]*W[a]*(Cx*(f38 + Ky*yij + Qy));
	    I[57] += C[2]*W[a]*((f38 + Ky*yij + Qy));
	    I[58] += C[7]*W[a]*(Cx*Kz*(Dy*yij + Qy));
	    I[59] += C[7]*W[a]*(Cz*Kx*(Dy*yij + Qy));
	    I[60] += C[6]*W[a]*(Kx*(Dy*yij + Qy));
	    I[61] += C[6]*W[a]*(Kz*(Dy*yij + Qy));
	    I[62] += C[7]*W[a]*((f38 + Qy)*(B00 + f2 + Dx*xij));
	    I[63] += C[7]*W[a]*((B00 + f2 + Dx*xij)*(f19 + f49 + B00));
	    I[64] += C[7]*W[a]*((f49 + B00)*(f27 + Kx*xij + B00 + f2));
	    I[65] += C[7]*W[a]*(Cy*Dz*(f27 + Kx*xij + B00 + f2));
	    I[66] += C[6]*W[a]*(Dz*(f27 + Kx*xij + B00 + f2));
	    I[67] += C[7]*W[a]*(Qy*(f27 + Kx*xij + B00 + f2));
	    I[68] += C[6]*W[a]*(Dy*(f27 + Kx*xij + B00 + f2));
	    I[69] += C[7]*W[a]*(Cz*Dy*(f27 + Kx*xij + B00 + f2));
	    I[70] += C[3]*W[a]*(Cz*(f27 + Kx*xij + B00 + f2));
	    I[71] += C[2]*W[a]*((f27 + Kx*xij + B00 + f2));
	    I[72] += C[3]*W[a]*(Cy*(f27 + Kx*xij + B00 + f2));
	    I[73] += C[7]*W[a]*(Cy*Kz*(B00 + f2 + Dx*xij));
	    I[74] += C[6]*W[a]*(Kz*(B00 + f2 + Dx*xij));
	    I[75] += C[7]*W[a]*(Cz*Ky*(B00 + f2 + Dx*xij));
	    I[76] += C[6]*W[a]*(Ky*(B00 + f2 + Dx*xij));
	    I[77] += C[7]*W[a]*(Dx*Ky*f59);
	    I[78] += C[7]*W[a]*(Cz*Ix*f42);
	    I[79] += C[7]*W[a]*(Ix*Ky*(f49 + B00));
	    I[80] += C[7]*W[a]*(Cx*Ky*(f49 + B00 + Dz*zij));
	    I[81] += C[6]*W[a]*(Ky*(f49 + B00 + Dz*zij));
	    I[82] += C[7]*W[a]*(Cy*Kx*(f49 + B00 + Dz*zij));
	    I[83] += C[6]*W[a]*(Kx*(f49 + B00 + Dz*zij));
	    I[84] += C[7]*W[a]*(Iy*Kx*(f49 + B00));
	    I[85] += C[5]*W[a]*(Kx*(f49 + B00));
	    I[86] += C[5]*W[a]*(Dx*(f19 + f49 + B00));
	    I[87] += C[7]*W[a]*(Dx*Iy*(f19 + f49 + B00));
	    I[88] += C[3]*W[a]*(Iy*(f19 + f49 + B00));
	    I[89] += C[5]*W[a]*(Dy*(f19 + f49 + B00));
	    I[90] += C[7]*W[a]*(Dy*Ix*(f19 + f49 + B00));
	    I[91] += C[3]*W[a]*(Ix*(f19 + f49 + B00));
	    I[92] += C[1]*W[a]*((f19 + f49 + B00));
	    I[93] += C[5]*W[a]*(Ky*(f49 + B00));
	    I[94] += C[7]*W[a]*(Iz*Ky*(B00 + f2));
	    I[95] += C[5]*W[a]*(Ky*(B00 + f2));
	    I[96] += C[7]*W[a]*(Dz*Iy*(f27 + B00 + f2));
	    I[97] += C[5]*W[a]*(Dz*(f27 + B00 + f2));
	    I[98] += C[5]*W[a]*(Dy*(f27 + B00 + f2));
	    I[99] += C[7]*W[a]*(Dy*Iz*(f27 + B00 + f2));
	    I[100] += C[3]*W[a]*(Iz*(f27 + B00 + f2));
	    I[101] += C[1]*W[a]*((f27 + B00 + f2));
	    I[102] += C[3]*W[a]*(Iy*(f27 + B00 + f2));
	    I[103] += C[7]*W[a]*(Iy*Kz*(B00 + f2));
	    I[104] += C[5]*W[a]*(Kz*(B00 + f2));
	    I[105] += C[7]*W[a]*(Ix*Kz*Qy);
	    I[106] += C[7]*W[a]*(Dx*Kz*f36);
	    I[107] += C[7]*W[a]*(f36*(Dz*Kz + B01));
	    I[108] += C[7]*W[a]*(f31*(Dz*Kz + B01));
	    I[109] += C[6]*W[a]*(Ix*(Dz*Kz + B01));
	    I[110] += C[7]*W[a]*(Cy*Ix*(Dz*Kz + B01));
	    I[111] += C[5]*W[a]*(Cy*(Dz*Kz + B01));
	    I[112] += C[6]*W[a]*(Iy*(Dz*Kz + B01));
	    I[113] += C[7]*W[a]*(Cx*Iy*(Dz*Kz + B01));
	    I[114] += C[5]*W[a]*(Cx*(Dz*Kz + B01));
	    I[115] += C[4]*W[a]*((Dz*Kz + B01));
	    I[116] += C[5]*W[a]*(Dz*(f38 + Qy));
	    I[117] += C[7]*W[a]*(Dz*Ix*(f38 + Qy));
	    I[118] += C[3]*W[a]*(Ix*(f38 + Qy));
	    I[119] += C[5]*W[a]*(Dx*(f38 + Qy));
	    I[120] += C[7]*W[a]*(Dx*Iz*(f38 + Qy));
	    I[121] += C[3]*W[a]*(Iz*(f38 + Qy));
	    I[122] += C[7]*W[a]*(Iz*Kx*Qy);
	    I[123] += C[5]*W[a]*(Kx*Qy);
	    I[124] += C[1]*W[a]*((f38 + Qy));
	    I[125] += C[7]*W[a]*(Cx*Iz*f42);
	    I[126] += C[5]*W[a]*(Cx*f42);
	    I[127] += C[5]*W[a]*(Cz*f42);
	    I[128] += C[6]*W[a]*(Ix*f42);
	    I[129] += C[7]*W[a]*(f31*f42);
	    I[130] += C[6]*W[a]*(Iz*f42);
	    I[131] += C[7]*W[a]*(f42*f59);
	    I[132] += C[3]*W[a]*(Ky*f59);
	    I[133] += C[7]*W[a]*(Dz*Ky*f31);
	    I[134] += C[6]*W[a]*(Dz*Ix*Ky);
	    I[135] += C[2]*W[a]*(Ix*Ky);
	    I[136] += C[3]*W[a]*(Cz*Ix*Ky);
	    I[137] += C[1]*W[a]*(Cz*Ky);
	    I[138] += C[5]*W[a]*(Cz*Dx*Ky);
	    I[139] += C[4]*W[a]*(Dx*Ky);
	    I[140] += C[6]*W[a]*(Dx*Iz*Ky);
	    I[141] += C[2]*W[a]*(Iz*Ky);
	    I[142] += C[3]*W[a]*(Cx*Iz*Ky);
	    I[143] += C[1]*W[a]*(Cx*Ky);
	    I[144] += C[5]*W[a]*(Cx*Dz*Ky);
	    I[145] += C[4]*W[a]*(Dz*Ky);
	    I[146] += C[3]*W[a]*(Ky*f31);
	    I[147] += C[7]*W[a]*(Dy*Kz*f31);
	    I[148] += C[5]*W[a]*(Cx*Dy*Kz);
	    I[149] += C[6]*W[a]*(Dy*Ix*Kz);
	    I[150] += C[3]*W[a]*(Cy*Ix*Kz);
	    I[151] += C[5]*W[a]*(Cy*Dx*Kz);
	    I[152] += C[6]*W[a]*(Dx*Iy*Kz);
	    I[153] += C[3]*W[a]*(Cx*Iy*Kz);
	    I[154] += C[1]*W[a]*(Cx*Kz);
	    I[155] += C[1]*W[a]*(Cy*Kz);
	    I[156] += C[2]*W[a]*(Ix*Kz);
	    I[157] += C[3]*W[a]*(Kz*f31);
	    I[158] += C[2]*W[a]*(Iy*Kz);
	    I[159] += C[4]*W[a]*(Dx*Kz);
	    I[160] += C[4]*W[a]*(Dy*Kz);
	    I[161] += C[5]*W[a]*(Kz*Qy);
	    I[162] += C[3]*W[a]*(Kz*f36);
	    I[163] += C[7]*W[a]*(f35*f36);
	    I[164] += C[6]*W[a]*(Iy*f35);
	    I[165] += C[7]*W[a]*(Cz*Iy*f35);
	    I[166] += C[5]*W[a]*(Cz*f35);
	    I[167] += C[6]*W[a]*(Iz*f35);
	    I[168] += C[7]*W[a]*(Cy*Iz*f35);
	    I[169] += C[5]*W[a]*(Cy*f35);
	    I[170] += C[7]*W[a]*(f35*f59);
	    I[171] += C[7]*W[a]*(Dy*Kx*f59);
	    I[172] += C[3]*W[a]*(Kx*f59);
	    I[173] += C[3]*W[a]*(Kx*f36);
	    I[174] += C[7]*W[a]*(Dz*Kx*f36);
	    I[175] += C[5]*W[a]*(Cy*Dz*Kx);
	    I[176] += C[4]*W[a]*(Dz*Kx);
	    I[177] += C[6]*W[a]*(Dz*Iy*Kx);
	    I[178] += C[2]*W[a]*(Iy*Kx);
	    I[179] += C[3]*W[a]*(Cz*Iy*Kx);
	    I[180] += C[1]*W[a]*(Cz*Kx);
	    I[181] += C[5]*W[a]*(Cz*Dy*Kx);
	    I[182] += C[4]*W[a]*(Dy*Kx);
	    I[183] += C[6]*W[a]*(Dy*Iz*Kx);
	    I[184] += C[2]*W[a]*(Iz*Kx);
	    I[185] += C[3]*W[a]*(Cy*Iz*Kx);
	    I[186] += C[1]*W[a]*(Cy*Kx);
	    I[187] += C[0]*W[a]*(Kx);
	    I[188] += C[4]*W[a]*(f35);
	    I[189] += C[0]*W[a]*(Ky);
	    I[190] += C[4]*W[a]*(f42);
	    I[191] += C[0]*W[a]*(Kz);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[192]) {
	double T[192];
	for (int i = 0; i < 192; ++i) {
	    T[i] = I[i];
	}
	I[21] = T[0];
	I[22] = T[1];
	I[23] = T[2];
	I[20] = T[3];
	I[25] = T[4];
	I[29] = T[5];
	I[17] = T[6];
	I[53] = T[7];
	I[37] = T[8];
	I[5] = T[9];
	I[85] = T[10];
	I[149] = T[11];
	I[106] = T[12];
	I[107] = T[13];
	I[105] = T[14];
	I[104] = T[15];
	I[110] = T[16];
	I[102] = T[17];
	I[98] = T[18];
	I[122] = T[19];
	I[90] = T[20];
	I[74] = T[21];
	I[42] = T[22];
	I[170] = T[23];
	I[191] = T[24];
	I[190] = T[25];
	I[189] = T[26];
	I[188] = T[27];
	I[187] = T[28];
	I[183] = T[29];
	I[179] = T[30];
	I[175] = T[31];
	I[159] = T[32];
	I[143] = T[33];
	I[127] = T[34];
	I[63] = T[35];
	I[126] = T[36];
	I[174] = T[37];
	I[157] = T[38];
	I[156] = T[39];
	I[158] = T[40];
	I[142] = T[41];
	I[172] = T[42];
	I[173] = T[43];
	I[141] = T[44];
	I[140] = T[45];
	I[61] = T[46];
	I[41] = T[47];
	I[171] = T[48];
	I[123] = T[49];
	I[89] = T[50];
	I[88] = T[51];
	I[91] = T[52];
	I[75] = T[53];
	I[120] = T[54];
	I[121] = T[55];
	I[73] = T[56];
	I[72] = T[57];
	I[169] = T[58];
	I[43] = T[59];
	I[40] = T[60];
	I[168] = T[61];
	I[86] = T[62];
	I[151] = T[63];
	I[55] = T[64];
	I[54] = T[65];
	I[52] = T[66];
	I[38] = T[67];
	I[36] = T[68];
	I[39] = T[69];
	I[7] = T[70];
	I[4] = T[71];
	I[6] = T[72];
	I[150] = T[73];
	I[148] = T[74];
	I[87] = T[75];
	I[84] = T[76];
	I[95] = T[77];
	I[103] = T[78];
	I[119] = T[79];
	I[125] = T[80];
	I[124] = T[81];
	I[62] = T[82];
	I[60] = T[83];
	I[59] = T[84];
	I[51] = T[85];
	I[147] = T[86];
	I[155] = T[87];
	I[139] = T[88];
	I[163] = T[89];
	I[167] = T[90];
	I[135] = T[91];
	I[131] = T[92];
	I[115] = T[93];
	I[93] = T[94];
	I[81] = T[95];
	I[57] = T[96];
	I[49] = T[97];
	I[33] = T[98];
	I[45] = T[99];
	I[13] = T[100];
	I[1] = T[101];
	I[9] = T[102];
	I[153] = T[103];
	I[145] = T[104];
	I[166] = T[105];
	I[154] = T[106];
	I[186] = T[107];
	I[181] = T[108];
	I[180] = T[109];
	I[182] = T[110];
	I[178] = T[111];
	I[184] = T[112];
	I[185] = T[113];
	I[177] = T[114];
	I[176] = T[115];
	I[114] = T[116];
	I[118] = T[117];
	I[70] = T[118];
	I[82] = T[119];
	I[94] = T[120];
	I[78] = T[121];
	I[46] = T[122];
	I[34] = T[123];
	I[66] = T[124];
	I[109] = T[125];
	I[97] = T[126];
	I[99] = T[127];
	I[100] = T[128];
	I[101] = T[129];
	I[108] = T[130];
	I[111] = T[131];
	I[79] = T[132];
	I[117] = T[133];
	I[116] = T[134];
	I[68] = T[135];
	I[71] = T[136];
	I[67] = T[137];
	I[83] = T[138];
	I[80] = T[139];
	I[92] = T[140];
	I[76] = T[141];
	I[77] = T[142];
	I[65] = T[143];
	I[113] = T[144];
	I[112] = T[145];
	I[69] = T[146];
	I[165] = T[147];
	I[161] = T[148];
	I[164] = T[149];
	I[134] = T[150];
	I[146] = T[151];
	I[152] = T[152];
	I[137] = T[153];
	I[129] = T[154];
	I[130] = T[155];
	I[132] = T[156];
	I[133] = T[157];
	I[136] = T[158];
	I[144] = T[159];
	I[160] = T[160];
	I[162] = T[161];
	I[138] = T[162];
	I[26] = T[163];
	I[24] = T[164];
	I[27] = T[165];
	I[19] = T[166];
	I[28] = T[167];
	I[30] = T[168];
	I[18] = T[169];
	I[31] = T[170];
	I[47] = T[171];
	I[15] = T[172];
	I[10] = T[173];
	I[58] = T[174];
	I[50] = T[175];
	I[48] = T[176];
	I[56] = T[177];
	I[8] = T[178];
	I[11] = T[179];
	I[3] = T[180];
	I[35] = T[181];
	I[32] = T[182];
	I[44] = T[183];
	I[12] = T[184];
	I[14] = T[185];
	I[2] = T[186];
	I[0] = T[187];
	I[16] = T[188];
	I[64] = T[189];
	I[96] = T[190];
	I[128] = T[191];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::S, rysq::SP, rysq::P> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[2], double (&I)[120]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;


	double xkl = rkl[0];
	double ykl = rkl[1];
	double zkl = rkl[2];

// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Kx = (xkl + Dx);
	    double Ky = (ykl + Dy);
	    double Kz = (Dz + zkl);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f0 = (3*B10 + pow(Cx,2));
	    double f1 = (Cx*Kx*(3*B10 + pow(Cx,2)) + 3*B00*Px);
	    double f10 = (2*B00*Cx*(xkl + 2*Dx) + 2*pow(B00,2) + Px*(B01 + Dx*Kx));
	    double f11 = (B01 + Dx*Kx);
	    double f14 = (Kx*Px + 2*B00*Cx);
	    double f16 = (Dy*Py + 2*B00*Cy);
	    double f17 = (2*pow(B00,2) + Pz*(Dz*Kz + B01) + 2*B00*Cz*(2*Dz + zkl));
	    double f2 = (Dz*Pz + 2*B00*Cz);
	    double f21 = (Dx*Px + 2*B00*Cx);
	    double f22 = (Dz*Kz + B01);
	    double f23 = (B01 + Dy*Ky);
	    double f25 = (B00 + Cx*Kx);
	    double f3 = (Py*(B01 + Dy*Ky) + 2*pow(B00,2) + 2*B00*Cy*(ykl + 2*Dy));
	    double f31 = (3*B10 + pow(Cz,2));
	    double f32 = (3*B10 + pow(Cy,2));
	    double f36 = (B01*Cz + B00*(2*Dz + zkl) + Cz*Dz*Kz);
	    double f4 = (Cy*Dy*Ky + B01*Cy + B00*(ykl + 2*Dy));
	    double f5 = (3*B00*Py + Cy*Ky*(3*B10 + pow(Cy,2)));
	    double f6 = (Cy*Dy*(3*B10 + pow(Cy,2)) + 3*B00*Py);
	    double f7 = 3*B00*B10;
	    double f8 = (B01*Cx + Cx*Dx*Kx + B00*(xkl + 2*Dx));

	    I[0] += C[1]*W[a]*((Dy*(2*f7 + Dy*pow(Cy,3) + 3*B10*Cy*Dy) + 6*Cy*pow(B00,2) + ykl*(f7 + Dy*pow(Cy,3) + 3*B10*Cy*Dy) + B01*(3*B10*Cy + pow(Cy,3)) + 3*B00*pow(Cy,2)*(ykl + 2*Dy)));
	    I[1] += C[1]*W[a]*((xkl*(Dx*pow(Cx,3) + 3*B10*Cx*Dx + f7) + Dx*(Dx*pow(Cx,3) + 3*B10*Cx*Dx + 2*f7) + 3*B00*pow(Cx,2)*(xkl + 2*Dx) + 6*Cx*pow(B00,2) + B01*(3*B10*Cx + pow(Cx,3))));
	    I[2] += C[1]*W[a]*(Kz*(Dx*pow(Cx,3) + 3*B00*pow(Cx,2) + 3*B10*Cx*Dx + f7));
	    I[3] += C[1]*W[a]*(Ky*(Dx*pow(Cx,3) + 3*B00*pow(Cx,2) + 3*B10*Cx*Dx + f7));
	    I[4] += C[1]*W[a]*(Ky*(3*B10*Cz*Dz + f7 + 3*B00*pow(Cz,2) + Dz*pow(Cz,3)));
	    I[5] += C[1]*W[a]*(Kx*(3*B10*Cz*Dz + f7 + 3*B00*pow(Cz,2) + Dz*pow(Cz,3)));
	    I[6] += C[1]*W[a]*((Dz*(3*B10*Cz*Dz + 2*f7 + Dz*pow(Cz,3)) + 3*B00*pow(Cz,2)*(2*Dz + zkl) + zkl*(3*B10*Cz*Dz + f7 + Dz*pow(Cz,3)) + B01*(3*B10*Cz + pow(Cz,3)) + 6*Cz*pow(B00,2)));
	    I[7] += C[1]*W[a]*(Dy*(Cz*Kz*f31 + f7 + 3*B00*pow(Cz,2)));
	    I[8] += C[0]*W[a]*((Cz*Kz*f31 + f7 + 3*B00*pow(Cz,2)));
	    I[9] += C[1]*W[a]*(Dx*(Cz*Kz*f31 + f7 + 3*B00*pow(Cz,2)));
	    I[10] += C[1]*W[a]*(Cy*Dx*(Pz*zkl + f2));
	    I[11] += C[1]*W[a]*(Cx*Dy*(Pz*zkl + f2));
	    I[12] += C[0]*W[a]*(Cx*(Pz*zkl + f2));
	    I[13] += C[1]*W[a]*(Qx*(Pz*zkl + f2));
	    I[14] += C[1]*W[a]*(Qy*(Pz*zkl + f2));
	    I[15] += C[0]*W[a]*(Cy*(Pz*zkl + f2));
	    I[16] += C[1]*W[a]*(Qz*(f16 + Py*ykl));
	    I[17] += C[1]*W[a]*(Qx*(f16 + Py*ykl));
	    I[18] += C[1]*W[a]*(Cz*Dx*(f16 + Py*ykl));
	    I[19] += C[0]*W[a]*(Cz*(f16 + Py*ykl));
	    I[20] += C[1]*W[a]*(Cx*Dz*(f16 + Py*ykl));
	    I[21] += C[0]*W[a]*(Cx*(f16 + Py*ykl));
	    I[22] += C[1]*W[a]*(Cx*Qz*(Cy*ykl + Qy));
	    I[23] += C[0]*W[a]*(Cx*Cz*(Cy*ykl + Qy));
	    I[24] += C[1]*W[a]*(Cz*Qx*(Cy*ykl + Qy));
	    I[25] += C[1]*W[a]*(Cy*Qx*(Cz*zkl + Qz));
	    I[26] += C[0]*W[a]*(Cx*Cy*(Cz*zkl + Qz));
	    I[27] += C[1]*W[a]*(Cx*Qy*(Cz*zkl + Qz));
	    I[28] += C[1]*W[a]*(Ky*Px*Qz);
	    I[29] += C[1]*W[a]*(Cx*Ky*f2);
	    I[30] += C[1]*W[a]*(Kz*Py*Qx);
	    I[31] += C[1]*W[a]*(Kz*Px*Qy);
	    I[32] += C[1]*W[a]*(Cx*Kz*f16);
	    I[33] += C[1]*W[a]*(Cx*Pz*f23);
	    I[34] += C[1]*W[a]*(Cx*f0*f23);
	    I[35] += C[1]*W[a]*(Cx*f0*f22);
	    I[36] += C[1]*W[a]*(Cx*Cz*f4);
	    I[37] += C[1]*W[a]*(Cz*Qy*f25);
	    I[38] += C[1]*W[a]*(Cy*Qz*f25);
	    I[39] += C[1]*W[a]*(Kx*Py*Qz);
	    I[40] += C[1]*W[a]*(Cy*Kx*f2);
	    I[41] += C[1]*W[a]*(f2*(Cy*ykl + Qy));
	    I[42] += C[1]*W[a]*(Dx*Pz*(Cy*ykl + Qy));
	    I[43] += C[0]*W[a]*(Pz*(Cy*ykl + Qy));
	    I[44] += C[1]*W[a]*(Dz*Px*(Cy*ykl + Qy));
	    I[45] += C[0]*W[a]*(Px*(Cy*ykl + Qy));
	    I[46] += C[1]*W[a]*(f21*(Cy*ykl + Qy));
	    I[47] += C[1]*W[a]*(Cy*Kz*f21);
	    I[48] += C[0]*W[a]*(Cy*Kz*Px);
	    I[49] += C[1]*W[a]*(Cy*Dx*Kz*f32);
	    I[50] += C[0]*W[a]*(Cy*Kz*f32);
	    I[51] += C[1]*W[a]*(Cy*f22*f32);
	    I[52] += C[1]*W[a]*(Cy*Px*f22);
	    I[53] += C[1]*W[a]*(Cy*Cz*f8);
	    I[54] += C[0]*W[a]*(Cy*Cz*f25);
	    I[55] += C[1]*W[a]*(Cz*Kx*f16);
	    I[56] += C[1]*W[a]*(f16*(Cz*zkl + Qz));
	    I[57] += C[1]*W[a]*(Dx*Py*(Cz*zkl + Qz));
	    I[58] += C[0]*W[a]*(Py*(Cz*zkl + Qz));
	    I[59] += C[1]*W[a]*(Dy*Px*(Cz*zkl + Qz));
	    I[60] += C[0]*W[a]*(Px*(Cz*zkl + Qz));
	    I[61] += C[1]*W[a]*(f21*(Cz*zkl + Qz));
	    I[62] += C[1]*W[a]*(Cz*Ky*f21);
	    I[63] += C[0]*W[a]*(Cz*Ky*Px);
	    I[64] += C[1]*W[a]*(Cz*Px*f23);
	    I[65] += C[1]*W[a]*(Cz*f23*f31);
	    I[66] += C[1]*W[a]*(Cz*f11*f31);
	    I[67] += C[1]*W[a]*(Cz*Py*f11);
	    I[68] += C[0]*W[a]*(Cz*Kx*Py);
	    I[69] += C[1]*W[a]*(Cz*Dy*Kx*f31);
	    I[70] += C[0]*W[a]*(Cz*Kx*f31);
	    I[71] += C[1]*W[a]*(Cz*Dx*Ky*f31);
	    I[72] += C[0]*W[a]*(Cz*Ky*f31);
	    I[73] += C[1]*W[a]*(Ky*Pz*Qx);
	    I[74] += C[0]*W[a]*(Cx*Ky*Pz);
	    I[75] += C[1]*W[a]*(Cx*Dz*Ky*f0);
	    I[76] += C[0]*W[a]*(Cx*Ky*f0);
	    I[77] += C[1]*W[a]*(Cx*Dy*Kz*f0);
	    I[78] += C[0]*W[a]*(Cx*Kz*f0);
	    I[79] += C[0]*W[a]*(Cx*Kz*Py);
	    I[80] += C[1]*W[a]*(Cx*Py*f22);
	    I[81] += C[1]*W[a]*(Cx*Cy*f36);
	    I[82] += C[1]*W[a]*(Cy*Pz*f11);
	    I[83] += C[1]*W[a]*(Cy*f11*f32);
	    I[84] += C[1]*W[a]*(Cy*Dz*Kx*f32);
	    I[85] += C[0]*W[a]*(Cy*Kx*f32);
	    I[86] += C[0]*W[a]*(Cy*Kx*Pz);
	    I[87] += C[1]*W[a]*(Kx*Pz*Qy);
	    I[88] += C[1]*W[a]*(Qy*f14);
	    I[89] += C[1]*W[a]*(Qz*f14);
	    I[90] += C[1]*W[a]*(f16*f25);
	    I[91] += C[1]*W[a]*(f2*f25);
	    I[92] += C[1]*W[a]*(Kx*f6);
	    I[93] += C[1]*W[a]*(Kz*f6);
	    I[94] += C[1]*W[a]*(Pz*f8);
	    I[95] += C[1]*W[a]*(Py*f8);
	    I[96] += C[1]*W[a]*(Dz*Py*f25);
	    I[97] += C[0]*W[a]*(Py*f25);
	    I[98] += C[1]*W[a]*(Dy*Pz*f25);
	    I[99] += C[0]*W[a]*(Pz*f25);
	    I[100] += C[1]*W[a]*(Pz*f4);
	    I[101] += C[1]*W[a]*(Px*f4);
	    I[102] += C[1]*W[a]*(Px*f36);
	    I[103] += C[1]*W[a]*(Py*f36);
	    I[104] += C[1]*W[a]*(Cz*f10);
	    I[105] += C[1]*W[a]*(Cy*f10);
	    I[106] += C[1]*W[a]*(Cy*Dz*f14);
	    I[107] += C[0]*W[a]*(Cy*f14);
	    I[108] += C[1]*W[a]*(Cz*Dy*f14);
	    I[109] += C[0]*W[a]*(Cz*f14);
	    I[110] += C[1]*W[a]*(Cz*f3);
	    I[111] += C[1]*W[a]*(Cx*f3);
	    I[112] += C[1]*W[a]*(Cx*f17);
	    I[113] += C[1]*W[a]*(Cy*f17);
	    I[114] += C[1]*W[a]*(Dx*f5);
	    I[115] += C[1]*W[a]*(Dz*f5);
	    I[116] += C[1]*W[a]*(Dz*f1);
	    I[117] += C[1]*W[a]*(Dy*f1);
	    I[118] += C[0]*W[a]*(f1);
	    I[119] += C[0]*W[a]*(f5);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[120]) {
	double T[120];
	for (int i = 0; i < 120; ++i) {
	    T[i] = I[i];
	}
	I[61] = T[0];
	I[10] = T[1];
	I[90] = T[2];
	I[50] = T[3];
	I[72] = T[4];
	I[32] = T[5];
	I[112] = T[6];
	I[102] = T[7];
	I[82] = T[8];
	I[92] = T[9];
	I[98] = T[10];
	I[107] = T[11];
	I[87] = T[12];
	I[97] = T[13];
	I[108] = T[14];
	I[88] = T[15];
	I[76] = T[16];
	I[55] = T[17];
	I[56] = T[18];
	I[46] = T[19];
	I[75] = T[20];
	I[45] = T[21];
	I[79] = T[22];
	I[49] = T[23];
	I[59] = T[24];
	I[99] = T[25];
	I[89] = T[26];
	I[109] = T[27];
	I[74] = T[28];
	I[77] = T[29];
	I[95] = T[30];
	I[103] = T[31];
	I[105] = T[32];
	I[67] = T[33];
	I[60] = T[34];
	I[110] = T[35];
	I[69] = T[36];
	I[29] = T[37];
	I[39] = T[38];
	I[36] = T[39];
	I[38] = T[40];
	I[78] = T[41];
	I[58] = T[42];
	I[48] = T[43];
	I[73] = T[44];
	I[43] = T[45];
	I[53] = T[46];
	I[93] = T[47];
	I[83] = T[48];
	I[91] = T[49];
	I[81] = T[50];
	I[111] = T[51];
	I[113] = T[52];
	I[19] = T[53];
	I[9] = T[54];
	I[26] = T[55];
	I[106] = T[56];
	I[96] = T[57];
	I[86] = T[58];
	I[104] = T[59];
	I[84] = T[60];
	I[94] = T[61];
	I[54] = T[62];
	I[44] = T[63];
	I[64] = T[64];
	I[62] = T[65];
	I[12] = T[66];
	I[16] = T[67];
	I[6] = T[68];
	I[22] = T[69];
	I[2] = T[70];
	I[52] = T[71];
	I[42] = T[72];
	I[57] = T[73];
	I[47] = T[74];
	I[70] = T[75];
	I[40] = T[76];
	I[100] = T[77];
	I[80] = T[78];
	I[85] = T[79];
	I[115] = T[80];
	I[119] = T[81];
	I[18] = T[82];
	I[11] = T[83];
	I[31] = T[84];
	I[1] = T[85];
	I[8] = T[86];
	I[28] = T[87];
	I[23] = T[88];
	I[34] = T[89];
	I[25] = T[90];
	I[37] = T[91];
	I[21] = T[92];
	I[101] = T[93];
	I[17] = T[94];
	I[15] = T[95];
	I[35] = T[96];
	I[5] = T[97];
	I[27] = T[98];
	I[7] = T[99];
	I[68] = T[100];
	I[63] = T[101];
	I[114] = T[102];
	I[116] = T[103];
	I[14] = T[104];
	I[13] = T[105];
	I[33] = T[106];
	I[3] = T[107];
	I[24] = T[108];
	I[4] = T[109];
	I[66] = T[110];
	I[65] = T[111];
	I[117] = T[112];
	I[118] = T[113];
	I[51] = T[114];
	I[71] = T[115];
	I[30] = T[116];
	I[20] = T[117];
	I[0] = T[118];
	I[41] = T[119];
    }

};

template<>
struct impl<meta::braket<rysq::SP, rysq::SP, rysq::P, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[4], double (&I)[48]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qy = (Cy*Dy + B00);
	    double f1 = Cz*Dz;
	    double f14 = (B10 + Cz*Iz);
	    double f15 = (Cy*Iy + B10);
	    double f4 = Cx*Dx;
	    double f8 = Cy*Dy;
	    double f9 = (Cx*Ix + B10);

	    I[0] += C[3]*W[a]*((B00*(yij + 2*Cy) + f8*yij + Dy*Py));
	    I[1] += C[3]*W[a]*((f4*xij + Dx*Px + B00*(xij + 2*Cx)));
	    I[2] += C[3]*W[a]*(Cy*(B00 + f4 + Dx*xij));
	    I[3] += C[3]*W[a]*(Cz*(B00 + f4 + Dx*xij));
	    I[4] += C[2]*W[a]*((B00 + f4 + Dx*xij));
	    I[5] += C[3]*W[a]*((Dz*Pz + B00*(2*Cz + zij) + f1*zij));
	    I[6] += C[3]*W[a]*(Cy*(B00 + Dz*zij + f1));
	    I[7] += C[3]*W[a]*(Cx*(B00 + Dz*zij + f1));
	    I[8] += C[2]*W[a]*((B00 + Dz*zij + f1));
	    I[9] += C[3]*W[a]*(Cz*(Dy*yij + Qy));
	    I[10] += C[3]*W[a]*(Cx*(Dy*yij + Qy));
	    I[11] += C[2]*W[a]*((Dy*yij + Qy));
	    I[12] += C[3]*W[a]*(Ix*Qy);
	    I[13] += C[3]*W[a]*(Iz*Qy);
	    I[14] += C[3]*W[a]*(Iz*(B00 + f4));
	    I[15] += C[1]*W[a]*((B00 + f4));
	    I[16] += C[3]*W[a]*(Iy*(B00 + f4));
	    I[17] += C[3]*W[a]*(Iy*(B00 + f1));
	    I[18] += C[3]*W[a]*(Ix*(B00 + f1));
	    I[19] += C[1]*W[a]*((B00 + f1));
	    I[20] += C[3]*W[a]*(Dx*f14);
	    I[21] += C[3]*W[a]*(Dy*f14);
	    I[22] += C[3]*W[a]*(Dy*f9);
	    I[23] += C[2]*W[a]*(Dy*Ix);
	    I[24] += C[3]*W[a]*(Cz*Dy*Ix);
	    I[25] += C[1]*W[a]*(Cz*Dy);
	    I[26] += C[2]*W[a]*(Dy*Iz);
	    I[27] += C[3]*W[a]*(Cx*Dy*Iz);
	    I[28] += C[1]*W[a]*(Cx*Dy);
	    I[29] += C[3]*W[a]*(Cx*Dz*Iy);
	    I[30] += C[1]*W[a]*(Cx*Dz);
	    I[31] += C[1]*W[a]*(Cy*Dz);
	    I[32] += C[3]*W[a]*(Cy*Dz*Ix);
	    I[33] += C[2]*W[a]*(Dz*Ix);
	    I[34] += C[3]*W[a]*(Dz*f9);
	    I[35] += C[2]*W[a]*(Dz*Iy);
	    I[36] += C[3]*W[a]*(Dz*f15);
	    I[37] += C[3]*W[a]*(Dx*f15);
	    I[38] += C[2]*W[a]*(Dx*Iy);
	    I[39] += C[3]*W[a]*(Cz*Dx*Iy);
	    I[40] += C[1]*W[a]*(Cz*Dx);
	    I[41] += C[2]*W[a]*(Dx*Iz);
	    I[42] += C[3]*W[a]*(Cy*Dx*Iz);
	    I[43] += C[1]*W[a]*(Cy*Dx);
	    I[44] += C[0]*W[a]*(Dx);
	    I[45] += C[0]*W[a]*(Dy);
	    I[46] += C[1]*W[a]*(Qy);
	    I[47] += C[0]*W[a]*(Dz);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[48]) {
	double T[48];
	for (int i = 0; i < 48; ++i) {
	    T[i] = I[i];
	}
	I[26] = T[0];
	I[5] = T[1];
	I[6] = T[2];
	I[7] = T[3];
	I[4] = T[4];
	I[47] = T[5];
	I[46] = T[6];
	I[45] = T[7];
	I[44] = T[8];
	I[27] = T[9];
	I[25] = T[10];
	I[24] = T[11];
	I[22] = T[12];
	I[30] = T[13];
	I[13] = T[14];
	I[1] = T[15];
	I[9] = T[16];
	I[43] = T[17];
	I[39] = T[18];
	I[35] = T[19];
	I[15] = T[20];
	I[31] = T[21];
	I[21] = T[22];
	I[20] = T[23];
	I[23] = T[24];
	I[19] = T[25];
	I[28] = T[26];
	I[29] = T[27];
	I[17] = T[28];
	I[41] = T[29];
	I[33] = T[30];
	I[34] = T[31];
	I[38] = T[32];
	I[36] = T[33];
	I[37] = T[34];
	I[40] = T[35];
	I[42] = T[36];
	I[10] = T[37];
	I[8] = T[38];
	I[11] = T[39];
	I[3] = T[40];
	I[12] = T[41];
	I[14] = T[42];
	I[2] = T[43];
	I[0] = T[44];
	I[16] = T[45];
	I[18] = T[46];
	I[32] = T[47];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::P, rysq::SP, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[2], double (&I)[72]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f1 = (B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10));
	    double f10 = (Cy*Iy + B10);
	    double f11 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f13 = (Dx*Px + 2*B00*Cx);
	    double f14 = (Dx*Ix + B00);
	    double f21 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f24 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));
	    double f3 = 3*B00*B10;
	    double f4 = (Dy*Py + 2*B00*Cy);
	    double f5 = (Dz*Pz + 2*B00*Cz);
	    double f7 = (Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx));
	    double f9 = (Dy*Iy + B00);

	    I[0] += C[1]*W[a]*((Dx*Ix*pow(Cx,2) + B00*Cx*(3*Cx + 2*xij) + B10*Dx*(3*Cx + xij) + f3));
	    I[1] += C[1]*W[a]*((Dy*Iy*pow(Cy,2) + f3 + B00*Cy*(3*Cy + 2*yij) + B10*Dy*(3*Cy + yij)));
	    I[2] += C[1]*W[a]*((B00*Cz*(3*Cz + 2*zij) + Dz*Iz*pow(Cz,2) + B10*Dz*(3*Cz + zij) + f3));
	    I[3] += C[1]*W[a]*(Cx*(f5 + Qz*zij));
	    I[4] += C[1]*W[a]*(Cy*(f5 + Qz*zij));
	    I[5] += C[1]*W[a]*(Cx*Cy*(Dz*zij + Qz));
	    I[6] += C[1]*W[a]*(Py*(Dz*zij + Qz));
	    I[7] += C[1]*W[a]*(Px*(Dz*zij + Qz));
	    I[8] += C[1]*W[a]*(Qz*(Px + Cx*xij));
	    I[9] += C[1]*W[a]*(Qy*(Px + Cx*xij));
	    I[10] += C[1]*W[a]*(Cy*Dz*(Px + Cx*xij));
	    I[11] += C[0]*W[a]*(Cy*(Px + Cx*xij));
	    I[12] += C[1]*W[a]*(Cz*Dy*(Px + Cx*xij));
	    I[13] += C[0]*W[a]*(Cz*(Px + Cx*xij));
	    I[14] += C[1]*W[a]*(Cx*Cz*f9);
	    I[15] += C[1]*W[a]*(Dy*Iz*Px);
	    I[16] += C[1]*W[a]*(Cy*Iz*Qx);
	    I[17] += C[0]*W[a]*(Cx*Cy*Iz);
	    I[18] += C[1]*W[a]*(Cx*Iz*Qy);
	    I[19] += C[1]*W[a]*(Cz*Ix*Qy);
	    I[20] += C[1]*W[a]*(Qy*(Cz*zij + Pz));
	    I[21] += C[1]*W[a]*(Qx*(Cz*zij + Pz));
	    I[22] += C[1]*W[a]*(Cx*Dy*(Cz*zij + Pz));
	    I[23] += C[0]*W[a]*(Cx*(Cz*zij + Pz));
	    I[24] += C[1]*W[a]*(Cy*Dx*(Cz*zij + Pz));
	    I[25] += C[0]*W[a]*(Cy*(Cz*zij + Pz));
	    I[26] += C[1]*W[a]*(Cy*Cz*f14);
	    I[27] += C[0]*W[a]*(Cy*Cz*Ix);
	    I[28] += C[1]*W[a]*(Cy*Ix*Qz);
	    I[29] += C[1]*W[a]*(Cx*Iy*Qz);
	    I[30] += C[0]*W[a]*(Cx*Cz*Iy);
	    I[31] += C[1]*W[a]*(Cz*Iy*Qx);
	    I[32] += C[1]*W[a]*(Qx*f10);
	    I[33] += C[1]*W[a]*(Iy*f13);
	    I[34] += C[1]*W[a]*(Iz*f13);
	    I[35] += C[1]*W[a]*(Ix*f4);
	    I[36] += C[1]*W[a]*(Cx*f1);
	    I[37] += C[1]*W[a]*(Cz*f1);
	    I[38] += C[1]*W[a]*(Px*f9);
	    I[39] += C[1]*W[a]*(Pz*f9);
	    I[40] += C[1]*W[a]*(Pz*f14);
	    I[41] += C[1]*W[a]*(Py*f14);
	    I[42] += C[1]*W[a]*(Dx*Iz*Py);
	    I[43] += C[0]*W[a]*(Iz*Py);
	    I[44] += C[1]*W[a]*(Iz*f4);
	    I[45] += C[0]*W[a]*(Iz*Px);
	    I[46] += C[1]*W[a]*(Dz*Iy*Px);
	    I[47] += C[0]*W[a]*(Iy*Px);
	    I[48] += C[1]*W[a]*(Dx*Iy*Pz);
	    I[49] += C[0]*W[a]*(Iy*Pz);
	    I[50] += C[1]*W[a]*(Dy*Ix*Pz);
	    I[51] += C[0]*W[a]*(Ix*Pz);
	    I[52] += C[1]*W[a]*(Dz*Ix*Py);
	    I[53] += C[0]*W[a]*(Ix*Py);
	    I[54] += C[1]*W[a]*(Ix*f5);
	    I[55] += C[1]*W[a]*(Iy*f5);
	    I[56] += C[1]*W[a]*(Cy*f7);
	    I[57] += C[1]*W[a]*(Cz*f7);
	    I[58] += C[1]*W[a]*(Cz*Dx*f10);
	    I[59] += C[0]*W[a]*(Cz*f10);
	    I[60] += C[1]*W[a]*(Cx*Dz*f10);
	    I[61] += C[0]*W[a]*(Cx*f10);
	    I[62] += C[1]*W[a]*(Qz*f10);
	    I[63] += C[1]*W[a]*(Dz*f11);
	    I[64] += C[1]*W[a]*(Dz*f24);
	    I[65] += C[1]*W[a]*(Dx*f24);
	    I[66] += C[1]*W[a]*(Dx*f21);
	    I[67] += C[1]*W[a]*(Dy*f21);
	    I[68] += C[1]*W[a]*(Dy*f11);
	    I[69] += C[0]*W[a]*(f11);
	    I[70] += C[0]*W[a]*(f24);
	    I[71] += C[0]*W[a]*(f21);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[72]) {
	double T[72];
	for (int i = 0; i < 72; ++i) {
	    T[i] = I[i];
	}
	I[18] = T[0];
	I[43] = T[1];
	I[68] = T[2];
	I[70] = T[3];
	I[71] = T[4];
	I[69] = T[5];
	I[67] = T[6];
	I[66] = T[7];
	I[58] = T[8];
	I[39] = T[9];
	I[57] = T[10];
	I[3] = T[11];
	I[40] = T[12];
	I[4] = T[13];
	I[46] = T[14];
	I[48] = T[15];
	I[33] = T[16];
	I[15] = T[17];
	I[51] = T[18];
	I[41] = T[19];
	I[53] = T[20];
	I[34] = T[21];
	I[52] = T[22];
	I[16] = T[23];
	I[35] = T[24];
	I[17] = T[25];
	I[23] = T[26];
	I[5] = T[27];
	I[59] = T[28];
	I[64] = T[29];
	I[10] = T[30];
	I[28] = T[31];
	I[27] = T[32];
	I[24] = T[33];
	I[30] = T[34];
	I[37] = T[35];
	I[45] = T[36];
	I[47] = T[37];
	I[42] = T[38];
	I[44] = T[39];
	I[20] = T[40];
	I[19] = T[41];
	I[31] = T[42];
	I[13] = T[43];
	I[49] = T[44];
	I[12] = T[45];
	I[60] = T[46];
	I[6] = T[47];
	I[26] = T[48];
	I[8] = T[49];
	I[38] = T[50];
	I[2] = T[51];
	I[55] = T[52];
	I[1] = T[53];
	I[56] = T[54];
	I[62] = T[55];
	I[21] = T[56];
	I[22] = T[57];
	I[29] = T[58];
	I[11] = T[59];
	I[63] = T[60];
	I[9] = T[61];
	I[65] = T[62];
	I[54] = T[63];
	I[61] = T[64];
	I[25] = T[65];
	I[32] = T[66];
	I[50] = T[67];
	I[36] = T[68];
	I[0] = T[69];
	I[7] = T[70];
	I[14] = T[71];
    }

};

template<>
struct impl<meta::braket<rysq::SP, rysq::SP, rysq::S, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[4], double (&I)[16]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 

// 	for (int a = 0; a < N; ++a) {
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);


	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);

	    I[0] += C[3]*W[a]*((Cy*Iy + B10));
	    I[1] += C[3]*W[a]*((B10 + Cz*Iz));
	    I[2] += C[3]*W[a]*((Cx*Ix + B10));
	    I[3] += C[3]*W[a]*(Cx*Iz);
	    I[4] += C[3]*W[a]*(Cy*Iz);
	    I[5] += C[3]*W[a]*(Cy*Ix);
	    I[6] += C[3]*W[a]*(Cz*Ix);
	    I[7] += C[3]*W[a]*(Cz*Iy);
	    I[8] += C[3]*W[a]*(Cx*Iy);
	    I[9] += C[1]*W[a]*(Cx);
	    I[10] += C[1]*W[a]*(Cy);
	    I[11] += C[1]*W[a]*(Cz);
	    I[12] += C[2]*W[a]*(Ix);
	    I[13] += C[2]*W[a]*(Iy);
	    I[14] += C[2]*W[a]*(Iz);
	    I[15] += C[0]*W[a]*(1);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[16]) {
	double T[16];
	for (int i = 0; i < 16; ++i) {
	    T[i] = I[i];
	}
	I[10] = T[0];
	I[15] = T[1];
	I[5] = T[2];
	I[13] = T[3];
	I[14] = T[4];
	I[6] = T[5];
	I[7] = T[6];
	I[11] = T[7];
	I[9] = T[8];
	I[1] = T[9];
	I[2] = T[10];
	I[3] = T[11];
	I[4] = T[12];
	I[8] = T[13];
	I[12] = T[14];
	I[0] = T[15];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::D, rysq::SP, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[2], double (&I)[144]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f1 = (Dz*Pz + 2*B00*Cz);
	    double f11 = (Dz*(Iz*pow(Cz,2) + B10*(3*Cz + zij)) + 3*B00*Pz + 2*B00*Cz*zij);
	    double f12 = (B00 + Dz*Iz);
	    double f15 = (B00*(3*B10 + Ix*(3*Cx + xij)) + Dx*(B10*(3*Cx + 2*xij) + Cx*pow(Ix,2)));
	    double f17 = (Dx*(B10*(3*Cx + xij) + Ix*pow(Cx,2)) + 2*B00*Cx*xij + 3*B00*Px);
	    double f18 = (B10 + pow(Iy,2));
	    double f2 = (B10*(3*Cx + 2*xij) + Cx*pow(Ix,2));
	    double f22 = (2*B00*Cy*yij + 3*B00*Py + Dy*(B10*(3*Cy + yij) + Iy*pow(Cy,2)));
	    double f25 = (B00*(Iz*(3*Cz + zij) + 3*B10) + Dz*(B10*(3*Cz + 2*zij) + Cz*pow(Iz,2)));
	    double f27 = (B10*(3*Cz + 2*zij) + Cz*pow(Iz,2));
	    double f28 = (Cy*Iy + B10);
	    double f29 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f3 = (3*pow(B10,2) + pow(Cz,2)*pow(Iz,2) + B10*(6*Cz*zij + pow(zij,2) + 6*pow(Cz,2)));
	    double f30 = (B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10));
	    double f35 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));
	    double f36 = (Dy*Py + 2*B00*Cy);
	    double f37 = (Dy*(Cy*pow(Iy,2) + B10*(3*Cy + 2*yij)) + B00*(3*B10 + Iy*(3*Cy + yij)));
	    double f38 = 3*pow(B10,2);
	    double f4 = (3*pow(B10,2) + B10*(6*Cx*xij + pow(xij,2) + 6*pow(Cx,2)) + pow(Cx,2)*pow(Ix,2));
	    double f41 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f5 = (3*pow(B10,2) + B10*(6*Cy*yij + pow(yij,2) + 6*pow(Cy,2)) + pow(Cy,2)*pow(Iy,2));
	    double f6 = (2*B00*Iy + Dy*(B10 + pow(Iy,2)));
	    double f7 = (Dx*Px + 2*B00*Cx);
	    double f8 = (Cy*pow(Iy,2) + B10*(3*Cy + 2*yij));

	    I[0] += C[1]*W[a]*((Px + Cx*xij)*(f1 + Qz*zij));
	    I[1] += C[1]*W[a]*(Cy*Ix*(f1 + Qz*zij));
	    I[2] += C[1]*W[a]*(Cx*Cy*(zij*(2*B00 + Dz*(2*Cz + zij)) + f1));
	    I[3] += C[1]*W[a]*(Cx*Iy*(f1 + Qz*zij));
	    I[4] += C[1]*W[a]*((2*B00*(xij + 2*Cx)*(Cx*Ix + 3*B10) + Dx*(f38 + B10*(6*Cx*xij + pow(xij,2) + 6*pow(Cx,2)) + pow(Cx,2)*pow(Ix,2))));
	    I[5] += C[1]*W[a]*((2*B00*(yij + 2*Cy)*(3*B10 + Cy*Iy) + Dy*(B10*(6*Cy*yij + pow(yij,2) + 6*pow(Cy,2)) + f38 + pow(Cy,2)*pow(Iy,2))));
	    I[6] += C[1]*W[a]*((2*B00*(2*Cz + zij)*(3*B10 + Cz*Iz) + Dz*(f38 + pow(Cz,2)*pow(Iz,2) + B10*(6*Cz*zij + pow(zij,2) + 6*pow(Cz,2)))));
	    I[7] += C[1]*W[a]*(Py*(zij*(2*B00 + Dz*(2*Cz + zij)) + f1));
	    I[8] += C[1]*W[a]*(Px*(zij*(2*B00 + Dz*(2*Cz + zij)) + f1));
	    I[9] += C[1]*W[a]*(f28*(f1 + Qz*zij));
	    I[10] += C[1]*W[a]*(Cy*(Cz*zij + Pz)*(Dx*xij + Qx));
	    I[11] += C[1]*W[a]*(Iz*Py*(Dx*xij + Qx));
	    I[12] += C[1]*W[a]*(Cy*Iz*(Qx*xij + f7));
	    I[13] += C[1]*W[a]*(Cz*Iy*(Qx*xij + f7));
	    I[14] += C[1]*W[a]*(Cy*Cz*(xij*(2*B00 + Dx*(xij + 2*Cx)) + f7));
	    I[15] += C[1]*W[a]*(f28*(Qx*xij + f7));
	    I[16] += C[1]*W[a]*(Py*(xij*(2*B00 + Dx*(xij + 2*Cx)) + f7));
	    I[17] += C[1]*W[a]*(Pz*(xij*(2*B00 + Dx*(xij + 2*Cx)) + f7));
	    I[18] += C[1]*W[a]*((Cz*zij + Pz)*(Qx*xij + f7));
	    I[19] += C[1]*W[a]*(Cx*(Cz*zij + Pz)*(Dy*yij + Qy));
	    I[20] += C[1]*W[a]*(Ix*Pz*(Dy*yij + Qy));
	    I[21] += C[1]*W[a]*(Iz*Px*(Dy*yij + Qy));
	    I[22] += C[1]*W[a]*(f29*(Dy*yij + Qy));
	    I[23] += C[1]*W[a]*(f41*(Dy*yij + Qy));
	    I[24] += C[1]*W[a]*(Cz*(Px + Cx*xij)*(Dy*yij + Qy));
	    I[25] += C[1]*W[a]*(Dy*(Px + Cx*xij)*(Cz*zij + Pz));
	    I[26] += C[0]*W[a]*((Px + Cx*xij)*(Cz*zij + Pz));
	    I[27] += C[1]*W[a]*(Cy*f12*(Px + Cx*xij));
	    I[28] += C[1]*W[a]*(Iy*Qz*(Px + Cx*xij));
	    I[29] += C[1]*W[a]*(Cy*Qz*(xij*(xij + 2*Cx) + Px));
	    I[30] += C[0]*W[a]*(Cy*Iz*(Px + Cx*xij));
	    I[31] += C[1]*W[a]*(Iz*Qy*(Px + Cx*xij));
	    I[32] += C[1]*W[a]*(Cz*Qy*(xij*(xij + 2*Cx) + Px));
	    I[33] += C[0]*W[a]*(Cz*Iy*(Px + Cx*xij));
	    I[34] += C[0]*W[a]*(Cy*Cz*(xij*(xij + 2*Cx) + Px));
	    I[35] += C[1]*W[a]*(f30*(Px + Cx*xij));
	    I[36] += C[1]*W[a]*(f36*(xij*(xij + 2*Cx) + Px));
	    I[37] += C[1]*W[a]*(Dy*Pz*(xij*(xij + 2*Cx) + Px));
	    I[38] += C[0]*W[a]*(Pz*(xij*(xij + 2*Cx) + Px));
	    I[39] += C[0]*W[a]*(f28*(Px + Cx*xij));
	    I[40] += C[1]*W[a]*(Dz*f28*(Px + Cx*xij));
	    I[41] += C[1]*W[a]*(Dz*Py*(xij*(xij + 2*Cx) + Px));
	    I[42] += C[0]*W[a]*(Py*(xij*(xij + 2*Cx) + Px));
	    I[43] += C[1]*W[a]*(f1*(xij*(xij + 2*Cx) + Px));
	    I[44] += C[1]*W[a]*(Cx*Cz*f6);
	    I[45] += C[1]*W[a]*(Ix*Iz*f36);
	    I[46] += C[1]*W[a]*(Cx*Iz*f30);
	    I[47] += C[1]*W[a]*(Cx*Dz*f8);
	    I[48] += C[1]*W[a]*(Cx*Qz*f18);
	    I[49] += C[0]*W[a]*(Cx*Cz*f18);
	    I[50] += C[1]*W[a]*(Cz*Qx*f18);
	    I[51] += C[1]*W[a]*(Cz*f28*(Dx*xij + Qx));
	    I[52] += C[1]*W[a]*(Iy*Pz*(Dx*xij + Qx));
	    I[53] += C[1]*W[a]*(f41*(Dx*xij + Qx));
	    I[54] += C[1]*W[a]*(f35*(Dx*xij + Qx));
	    I[55] += C[1]*W[a]*(Dx*Iz*f35);
	    I[56] += C[1]*W[a]*(Dz*Ix*f35);
	    I[57] += C[1]*W[a]*(Ix*Qz*f28);
	    I[58] += C[1]*W[a]*(Ix*Py*f12);
	    I[59] += C[0]*W[a]*(Ix*Iz*Py);
	    I[60] += C[1]*W[a]*(Iy*Iz*f7);
	    I[61] += C[0]*W[a]*(Iy*Iz*Px);
	    I[62] += C[1]*W[a]*(Iy*Px*f12);
	    I[63] += C[1]*W[a]*(Ix*Iy*f1);
	    I[64] += C[0]*W[a]*(Ix*Iy*Pz);
	    I[65] += C[1]*W[a]*(Ix*Qy*(Cz*zij + Pz));
	    I[66] += C[1]*W[a]*(Cx*Qy*(zij*(2*Cz + zij) + Pz));
	    I[67] += C[0]*W[a]*(Cx*Iy*(Cz*zij + Pz));
	    I[68] += C[1]*W[a]*(Iy*Qx*(Cz*zij + Pz));
	    I[69] += C[1]*W[a]*(Cy*Qx*(zij*(2*Cz + zij) + Pz));
	    I[70] += C[0]*W[a]*(Cy*Ix*(Cz*zij + Pz));
	    I[71] += C[0]*W[a]*(Cx*Cy*(zij*(2*Cz + zij) + Pz));
	    I[72] += C[1]*W[a]*(f36*(zij*(2*Cz + zij) + Pz));
	    I[73] += C[1]*W[a]*(f7*(zij*(2*Cz + zij) + Pz));
	    I[74] += C[0]*W[a]*(f28*(Cz*zij + Pz));
	    I[75] += C[1]*W[a]*(Dx*f28*(Cz*zij + Pz));
	    I[76] += C[1]*W[a]*(Dx*Py*(zij*(2*Cz + zij) + Pz));
	    I[77] += C[0]*W[a]*(Py*(zij*(2*Cz + zij) + Pz));
	    I[78] += C[1]*W[a]*(Dy*Px*(zij*(2*Cz + zij) + Pz));
	    I[79] += C[0]*W[a]*(Px*(zij*(2*Cz + zij) + Pz));
	    I[80] += C[1]*W[a]*(f30*(Cz*zij + Pz));
	    I[81] += C[1]*W[a]*(Cz*Ix*f30);
	    I[82] += C[0]*W[a]*(Cz*Ix*f28);
	    I[83] += C[1]*W[a]*(Iz*Qx*f28);
	    I[84] += C[0]*W[a]*(Cx*Iz*f28);
	    I[85] += C[1]*W[a]*(Cx*f12*f28);
	    I[86] += C[1]*W[a]*(Cy*f15);
	    I[87] += C[1]*W[a]*(Cz*f15);
	    I[88] += C[1]*W[a]*(Px*f6);
	    I[89] += C[1]*W[a]*(Pz*f6);
	    I[90] += C[1]*W[a]*(Cx*f37);
	    I[91] += C[1]*W[a]*(Cz*f37);
	    I[92] += C[1]*W[a]*(Qy*f27);
	    I[93] += C[1]*W[a]*(Qy*f2);
	    I[94] += C[1]*W[a]*(Qz*f2);
	    I[95] += C[1]*W[a]*(f18*f7);
	    I[96] += C[1]*W[a]*(Dx*Pz*f18);
	    I[97] += C[0]*W[a]*(Pz*f18);
	    I[98] += C[1]*W[a]*(Dz*Px*f18);
	    I[99] += C[0]*W[a]*(Px*f18);
	    I[100] += C[1]*W[a]*(f1*f18);
	    I[101] += C[1]*W[a]*(Qx*f8);
	    I[102] += C[1]*W[a]*(Qx*f27);
	    I[103] += C[1]*W[a]*(Cx*Dy*f27);
	    I[104] += C[0]*W[a]*(Cx*f27);
	    I[105] += C[1]*W[a]*(Cy*Dx*f27);
	    I[106] += C[0]*W[a]*(Cy*f27);
	    I[107] += C[1]*W[a]*(Cy*Dz*f2);
	    I[108] += C[0]*W[a]*(Cy*f2);
	    I[109] += C[1]*W[a]*(Cz*Dy*f2);
	    I[110] += C[0]*W[a]*(Cz*f2);
	    I[111] += C[1]*W[a]*(Cz*Dx*f8);
	    I[112] += C[0]*W[a]*(Cz*f8);
	    I[113] += C[1]*W[a]*(Qz*f8);
	    I[114] += C[0]*W[a]*(Cx*f8);
	    I[115] += C[1]*W[a]*(Cx*f25);
	    I[116] += C[1]*W[a]*(Cy*f25);
	    I[117] += C[1]*W[a]*(f12*f29);
	    I[118] += C[1]*W[a]*(Ix*f22);
	    I[119] += C[1]*W[a]*(Iz*f22);
	    I[120] += C[1]*W[a]*(Iz*f17);
	    I[121] += C[1]*W[a]*(Iy*f17);
	    I[122] += C[1]*W[a]*(Dx*Iy*f41);
	    I[123] += C[0]*W[a]*(Iy*f41);
	    I[124] += C[1]*W[a]*(Dy*Ix*f41);
	    I[125] += C[0]*W[a]*(Ix*f41);
	    I[126] += C[1]*W[a]*(Ix*f11);
	    I[127] += C[0]*W[a]*(Ix*f35);
	    I[128] += C[1]*W[a]*(f12*f35);
	    I[129] += C[0]*W[a]*(Iz*f35);
	    I[130] += C[1]*W[a]*(Dy*Iz*f29);
	    I[131] += C[0]*W[a]*(Iz*f29);
	    I[132] += C[1]*W[a]*(Dz*Iy*f29);
	    I[133] += C[0]*W[a]*(Iy*f29);
	    I[134] += C[1]*W[a]*(Iy*f11);
	    I[135] += C[1]*W[a]*(Dz*f4);
	    I[136] += C[1]*W[a]*(Dz*f5);
	    I[137] += C[1]*W[a]*(Dx*f5);
	    I[138] += C[1]*W[a]*(Dx*f3);
	    I[139] += C[1]*W[a]*(Dy*f3);
	    I[140] += C[1]*W[a]*(Dy*f4);
	    I[141] += C[0]*W[a]*(f4);
	    I[142] += C[0]*W[a]*(f5);
	    I[143] += C[0]*W[a]*(f3);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[144]) {
	double T[144];
	for (int i = 0; i < 144; ++i) {
	    T[i] = I[i];
	}
	I[136] = T[0];
	I[137] = T[1];
	I[123] = T[2];
	I[142] = T[3];
	I[36] = T[4];
	I[79] = T[5];
	I[122] = T[6];
	I[121] = T[7];
	I[120] = T[8];
	I[143] = T[9];
	I[65] = T[10];
	I[61] = T[11];
	I[63] = T[12];
	I[58] = T[13];
	I[41] = T[14];
	I[57] = T[15];
	I[37] = T[16];
	I[38] = T[17];
	I[64] = T[18];
	I[106] = T[19];
	I[92] = T[20];
	I[102] = T[21];
	I[90] = T[22];
	I[104] = T[23];
	I[94] = T[24];
	I[100] = T[25];
	I[28] = T[26];
	I[135] = T[27];
	I[130] = T[28];
	I[113] = T[29];
	I[27] = T[30];
	I[99] = T[31];
	I[77] = T[32];
	I[22] = T[33];
	I[5] = T[34];
	I[93] = T[35];
	I[73] = T[36];
	I[74] = T[37];
	I[2] = T[38];
	I[21] = T[39];
	I[129] = T[40];
	I[109] = T[41];
	I[1] = T[42];
	I[110] = T[43];
	I[82] = T[44];
	I[97] = T[45];
	I[105] = T[46];
	I[117] = T[47];
	I[118] = T[48];
	I[10] = T[49];
	I[46] = T[50];
	I[59] = T[51];
	I[56] = T[52];
	I[62] = T[53];
	I[55] = T[54];
	I[67] = T[55];
	I[127] = T[56];
	I[131] = T[57];
	I[133] = T[58];
	I[25] = T[59];
	I[66] = T[60];
	I[30] = T[61];
	I[138] = T[62];
	I[128] = T[63];
	I[20] = T[64];
	I[101] = T[65];
	I[87] = T[66];
	I[34] = T[67];
	I[70] = T[68];
	I[51] = T[69];
	I[29] = T[70];
	I[15] = T[71];
	I[85] = T[72];
	I[48] = T[73];
	I[35] = T[74];
	I[71] = T[75];
	I[49] = T[76];
	I[13] = T[77];
	I[84] = T[78];
	I[12] = T[79];
	I[107] = T[80];
	I[95] = T[81];
	I[23] = T[82];
	I[69] = T[83];
	I[33] = T[84];
	I[141] = T[85];
	I[39] = T[86];
	I[40] = T[87];
	I[78] = T[88];
	I[80] = T[89];
	I[81] = T[90];
	I[83] = T[91];
	I[89] = T[92];
	I[75] = T[93];
	I[112] = T[94];
	I[42] = T[95];
	I[44] = T[96];
	I[8] = T[97];
	I[114] = T[98];
	I[6] = T[99];
	I[116] = T[100];
	I[45] = T[101];
	I[52] = T[102];
	I[88] = T[103];
	I[16] = T[104];
	I[53] = T[105];
	I[17] = T[106];
	I[111] = T[107];
	I[3] = T[108];
	I[76] = T[109];
	I[4] = T[110];
	I[47] = T[111];
	I[11] = T[112];
	I[119] = T[113];
	I[9] = T[114];
	I[124] = T[115];
	I[125] = T[116];
	I[132] = T[117];
	I[91] = T[118];
	I[103] = T[119];
	I[60] = T[120];
	I[54] = T[121];
	I[68] = T[122];
	I[32] = T[123];
	I[98] = T[124];
	I[26] = T[125];
	I[134] = T[126];
	I[19] = T[127];
	I[139] = T[128];
	I[31] = T[129];
	I[96] = T[130];
	I[24] = T[131];
	I[126] = T[132];
	I[18] = T[133];
	I[140] = T[134];
	I[108] = T[135];
	I[115] = T[136];
	I[43] = T[137];
	I[50] = T[138];
	I[86] = T[139];
	I[72] = T[140];
	I[0] = T[141];
	I[7] = T[142];
	I[14] = T[143];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::S, rysq::P, rysq::P> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[1], double (&I)[54]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;


	double xkl = rkl[0];
	double ykl = rkl[1];
	double zkl = rkl[2];

// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Kx = (xkl + Dx);
	    double Ky = (ykl + Dy);
	    double Kz = (Dz + zkl);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f0 = (Dy*Py + 2*B00*Cy);
	    double f1 = (Dz*Pz + 2*B00*Cz);
	    double f10 = B01*B10;
	    double f13 = (Dx*Px + 2*B00*Cx);
	    double f15 = (Dz*Kz + B01);
	    double f16 = (B00 + Cx*Kx);
	    double f17 = (B00 + Cy*Ky);
	    double f2 = (Cy*Dy*Ky + B01*Cy + B00*(ykl + 2*Dy));
	    double f20 = 2*pow(B00,2);
	    double f22 = (B01 + Dy*Ky);
	    double f4 = (B01*Cx + Cx*Dx*Kx + B00*(xkl + 2*Dx));
	    double f6 = (B01*Cz + B00*(2*Dz + zkl) + Cz*Dz*Kz);
	    double f7 = (B01 + Dx*Kx);
	    double f9 = (Kx*Px + 2*B00*Cx);

	    I[0] += C[0]*W[a]*((B01*pow(Cx,2) + f10 + f20 + 2*B00*Cx*(xkl + 2*Dx) + Dx*Kx*Px));
	    I[1] += C[0]*W[a]*((f10 + f20 + Dy*Ky*Py + B01*pow(Cy,2) + 2*B00*Cy*(ykl + 2*Dy)));
	    I[2] += C[0]*W[a]*((f10 + f20 + B01*pow(Cz,2) + Dz*Kz*Pz + 2*B00*Cz*(2*Dz + zkl)));
	    I[3] += C[0]*W[a]*(Cy*Dx*(Cz*zkl + Qz));
	    I[4] += C[0]*W[a]*(Cx*Dy*(Cz*zkl + Qz));
	    I[5] += C[0]*W[a]*(Dy*(Pz*zkl + f1));
	    I[6] += C[0]*W[a]*(Dx*(Pz*zkl + f1));
	    I[7] += C[0]*W[a]*(Dz*(f0 + Py*ykl));
	    I[8] += C[0]*W[a]*(Dx*(f0 + Py*ykl));
	    I[9] += C[0]*W[a]*(Dx*Ky*Pz);
	    I[10] += C[0]*W[a]*(Cz*Dx*f17);
	    I[11] += C[0]*W[a]*(Cx*Cz*f22);
	    I[12] += C[0]*W[a]*(Cy*Kx*Qz);
	    I[13] += C[0]*W[a]*(Cy*Dz*f16);
	    I[14] += C[0]*W[a]*(Dz*Ky*Px);
	    I[15] += C[0]*W[a]*(Cx*Dz*f17);
	    I[16] += C[0]*W[a]*(Cx*Ky*Qz);
	    I[17] += C[0]*W[a]*(Cz*Ky*Qx);
	    I[18] += C[0]*W[a]*(Qx*(Cz*zkl + Qz));
	    I[19] += C[0]*W[a]*(Qy*(Cz*zkl + Qz));
	    I[20] += C[0]*W[a]*(Cz*Kx*Qy);
	    I[21] += C[0]*W[a]*(Dz*Kx*Py);
	    I[22] += C[0]*W[a]*(Dx*Kz*Py);
	    I[23] += C[0]*W[a]*(Cy*Kz*Qx);
	    I[24] += C[0]*W[a]*(Dy*Kz*Px);
	    I[25] += C[0]*W[a]*(Cx*Kz*Qy);
	    I[26] += C[0]*W[a]*(Cx*Cy*f15);
	    I[27] += C[0]*W[a]*(Cy*Cz*f7);
	    I[28] += C[0]*W[a]*(Cz*Dy*f16);
	    I[29] += C[0]*W[a]*(Dy*Kx*Pz);
	    I[30] += C[0]*W[a]*(Dy*f9);
	    I[31] += C[0]*W[a]*(Dz*f9);
	    I[32] += C[0]*W[a]*(Qy*f16);
	    I[33] += C[0]*W[a]*(Qz*f16);
	    I[34] += C[0]*W[a]*(Qx*f17);
	    I[35] += C[0]*W[a]*(Qz*f17);
	    I[36] += C[0]*W[a]*(Kx*f0);
	    I[37] += C[0]*W[a]*(Kx*f1);
	    I[38] += C[0]*W[a]*(Ky*f1);
	    I[39] += C[0]*W[a]*(Ky*f13);
	    I[40] += C[0]*W[a]*(Kz*f13);
	    I[41] += C[0]*W[a]*(Kz*f0);
	    I[42] += C[0]*W[a]*(Py*f7);
	    I[43] += C[0]*W[a]*(Pz*f7);
	    I[44] += C[0]*W[a]*(Pz*f22);
	    I[45] += C[0]*W[a]*(Px*f22);
	    I[46] += C[0]*W[a]*(Px*f15);
	    I[47] += C[0]*W[a]*(Py*f15);
	    I[48] += C[0]*W[a]*(Cy*f4);
	    I[49] += C[0]*W[a]*(Cz*f4);
	    I[50] += C[0]*W[a]*(Cz*f2);
	    I[51] += C[0]*W[a]*(Cx*f2);
	    I[52] += C[0]*W[a]*(Cx*f6);
	    I[53] += C[0]*W[a]*(Cy*f6);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[54]) {
	double T[54];
	for (int i = 0; i < 54; ++i) {
	    T[i] = I[i];
	}
	I[0] = T[0];
	I[25] = T[1];
	I[50] = T[2];
	I[41] = T[3];
	I[46] = T[4];
	I[44] = T[5];
	I[38] = T[6];
	I[31] = T[7];
	I[19] = T[8];
	I[20] = T[9];
	I[23] = T[10];
	I[28] = T[11];
	I[17] = T[12];
	I[15] = T[13];
	I[30] = T[14];
	I[33] = T[15];
	I[34] = T[16];
	I[22] = T[17];
	I[40] = T[18];
	I[47] = T[19];
	I[11] = T[20];
	I[13] = T[21];
	I[37] = T[22];
	I[39] = T[23];
	I[42] = T[24];
	I[45] = T[25];
	I[51] = T[26];
	I[5] = T[27];
	I[10] = T[28];
	I[8] = T[29];
	I[6] = T[30];
	I[12] = T[31];
	I[9] = T[32];
	I[16] = T[33];
	I[21] = T[34];
	I[35] = T[35];
	I[7] = T[36];
	I[14] = T[37];
	I[32] = T[38];
	I[18] = T[39];
	I[36] = T[40];
	I[43] = T[41];
	I[1] = T[42];
	I[2] = T[43];
	I[26] = T[44];
	I[24] = T[45];
	I[48] = T[46];
	I[49] = T[47];
	I[3] = T[48];
	I[4] = T[49];
	I[29] = T[50];
	I[27] = T[51];
	I[52] = T[52];
	I[53] = T[53];
    }

};

template<>
struct impl<meta::braket<rysq::SP, rysq::SP, rysq::P, rysq::P> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[4], double (&I)[144]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];

	double xkl = rkl[0];
	double ykl = rkl[1];
	double zkl = rkl[2];

// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Kx = (xkl + Dx);
	    double Ky = (ykl + Dy);
	    double Kz = (Dz + zkl);
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f1 = B00*zkl;
	    double f10 = B00*xkl;
	    double f13 = Cz*pow(Dz,2);
	    double f14 = (B00 + Dz*Iz);
	    double f18 = (B00 + Iy*Ky);
	    double f19 = (B00 + Cy*Ky);
	    double f2 = B01*B10;
	    double f20 = (Dy*Iy + B00);
	    double f21 = Cz*Dz*zkl;
	    double f22 = (Dz*(B10 + Cz*Iz) + B00*(2*Cz + zij));
	    double f23 = B01*Cy;
	    double f24 = (Cx*Ix + B10);
	    double f27 = (Kx*(Cx*Ix + B10) + B00*(xij + 2*Cx));
	    double f28 = 2*B00*Dy;
	    double f3 = (Dz*Kz + B01);
	    double f30 = (Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx));
	    double f31 = (B01 + Dx*Kx);
	    double f32 = (Cy*Iy + B10);
	    double f33 = (B00*(yij + 2*Cy) + Ky*(Cy*Iy + B10));
	    double f34 = (B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10));
	    double f36 = 2*B00*Dz;
	    double f37 = (B01 + Dy*Ky);
	    double f39 = B01*Cx;
	    double f4 = (B00 + Cx*Kx);
	    double f42 = Cx*pow(Dx,2);
	    double f45 = Cy*pow(Dy,2);
	    double f46 = 2*pow(B00,2);
	    double f5 = (B00*(ykl + 2*Dy) + Iy*(B01 + Dy*Ky));
	    double f51 = B01*Cz;
	    double f52 = (B10 + Cz*Iz);
	    double f6 = B00*ykl;
	    double f7 = Cy*Dy*ykl;
	    double f8 = Cx*Dx*xkl;
	    double f9 = 2*B00*Dx;

	    I[0] += C[3]*W[a]*((f46 + f23*yij + 2*Cy*(f28 + f6) + B01*pow(Cy,2) + Dy*Ky*pow(Cy,2) + f2 + B10*pow(Dy,2) + B10*Dy*ykl + yij*(f45 + f28 + f7 + f6)));
	    I[1] += C[3]*W[a]*((f46 + 2*Cx*(f10 + f9) + B10*Dx*xkl + f39*(Cx + xij) + B10*pow(Dx,2) + f2 + xij*(f10 + f42 + f9 + f8) + Dx*Kx*pow(Cx,2)));
	    I[2] += C[3]*W[a]*(Cy*(f10 + f42 + f9 + f8 + B01*Ix + Dx*Kx*xij));
	    I[3] += C[3]*W[a]*(Cz*(f10 + f42 + f9 + f8 + B01*Ix + Dx*Kx*xij));
	    I[4] += C[2]*W[a]*((f10 + f42 + f9 + f8 + B01*Ix + Dx*Kx*xij));
	    I[5] += C[3]*W[a]*(Iy*(f10 + f42 + f39 + f9 + f8));
	    I[6] += C[3]*W[a]*(Iz*(f10 + f42 + f39 + f9 + f8));
	    I[7] += C[1]*W[a]*((f10 + f42 + f39 + f9 + f8));
	    I[8] += C[3]*W[a]*((Dz*Kz*pow(Cz,2) + f46 + zij*(f13 + f21 + f36 + f1) + B10*Dz*zkl + B10*pow(Dz,2) + 2*Cz*(f36 + f1) + f2 + f51*(Cz + zij)));
	    I[9] += C[3]*W[a]*(Cy*(f13 + f21 + f36 + Dz*Kz*zij + f1 + B01*Iz));
	    I[10] += C[3]*W[a]*(Cx*(f13 + f21 + f36 + Dz*Kz*zij + f1 + B01*Iz));
	    I[11] += C[2]*W[a]*((f13 + f21 + f36 + Dz*Kz*zij + f1 + B01*Iz));
	    I[12] += C[3]*W[a]*(Iy*(f13 + f51 + f21 + f36 + f1));
	    I[13] += C[3]*W[a]*(Ix*(f13 + f51 + f21 + f36 + f1));
	    I[14] += C[1]*W[a]*((f13 + f51 + f21 + f36 + f1));
	    I[15] += C[3]*W[a]*(Iz*(f45 + f28 + f23 + f7 + f6));
	    I[16] += C[3]*W[a]*(Ix*(f45 + f28 + f23 + f7 + f6));
	    I[17] += C[1]*W[a]*((f45 + f28 + f23 + f7 + f6));
	    I[18] += C[3]*W[a]*(Dx*(f22 + zkl*(B10 + Cz*Iz)));
	    I[19] += C[3]*W[a]*(Dy*(f22 + zkl*(B10 + Cz*Iz)));
	    I[20] += C[3]*W[a]*(Dy*Ix*(Cz*zkl + Qz));
	    I[21] += C[3]*W[a]*(Qx*(f14 + Iz*zkl));
	    I[22] += C[3]*W[a]*(Cx*Dy*(f14 + Iz*zkl));
	    I[23] += C[2]*W[a]*(Dy*(f14 + Iz*zkl));
	    I[24] += C[3]*W[a]*(Qy*(f14 + Iz*zkl));
	    I[25] += C[2]*W[a]*(Dx*(f14 + Iz*zkl));
	    I[26] += C[3]*W[a]*(Cy*Dx*(f14 + Iz*zkl));
	    I[27] += C[3]*W[a]*(Cy*Dz*(Kx*xij + f4));
	    I[28] += C[2]*W[a]*(Dz*(Kx*xij + f4));
	    I[29] += C[3]*W[a]*(Qy*(Kx*xij + f4));
	    I[30] += C[3]*W[a]*(Cz*Dy*(Kx*xij + f4));
	    I[31] += C[2]*W[a]*(Dy*(Kx*xij + f4));
	    I[32] += C[3]*W[a]*(Qz*(Kx*xij + f4));
	    I[33] += C[3]*W[a]*((Cz*zkl + Qz)*(Dx*xij + Qx));
	    I[34] += C[3]*W[a]*(Dx*Iy*(Cz*zkl + Qz));
	    I[35] += C[1]*W[a]*(Dx*(Cz*zkl + Qz));
	    I[36] += C[1]*W[a]*(Dy*(Cz*zkl + Qz));
	    I[37] += C[3]*W[a]*(f20*(Cz*zkl + Qz));
	    I[38] += C[3]*W[a]*(Cx*Kz*f20);
	    I[39] += C[3]*W[a]*(Cy*Ix*f3);
	    I[40] += C[3]*W[a]*(Ix*Kz*Qy);
	    I[41] += C[3]*W[a]*(Dy*Kz*f24);
	    I[42] += C[2]*W[a]*(Dy*Ix*Kz);
	    I[43] += C[1]*W[a]*(Cx*Dy*Kz);
	    I[44] += C[3]*W[a]*(Dy*Kx*f52);
	    I[45] += C[3]*W[a]*(Dz*Kx*f32);
	    I[46] += C[3]*W[a]*(Dx*Kz*f32);
	    I[47] += C[3]*W[a]*(Cy*Kz*(Dx*xij + Qx));
	    I[48] += C[2]*W[a]*(Kz*(Dx*xij + Qx));
	    I[49] += C[3]*W[a]*(f19*(Dx*xij + Qx));
	    I[50] += C[3]*W[a]*(Cz*Ky*(Dx*xij + Qx));
	    I[51] += C[2]*W[a]*(Ky*(Dx*xij + Qx));
	    I[52] += C[3]*W[a]*(Dx*Ky*f52);
	    I[53] += C[3]*W[a]*(Dx*f33);
	    I[54] += C[3]*W[a]*(Cx*f5);
	    I[55] += C[3]*W[a]*(Cz*f5);
	    I[56] += C[3]*W[a]*(Iz*Ky*Qx);
	    I[57] += C[2]*W[a]*(Dx*Iz*Ky);
	    I[58] += C[3]*W[a]*(Dx*Iz*f19);
	    I[59] += C[1]*W[a]*(Dx*f19);
	    I[60] += C[3]*W[a]*(Dz*Ix*f19);
	    I[61] += C[1]*W[a]*(Dz*f19);
	    I[62] += C[3]*W[a]*(Ky*f30);
	    I[63] += C[3]*W[a]*(Ix*Ky*Qz);
	    I[64] += C[1]*W[a]*(Ky*Qz);
	    I[65] += C[3]*W[a]*(Qx*f18);
	    I[66] += C[3]*W[a]*(Cx*Dz*f18);
	    I[67] += C[2]*W[a]*(Dz*f18);
	    I[68] += C[3]*W[a]*(Dz*f33);
	    I[69] += C[3]*W[a]*(Cy*Kx*f14);
	    I[70] += C[2]*W[a]*(Kx*f14);
	    I[71] += C[3]*W[a]*(Kx*f34);
	    I[72] += C[3]*W[a]*(Iy*Kx*Qz);
	    I[73] += C[1]*W[a]*(Kx*Qz);
	    I[74] += C[3]*W[a]*(Qz*f18);
	    I[75] += C[2]*W[a]*(Dx*f18);
	    I[76] += C[3]*W[a]*(Cz*Dx*f18);
	    I[77] += C[1]*W[a]*(Cz*Dx*Ky);
	    I[78] += C[0]*W[a]*(Dx*Ky);
	    I[79] += C[3]*W[a]*(Dz*Ky*f24);
	    I[80] += C[2]*W[a]*(Dz*Ix*Ky);
	    I[81] += C[0]*W[a]*(Dz*Ky);
	    I[82] += C[1]*W[a]*(Cx*Dz*Ky);
	    I[83] += C[3]*W[a]*(Cx*Ky*f14);
	    I[84] += C[2]*W[a]*(Ky*f14);
	    I[85] += C[3]*W[a]*(f14*f19);
	    I[86] += C[3]*W[a]*(f14*f4);
	    I[87] += C[3]*W[a]*(f20*f4);
	    I[88] += C[1]*W[a]*(Dz*f4);
	    I[89] += C[3]*W[a]*(Dz*Iy*f4);
	    I[90] += C[2]*W[a]*(Dz*Iy*Kx);
	    I[91] += C[1]*W[a]*(Cy*Dz*Kx);
	    I[92] += C[0]*W[a]*(Dz*Kx);
	    I[93] += C[3]*W[a]*(Dz*f27);
	    I[94] += C[3]*W[a]*(Dy*f27);
	    I[95] += C[1]*W[a]*(Dy*f4);
	    I[96] += C[3]*W[a]*(Dy*Iz*f4);
	    I[97] += C[2]*W[a]*(Dy*Iz*Kx);
	    I[98] += C[3]*W[a]*(Iz*Kx*Qy);
	    I[99] += C[1]*W[a]*(Kx*Qy);
	    I[100] += C[2]*W[a]*(Kx*f20);
	    I[101] += C[3]*W[a]*(Cz*Kx*f20);
	    I[102] += C[1]*W[a]*(Cz*Dy*Kx);
	    I[103] += C[0]*W[a]*(Dy*Kx);
	    I[104] += C[3]*W[a]*(Kx*f22);
	    I[105] += C[3]*W[a]*(Ky*f22);
	    I[106] += C[1]*W[a]*(Ky*Qx);
	    I[107] += C[3]*W[a]*(Iy*Kz*Qx);
	    I[108] += C[2]*W[a]*(Dx*Iy*Kz);
	    I[109] += C[1]*W[a]*(Cy*Dx*Kz);
	    I[110] += C[0]*W[a]*(Dx*Kz);
	    I[111] += C[1]*W[a]*(Kz*Qx);
	    I[112] += C[3]*W[a]*(Kz*f30);
	    I[113] += C[0]*W[a]*(Dy*Kz);
	    I[114] += C[1]*W[a]*(Kz*Qy);
	    I[115] += C[2]*W[a]*(Kz*f20);
	    I[116] += C[3]*W[a]*(Kz*f34);
	    I[117] += C[3]*W[a]*(f31*f52);
	    I[118] += C[3]*W[a]*(f37*f52);
	    I[119] += C[3]*W[a]*(f24*f37);
	    I[120] += C[2]*W[a]*(Ix*f37);
	    I[121] += C[3]*W[a]*(Cz*Ix*f37);
	    I[122] += C[1]*W[a]*(Cz*f37);
	    I[123] += C[2]*W[a]*(Iz*f37);
	    I[124] += C[3]*W[a]*(Cx*Iz*f37);
	    I[125] += C[1]*W[a]*(Cx*f37);
	    I[126] += C[3]*W[a]*(Cx*Iy*f3);
	    I[127] += C[1]*W[a]*(Cx*f3);
	    I[128] += C[1]*W[a]*(Cy*f3);
	    I[129] += C[2]*W[a]*(Ix*f3);
	    I[130] += C[3]*W[a]*(f24*f3);
	    I[131] += C[2]*W[a]*(Iy*f3);
	    I[132] += C[3]*W[a]*(f3*f32);
	    I[133] += C[3]*W[a]*(f31*f32);
	    I[134] += C[2]*W[a]*(Iy*f31);
	    I[135] += C[3]*W[a]*(Cz*Iy*f31);
	    I[136] += C[1]*W[a]*(Cz*f31);
	    I[137] += C[2]*W[a]*(Iz*f31);
	    I[138] += C[3]*W[a]*(Cy*Iz*f31);
	    I[139] += C[1]*W[a]*(Cy*f31);
	    I[140] += C[0]*W[a]*(f31);
	    I[141] += C[0]*W[a]*(f37);
	    I[142] += C[2]*W[a]*(f5);
	    I[143] += C[0]*W[a]*(f3);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[144]) {
	double T[144];
	for (int i = 0; i < 144; ++i) {
	    T[i] = I[i];
	}
	I[74] = T[0];
	I[5] = T[1];
	I[6] = T[2];
	I[7] = T[3];
	I[4] = T[4];
	I[9] = T[5];
	I[13] = T[6];
	I[1] = T[7];
	I[143] = T[8];
	I[142] = T[9];
	I[141] = T[10];
	I[140] = T[11];
	I[139] = T[12];
	I[135] = T[13];
	I[131] = T[14];
	I[78] = T[15];
	I[70] = T[16];
	I[66] = T[17];
	I[111] = T[18];
	I[127] = T[19];
	I[119] = T[20];
	I[109] = T[21];
	I[125] = T[22];
	I[124] = T[23];
	I[126] = T[24];
	I[108] = T[25];
	I[110] = T[26];
	I[38] = T[27];
	I[36] = T[28];
	I[22] = T[29];
	I[23] = T[30];
	I[20] = T[31];
	I[39] = T[32];
	I[103] = T[33];
	I[107] = T[34];
	I[99] = T[35];
	I[115] = T[36];
	I[123] = T[37];
	I[121] = T[38];
	I[134] = T[39];
	I[118] = T[40];
	I[117] = T[41];
	I[116] = T[42];
	I[113] = T[43];
	I[31] = T[44];
	I[42] = T[45];
	I[106] = T[46];
	I[102] = T[47];
	I[100] = T[48];
	I[54] = T[49];
	I[55] = T[50];
	I[52] = T[51];
	I[63] = T[52];
	I[58] = T[53];
	I[73] = T[54];
	I[75] = T[55];
	I[61] = T[56];
	I[60] = T[57];
	I[62] = T[58];
	I[50] = T[59];
	I[86] = T[60];
	I[82] = T[61];
	I[53] = T[62];
	I[87] = T[63];
	I[83] = T[64];
	I[57] = T[65];
	I[89] = T[66];
	I[88] = T[67];
	I[90] = T[68];
	I[46] = T[69];
	I[44] = T[70];
	I[26] = T[71];
	I[43] = T[72];
	I[35] = T[73];
	I[91] = T[74];
	I[56] = T[75];
	I[59] = T[76];
	I[51] = T[77];
	I[48] = T[78];
	I[85] = T[79];
	I[84] = T[80];
	I[80] = T[81];
	I[81] = T[82];
	I[93] = T[83];
	I[92] = T[84];
	I[94] = T[85];
	I[45] = T[86];
	I[25] = T[87];
	I[33] = T[88];
	I[41] = T[89];
	I[40] = T[90];
	I[34] = T[91];
	I[32] = T[92];
	I[37] = T[93];
	I[21] = T[94];
	I[17] = T[95];
	I[29] = T[96];
	I[28] = T[97];
	I[30] = T[98];
	I[18] = T[99];
	I[24] = T[100];
	I[27] = T[101];
	I[19] = T[102];
	I[16] = T[103];
	I[47] = T[104];
	I[95] = T[105];
	I[49] = T[106];
	I[105] = T[107];
	I[104] = T[108];
	I[98] = T[109];
	I[96] = T[110];
	I[97] = T[111];
	I[101] = T[112];
	I[112] = T[113];
	I[114] = T[114];
	I[120] = T[115];
	I[122] = T[116];
	I[15] = T[117];
	I[79] = T[118];
	I[69] = T[119];
	I[68] = T[120];
	I[71] = T[121];
	I[67] = T[122];
	I[76] = T[123];
	I[77] = T[124];
	I[65] = T[125];
	I[137] = T[126];
	I[129] = T[127];
	I[130] = T[128];
	I[132] = T[129];
	I[133] = T[130];
	I[136] = T[131];
	I[138] = T[132];
	I[10] = T[133];
	I[8] = T[134];
	I[11] = T[135];
	I[3] = T[136];
	I[12] = T[137];
	I[14] = T[138];
	I[2] = T[139];
	I[0] = T[140];
	I[64] = T[141];
	I[72] = T[142];
	I[128] = T[143];
    }

};

template<>
struct impl<meta::braket<rysq::SP, rysq::P, rysq::SP, rysq::P> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[4], double (&I)[144]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];

	double xkl = rkl[0];
	double ykl = rkl[1];
	double zkl = rkl[2];

// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Kx = (xkl + Dx);
	    double Ky = (ykl + Dy);
	    double Kz = (Dz + zkl);
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f1 = B01*B10;
	    double f14 = (Cx*Ix + B10);
	    double f2 = (B00 + Cx*Kx);
	    double f21 = (B01 + Dx*Kx);
	    double f22 = (Cy*Iy + B10);
	    double f26 = (B01 + Dy*Ky);
	    double f3 = (B00*(ykl + 2*Dy) + Iy*(B01 + Dy*Ky));
	    double f30 = (Dz*Kz + B01);
	    double f32 = 2*pow(B00,2);
	    double f37 = (B10 + Cz*Iz);

	    I[0] += C[3]*W[a]*((f32 + Dy*Ky*(Cy*Iy + B10) + B01*Cy*Iy + B00*(yij + 2*Cy)*(ykl + 2*Dy) + f1));
	    I[1] += C[3]*W[a]*((B00 + Cz*Kz)*(Dx*Ix + B00));
	    I[2] += C[3]*W[a]*((B00 + Cy*Ky)*(Dx*Ix + B00));
	    I[3] += C[3]*W[a]*(Cz*Dx*(B00 + Iy*Ky));
	    I[4] += C[3]*W[a]*((B00 + Cz*Kz)*(Dy*Iy + B00));
	    I[5] += C[3]*W[a]*((Dz*Kz*(B10 + Cz*Iz) + f32 + B00*(2*Cz + zij)*(2*Dz + zkl) + f1 + B01*Cz*Iz));
	    I[6] += C[3]*W[a]*(Cx*(B00*(2*Dz + zkl) + Iz*(Dz*Kz + B01)));
	    I[7] += C[3]*W[a]*(Cy*(B00*(2*Dz + zkl) + Iz*(Dz*Kz + B01)));
	    I[8] += C[2]*W[a]*((B00*(2*Dz + zkl) + Iz*(Dz*Kz + B01)));
	    I[9] += C[3]*W[a]*(Dx*Iy*(B00 + Cz*Kz));
	    I[10] += C[3]*W[a]*(Iy*Kx*Qz);
	    I[11] += C[3]*W[a]*(Iy*(B01*Cz + B00*(2*Dz + zkl) + Cz*Dz*Kz));
	    I[12] += C[3]*W[a]*(Ix*(B01*Cz + B00*(2*Dz + zkl) + Cz*Dz*Kz));
	    I[13] += C[3]*W[a]*(Dy*Ix*(B00 + Cz*Kz));
	    I[14] += C[1]*W[a]*(Ix*(B00 + Cz*Kz));
	    I[15] += C[1]*W[a]*(Iy*(B00 + Cz*Kz));
	    I[16] += C[3]*W[a]*(f2*(B00 + Dz*Iz));
	    I[17] += C[3]*W[a]*(Qz*(Kx*xij + f2));
	    I[18] += C[3]*W[a]*(Dz*(Kx*(Cx*Ix + B10) + B00*(xij + 2*Cx)));
	    I[19] += C[3]*W[a]*(Dy*(Kx*(Cx*Ix + B10) + B00*(xij + 2*Cx)));
	    I[20] += C[3]*W[a]*((B00*(xij + 2*Cx)*(xkl + 2*Dx) + f32 + B01*Cx*Ix + f1 + Dx*Kx*(Cx*Ix + B10)));
	    I[21] += C[3]*W[a]*(Iz*(B01*Cx + Cx*Dx*Kx + B00*(xkl + 2*Dx)));
	    I[22] += C[3]*W[a]*(Iz*Ky*Qx);
	    I[23] += C[3]*W[a]*(Cz*Ky*(Dx*Ix + B00));
	    I[24] += C[3]*W[a]*(Ky*(Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx)));
	    I[25] += C[2]*W[a]*(Ky*(Dx*Ix + B00));
	    I[26] += C[3]*W[a]*(Iy*Kz*Qx);
	    I[27] += C[3]*W[a]*(Iy*(B01*Cx + Cx*Dx*Kx + B00*(xkl + 2*Dx)));
	    I[28] += C[3]*W[a]*(Cz*(Ix*(B01 + Dx*Kx) + B00*(xkl + 2*Dx)));
	    I[29] += C[2]*W[a]*((Ix*(B01 + Dx*Kx) + B00*(xkl + 2*Dx)));
	    I[30] += C[3]*W[a]*(Cy*(Ix*(B01 + Dx*Kx) + B00*(xkl + 2*Dx)));
	    I[31] += C[3]*W[a]*(Qy*(Kx*xij + f2));
	    I[32] += C[2]*W[a]*(Dy*(Kx*xij + f2));
	    I[33] += C[3]*W[a]*(Cz*Dy*(Kx*xij + f2));
	    I[34] += C[1]*W[a]*(Cz*(Kx*xij + f2));
	    I[35] += C[2]*W[a]*(Dz*(Kx*xij + f2));
	    I[36] += C[3]*W[a]*(Cy*Dz*(Kx*xij + f2));
	    I[37] += C[1]*W[a]*(Cy*(Kx*xij + f2));
	    I[38] += C[0]*W[a]*((Kx*xij + f2));
	    I[39] += C[1]*W[a]*((Kx*(Cx*Ix + B10) + B00*(xij + 2*Cx)));
	    I[40] += C[3]*W[a]*(Kz*(Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx)));
	    I[41] += C[2]*W[a]*(Kz*(Dx*Ix + B00));
	    I[42] += C[3]*W[a]*(Cy*Kz*(Dx*Ix + B00));
	    I[43] += C[3]*W[a]*(Ix*Kz*Qy);
	    I[44] += C[3]*W[a]*(Kz*(B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10)));
	    I[45] += C[3]*W[a]*(Cx*Kz*(Dy*Iy + B00));
	    I[46] += C[2]*W[a]*(Kz*(Dy*Iy + B00));
	    I[47] += C[3]*W[a]*(f2*(Dy*Iy + B00));
	    I[48] += C[3]*W[a]*(Cz*Kx*(Dy*Iy + B00));
	    I[49] += C[3]*W[a]*(Kx*(B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10)));
	    I[50] += C[2]*W[a]*(Kx*(Dy*Iy + B00));
	    I[51] += C[3]*W[a]*(Iz*Kx*Qy);
	    I[52] += C[3]*W[a]*(Dx*Iz*(B00 + Cy*Ky));
	    I[53] += C[1]*W[a]*(Iz*(B00 + Cy*Ky));
	    I[54] += C[3]*W[a]*(Iz*(Cy*Dy*Ky + B01*Cy + B00*(ykl + 2*Dy)));
	    I[55] += C[3]*W[a]*(Ix*(Cy*Dy*Ky + B01*Cy + B00*(ykl + 2*Dy)));
	    I[56] += C[1]*W[a]*(Ix*(B00 + Cy*Ky));
	    I[57] += C[3]*W[a]*(Dz*Ix*(B00 + Cy*Ky));
	    I[58] += C[3]*W[a]*(Cx*Dz*(B00 + Iy*Ky));
	    I[59] += C[2]*W[a]*(Dz*(B00 + Iy*Ky));
	    I[60] += C[3]*W[a]*(Dz*(B00*(yij + 2*Cy) + Ky*(Cy*Iy + B10)));
	    I[61] += C[1]*W[a]*((B00*(yij + 2*Cy) + Ky*(Cy*Iy + B10)));
	    I[62] += C[3]*W[a]*(Dx*(B00*(yij + 2*Cy) + Ky*(Cy*Iy + B10)));
	    I[63] += C[2]*W[a]*(Dx*(B00 + Iy*Ky));
	    I[64] += C[3]*W[a]*(Qx*(B00 + Iy*Ky));
	    I[65] += C[1]*W[a]*(Cx*(B00 + Iy*Ky));
	    I[66] += C[0]*W[a]*((B00 + Iy*Ky));
	    I[67] += C[1]*W[a]*(Cz*(B00 + Iy*Ky));
	    I[68] += C[3]*W[a]*(Qz*(B00 + Iy*Ky));
	    I[69] += C[3]*W[a]*((B00 + Cy*Ky)*(B00 + Dz*Iz));
	    I[70] += C[3]*W[a]*(Qy*(Iz*Kz + B00));
	    I[71] += C[3]*W[a]*(Cx*Dy*(Iz*Kz + B00));
	    I[72] += C[2]*W[a]*(Dy*(Iz*Kz + B00));
	    I[73] += C[3]*W[a]*(Dy*(Kz*(B10 + Cz*Iz) + B00*(2*Cz + zij)));
	    I[74] += C[1]*W[a]*((Kz*(B10 + Cz*Iz) + B00*(2*Cz + zij)));
	    I[75] += C[3]*W[a]*(Dx*(Kz*(B10 + Cz*Iz) + B00*(2*Cz + zij)));
	    I[76] += C[3]*W[a]*(Cy*Dx*(Iz*Kz + B00));
	    I[77] += C[2]*W[a]*(Dx*(Iz*Kz + B00));
	    I[78] += C[3]*W[a]*(Qx*(Iz*Kz + B00));
	    I[79] += C[1]*W[a]*(Cx*(Iz*Kz + B00));
	    I[80] += C[0]*W[a]*((Iz*Kz + B00));
	    I[81] += C[1]*W[a]*(Cy*(Iz*Kz + B00));
	    I[82] += C[3]*W[a]*(Cy*Kx*(B00 + Dz*Iz));
	    I[83] += C[2]*W[a]*(Kx*(B00 + Dz*Iz));
	    I[84] += C[3]*W[a]*(Kx*(Dz*(B10 + Cz*Iz) + B00*(2*Cz + zij)));
	    I[85] += C[3]*W[a]*(Ky*(Dz*(B10 + Cz*Iz) + B00*(2*Cz + zij)));
	    I[86] += C[3]*W[a]*(Cx*Ky*(B00 + Dz*Iz));
	    I[87] += C[2]*W[a]*(Ky*(B00 + Dz*Iz));
	    I[88] += C[3]*W[a]*(Ix*Ky*Qz);
	    I[89] += C[2]*W[a]*(Dz*Ix*Ky);
	    I[90] += C[3]*W[a]*(Dz*Ky*f14);
	    I[91] += C[3]*W[a]*(Dx*Kz*f22);
	    I[92] += C[2]*W[a]*(Dx*Iy*Kz);
	    I[93] += C[1]*W[a]*(Cx*Iy*Kz);
	    I[94] += C[3]*W[a]*(Cx*Iy*f30);
	    I[95] += C[3]*W[a]*(Cx*Iz*f26);
	    I[96] += C[3]*W[a]*(Cy*Iz*f21);
	    I[97] += C[1]*W[a]*(Cy*Iz*Kx);
	    I[98] += C[3]*W[a]*(Dz*Kx*f22);
	    I[99] += C[1]*W[a]*(Kx*f22);
	    I[100] += C[2]*W[a]*(Dz*Iy*Kx);
	    I[101] += C[3]*W[a]*(Dz*Iy*f2);
	    I[102] += C[1]*W[a]*(Iy*f2);
	    I[103] += C[3]*W[a]*(Dy*Iz*f2);
	    I[104] += C[1]*W[a]*(Iz*f2);
	    I[105] += C[0]*W[a]*(Iz*Kx);
	    I[106] += C[2]*W[a]*(Dy*Iz*Kx);
	    I[107] += C[3]*W[a]*(Dy*Kx*f37);
	    I[108] += C[1]*W[a]*(Kx*f37);
	    I[109] += C[0]*W[a]*(Iy*Kx);
	    I[110] += C[1]*W[a]*(Cz*Iy*Kx);
	    I[111] += C[3]*W[a]*(Cz*Iy*f21);
	    I[112] += C[2]*W[a]*(Iy*f21);
	    I[113] += C[3]*W[a]*(f21*f22);
	    I[114] += C[2]*W[a]*(Iz*f21);
	    I[115] += C[3]*W[a]*(f21*f37);
	    I[116] += C[3]*W[a]*(Dx*Ky*f37);
	    I[117] += C[2]*W[a]*(Dx*Iz*Ky);
	    I[118] += C[1]*W[a]*(Cx*Iz*Ky);
	    I[119] += C[0]*W[a]*(Iz*Ky);
	    I[120] += C[1]*W[a]*(Ky*f37);
	    I[121] += C[1]*W[a]*(Ky*f14);
	    I[122] += C[3]*W[a]*(f14*f26);
	    I[123] += C[2]*W[a]*(Iz*f26);
	    I[124] += C[3]*W[a]*(f26*f37);
	    I[125] += C[2]*W[a]*(Ix*f26);
	    I[126] += C[3]*W[a]*(Cz*Ix*f26);
	    I[127] += C[1]*W[a]*(Cz*Ix*Ky);
	    I[128] += C[0]*W[a]*(Ix*Ky);
	    I[129] += C[2]*W[a]*(Dy*Ix*Kz);
	    I[130] += C[3]*W[a]*(Dy*Kz*f14);
	    I[131] += C[1]*W[a]*(Kz*f14);
	    I[132] += C[0]*W[a]*(Iy*Kz);
	    I[133] += C[1]*W[a]*(Kz*f22);
	    I[134] += C[0]*W[a]*(Ix*Kz);
	    I[135] += C[1]*W[a]*(Cy*Ix*Kz);
	    I[136] += C[3]*W[a]*(Cy*Ix*f30);
	    I[137] += C[2]*W[a]*(Ix*f30);
	    I[138] += C[3]*W[a]*(f14*f30);
	    I[139] += C[2]*W[a]*(Iy*f30);
	    I[140] += C[3]*W[a]*(f22*f30);
	    I[141] += C[3]*W[a]*(Cz*f3);
	    I[142] += C[3]*W[a]*(Cx*f3);
	    I[143] += C[2]*W[a]*(f3);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[144]) {
	double T[144];
	for (int i = 0; i < 144; ++i) {
	    T[i] = I[i];
	}
	I[78] = T[0];
	I[111] = T[1];
	I[62] = T[2];
	I[67] = T[3];
	I[127] = T[4];
	I[143] = T[5];
	I[141] = T[6];
	I[142] = T[7];
	I[140] = T[8];
	I[115] = T[9];
	I[43] = T[10];
	I[139] = T[11];
	I[135] = T[12];
	I[123] = T[13];
	I[99] = T[14];
	I[103] = T[15];
	I[45] = T[16];
	I[39] = T[17];
	I[37] = T[18];
	I[25] = T[19];
	I[13] = T[20];
	I[21] = T[21];
	I[69] = T[22];
	I[63] = T[23];
	I[61] = T[24];
	I[60] = T[25];
	I[113] = T[26];
	I[17] = T[27];
	I[15] = T[28];
	I[12] = T[29];
	I[14] = T[30];
	I[26] = T[31];
	I[24] = T[32];
	I[27] = T[33];
	I[3] = T[34];
	I[36] = T[35];
	I[38] = T[36];
	I[2] = T[37];
	I[0] = T[38];
	I[1] = T[39];
	I[109] = T[40];
	I[108] = T[41];
	I[110] = T[42];
	I[122] = T[43];
	I[126] = T[44];
	I[125] = T[45];
	I[124] = T[46];
	I[29] = T[47];
	I[31] = T[48];
	I[30] = T[49];
	I[28] = T[50];
	I[34] = T[51];
	I[70] = T[52];
	I[58] = T[53];
	I[82] = T[54];
	I[74] = T[55];
	I[50] = T[56];
	I[86] = T[57];
	I[89] = T[58];
	I[88] = T[59];
	I[90] = T[60];
	I[54] = T[61];
	I[66] = T[62];
	I[64] = T[63];
	I[65] = T[64];
	I[53] = T[65];
	I[52] = T[66];
	I[55] = T[67];
	I[91] = T[68];
	I[94] = T[69];
	I[130] = T[70];
	I[129] = T[71];
	I[128] = T[72];
	I[131] = T[73];
	I[107] = T[74];
	I[119] = T[75];
	I[118] = T[76];
	I[116] = T[77];
	I[117] = T[78];
	I[105] = T[79];
	I[104] = T[80];
	I[106] = T[81];
	I[46] = T[82];
	I[44] = T[83];
	I[47] = T[84];
	I[95] = T[85];
	I[93] = T[86];
	I[92] = T[87];
	I[87] = T[88];
	I[84] = T[89];
	I[85] = T[90];
	I[114] = T[91];
	I[112] = T[92];
	I[101] = T[93];
	I[137] = T[94];
	I[81] = T[95];
	I[22] = T[96];
	I[10] = T[97];
	I[42] = T[98];
	I[6] = T[99];
	I[40] = T[100];
	I[41] = T[101];
	I[5] = T[102];
	I[33] = T[103];
	I[9] = T[104];
	I[8] = T[105];
	I[32] = T[106];
	I[35] = T[107];
	I[11] = T[108];
	I[4] = T[109];
	I[7] = T[110];
	I[19] = T[111];
	I[16] = T[112];
	I[18] = T[113];
	I[20] = T[114];
	I[23] = T[115];
	I[71] = T[116];
	I[68] = T[117];
	I[57] = T[118];
	I[56] = T[119];
	I[59] = T[120];
	I[49] = T[121];
	I[73] = T[122];
	I[80] = T[123];
	I[83] = T[124];
	I[72] = T[125];
	I[75] = T[126];
	I[51] = T[127];
	I[48] = T[128];
	I[120] = T[129];
	I[121] = T[130];
	I[97] = T[131];
	I[100] = T[132];
	I[102] = T[133];
	I[96] = T[134];
	I[98] = T[135];
	I[134] = T[136];
	I[132] = T[137];
	I[133] = T[138];
	I[136] = T[139];
	I[138] = T[140];
	I[79] = T[141];
	I[77] = T[142];
	I[76] = T[143];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::P, rysq::P, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[1], double (&I)[90]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f0 = (3*B10 + pow(Cx,2));
	    double f10 = (3*B00*Pz + Cz*Dz*(3*B10 + pow(Cz,2)));
	    double f11 = (3*pow(B10,2) + 3*B10*Cz*(2*Cz + zij) + Iz*pow(Cz,3));
	    double f12 = (B00 + Dz*Iz);
	    double f14 = (Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx));
	    double f16 = (Dy*Iy + B00);
	    double f17 = (Cy*Iy + B10);
	    double f18 = (3*pow(B10,2) + 3*B10*Cx*(xij + 2*Cx) + Ix*pow(Cx,3));
	    double f2 = (Dz*Pz + 2*B00*Cz);
	    double f20 = (Dx*Px + 2*B00*Cx);
	    double f22 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f28 = (Cx*Dx*(3*B10 + pow(Cx,2)) + 3*B00*Px);
	    double f3 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f31 = (3*B10 + pow(Cz,2));
	    double f32 = (3*B10 + pow(Cy,2));
	    double f35 = (3*pow(B10,2) + Iy*pow(Cy,3) + 3*B10*Cy*(yij + 2*Cy));
	    double f4 = (Dz*(Iz*pow(Cz,2) + B10*(3*Cz + zij)) + 3*B00*Pz + 2*B00*Cz*zij);
	    double f5 = (B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10));
	    double f6 = (Dx*(B10*(3*Cx + xij) + Ix*pow(Cx,2)) + 2*B00*Cx*xij + 3*B00*Px);
	    double f7 = (Cy*Dy*(3*B10 + pow(Cy,2)) + 3*B00*Py);
	    double f8 = (Dy*Py + 2*B00*Cy);
	    double f9 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));

	    I[0] += C[0]*W[a]*(Cx*(f7 + yij*(Dy*Py + 2*B00*Cy)));
	    I[1] += C[0]*W[a]*(Cz*(f7 + yij*(Dy*Py + 2*B00*Cy)));
	    I[2] += C[0]*W[a]*((3*B00*Py*yij + Dy*(3*pow(B10,2) + Iy*pow(Cy,3) + 3*B10*Cy*(yij + 2*Cy)) + 4*B00*Cy*f32));
	    I[3] += C[0]*W[a]*(Dz*Py*(Px + Cx*xij));
	    I[4] += C[0]*W[a]*(Cy*Qz*(Px + Cx*xij));
	    I[5] += C[0]*W[a]*(Cz*Qy*(Px + Cx*xij));
	    I[6] += C[0]*W[a]*(Cx*Qy*(Cz*zij + Pz));
	    I[7] += C[0]*W[a]*(Cx*Cy*(f2 + Qz*zij));
	    I[8] += C[0]*W[a]*(Py*(f2 + Qz*zij));
	    I[9] += C[0]*W[a]*(Px*(f2 + Qz*zij));
	    I[10] += C[0]*W[a]*((3*B00*Pz*zij + Dz*(3*pow(B10,2) + 3*B10*Cz*(2*Cz + zij) + Iz*pow(Cz,3)) + 4*B00*Cz*f31));
	    I[11] += C[0]*W[a]*((3*B00*Px*xij + 4*B00*Cx*f0 + Dx*(3*pow(B10,2) + 3*B10*Cx*(xij + 2*Cx) + Ix*pow(Cx,3))));
	    I[12] += C[0]*W[a]*(Cx*Dy*Iz*f0);
	    I[13] += C[0]*W[a]*(Dy*Pz*(Px + Cx*xij));
	    I[14] += C[0]*W[a]*(Dy*Px*(Cz*zij + Pz));
	    I[15] += C[0]*W[a]*(Cy*Qx*(Cz*zij + Pz));
	    I[16] += C[0]*W[a]*(Cy*Pz*(Dx*xij + Qx));
	    I[17] += C[0]*W[a]*(Cy*f32*(Dx*xij + Qx));
	    I[18] += C[0]*W[a]*(Cz*f31*(Dx*xij + Qx));
	    I[19] += C[0]*W[a]*(Cz*Py*(Dx*xij + Qx));
	    I[20] += C[0]*W[a]*(Dx*Py*(Cz*zij + Pz));
	    I[21] += C[0]*W[a]*(f20*(Cz*zij + Pz));
	    I[22] += C[0]*W[a]*(f8*(Cz*zij + Pz));
	    I[23] += C[0]*W[a]*(Cx*Cz*f5);
	    I[24] += C[0]*W[a]*(Cy*Dz*Ix*f32);
	    I[25] += C[0]*W[a]*(Cy*Dz*f22);
	    I[26] += C[0]*W[a]*(Ix*Py*Qz);
	    I[27] += C[0]*W[a]*(Cy*Ix*f2);
	    I[28] += C[0]*W[a]*(f2*(Px + Cx*xij));
	    I[29] += C[0]*W[a]*(f8*(Px + Cx*xij));
	    I[30] += C[0]*W[a]*(Cx*Iz*f8);
	    I[31] += C[0]*W[a]*(Iz*Px*Qy);
	    I[32] += C[0]*W[a]*(Dz*Px*f17);
	    I[33] += C[0]*W[a]*(Iy*Px*Qz);
	    I[34] += C[0]*W[a]*(Cx*Dz*Iy*f0);
	    I[35] += C[0]*W[a]*(Cx*Dz*f9);
	    I[36] += C[0]*W[a]*(Cx*Iy*f2);
	    I[37] += C[0]*W[a]*(Cx*Qz*f17);
	    I[38] += C[0]*W[a]*(Cx*Dy*f3);
	    I[39] += C[0]*W[a]*(Cy*Dx*f3);
	    I[40] += C[0]*W[a]*(Cy*Iz*f20);
	    I[41] += C[0]*W[a]*(Cy*Dx*Iz*f32);
	    I[42] += C[0]*W[a]*(Cy*f12*f32);
	    I[43] += C[0]*W[a]*(Cy*Px*f12);
	    I[44] += C[0]*W[a]*(Cz*Px*f16);
	    I[45] += C[0]*W[a]*(Cz*f16*f31);
	    I[46] += C[0]*W[a]*(Cz*Dy*f22);
	    I[47] += C[0]*W[a]*(Cz*Dy*Ix*f31);
	    I[48] += C[0]*W[a]*(Cz*Ix*f8);
	    I[49] += C[0]*W[a]*(Ix*Pz*Qy);
	    I[50] += C[0]*W[a]*(Cx*Pz*f16);
	    I[51] += C[0]*W[a]*(Cx*f0*f16);
	    I[52] += C[0]*W[a]*(Cx*f0*f12);
	    I[53] += C[0]*W[a]*(Cx*Py*f12);
	    I[54] += C[0]*W[a]*(Iz*Py*Qx);
	    I[55] += C[0]*W[a]*(Cz*Qx*f17);
	    I[56] += C[0]*W[a]*(Iy*Pz*Qx);
	    I[57] += C[0]*W[a]*(Dx*Pz*f17);
	    I[58] += C[0]*W[a]*(Cz*Dx*f9);
	    I[59] += C[0]*W[a]*(Cz*Dx*Iy*f31);
	    I[60] += C[0]*W[a]*(Cz*Iy*f20);
	    I[61] += C[0]*W[a]*(Cy*Cz*f14);
	    I[62] += C[0]*W[a]*(Cy*f6);
	    I[63] += C[0]*W[a]*(Cz*f6);
	    I[64] += C[0]*W[a]*(Py*f14);
	    I[65] += C[0]*W[a]*(Pz*f14);
	    I[66] += C[0]*W[a]*(Px*f5);
	    I[67] += C[0]*W[a]*(Pz*f5);
	    I[68] += C[0]*W[a]*(Dx*f35);
	    I[69] += C[0]*W[a]*(Dx*f11);
	    I[70] += C[0]*W[a]*(Dy*f11);
	    I[71] += C[0]*W[a]*(Dy*f18);
	    I[72] += C[0]*W[a]*(Dz*f18);
	    I[73] += C[0]*W[a]*(Dz*f35);
	    I[74] += C[0]*W[a]*(Ix*f10);
	    I[75] += C[0]*W[a]*(Ix*f7);
	    I[76] += C[0]*W[a]*(Iz*f7);
	    I[77] += C[0]*W[a]*(Iz*f28);
	    I[78] += C[0]*W[a]*(Iy*f28);
	    I[79] += C[0]*W[a]*(Iy*f10);
	    I[80] += C[0]*W[a]*(Qz*f22);
	    I[81] += C[0]*W[a]*(Qy*f22);
	    I[82] += C[0]*W[a]*(Qy*f3);
	    I[83] += C[0]*W[a]*(Qx*f3);
	    I[84] += C[0]*W[a]*(Qx*f9);
	    I[85] += C[0]*W[a]*(Qz*f9);
	    I[86] += C[0]*W[a]*(f17*f20);
	    I[87] += C[0]*W[a]*(f17*f2);
	    I[88] += C[0]*W[a]*(Cx*f4);
	    I[89] += C[0]*W[a]*(Cy*f4);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[90]) {
	double T[90];
	for (int i = 0; i < 90; ++i) {
	    T[i] = I[i];
	}
	I[45] = T[0];
	I[46] = T[1];
	I[41] = T[2];
	I[65] = T[3];
	I[69] = T[4];
	I[39] = T[5];
	I[59] = T[6];
	I[89] = T[7];
	I[86] = T[8];
	I[84] = T[9];
	I[82] = T[10];
	I[0] = T[11];
	I[50] = T[12];
	I[37] = T[13];
	I[54] = T[14];
	I[29] = T[15];
	I[8] = T[16];
	I[1] = T[17];
	I[2] = T[18];
	I[6] = T[19];
	I[26] = T[20];
	I[24] = T[21];
	I[56] = T[22];
	I[49] = T[23];
	I[61] = T[24];
	I[63] = T[25];
	I[66] = T[26];
	I[68] = T[27];
	I[67] = T[28];
	I[35] = T[29];
	I[55] = T[30];
	I[53] = T[31];
	I[73] = T[32];
	I[74] = T[33];
	I[70] = T[34];
	I[75] = T[35];
	I[77] = T[36];
	I[79] = T[37];
	I[57] = T[38];
	I[28] = T[39];
	I[23] = T[40];
	I[21] = T[41];
	I[81] = T[42];
	I[83] = T[43];
	I[44] = T[44];
	I[42] = T[45];
	I[34] = T[46];
	I[32] = T[47];
	I[36] = T[48];
	I[38] = T[49];
	I[47] = T[50];
	I[40] = T[51];
	I[80] = T[52];
	I[85] = T[53];
	I[25] = T[54];
	I[19] = T[55];
	I[17] = T[56];
	I[18] = T[57];
	I[16] = T[58];
	I[12] = T[59];
	I[14] = T[60];
	I[9] = T[61];
	I[3] = T[62];
	I[4] = T[63];
	I[5] = T[64];
	I[7] = T[65];
	I[43] = T[66];
	I[48] = T[67];
	I[11] = T[68];
	I[22] = T[69];
	I[52] = T[70];
	I[30] = T[71];
	I[60] = T[72];
	I[71] = T[73];
	I[62] = T[74];
	I[31] = T[75];
	I[51] = T[76];
	I[20] = T[77];
	I[10] = T[78];
	I[72] = T[79];
	I[64] = T[80];
	I[33] = T[81];
	I[58] = T[82];
	I[27] = T[83];
	I[15] = T[84];
	I[76] = T[85];
	I[13] = T[86];
	I[78] = T[87];
	I[87] = T[88];
	I[88] = T[89];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::D, rysq::S, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[1], double (&I)[60]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 

// 	for (int a = 0; a < N; ++a) {
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);


	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double f0 = (3*B10 + pow(Cx,2));
	    double f10 = (Cy*Iy + B10);
	    double f11 = (Cy*pow(Iy,2) + B10*(3*Cy + 2*yij));
	    double f12 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f13 = (B10 + Cz*Iz);
	    double f14 = (3*pow(B10,2) + 3*B10*Cx*(xij + 2*Cx) + Ix*pow(Cx,3));
	    double f15 = (B10 + pow(Iy,2));
	    double f16 = (B10 + pow(Iz,2));
	    double f19 = (3*pow(B10,2) + 3*B10*Cz*(2*Cz + zij) + Iz*pow(Cz,3));
	    double f2 = (B10*(3*Cx + 2*xij) + Cx*pow(Ix,2));
	    double f21 = (3*B10 + pow(Cz,2));
	    double f22 = (3*B10 + pow(Cy,2));
	    double f23 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f26 = (3*pow(B10,2) + Iy*pow(Cy,3) + 3*B10*Cy*(yij + 2*Cy));
	    double f3 = (3*pow(B10,2) + pow(Cz,2)*pow(Iz,2) + B10*(6*Cz*zij + pow(zij,2) + 6*pow(Cz,2)));
	    double f4 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));
	    double f5 = (3*pow(B10,2) + B10*(6*Cx*xij + pow(xij,2) + 6*pow(Cx,2)) + pow(Cx,2)*pow(Ix,2));
	    double f7 = (B10*(3*Cz + 2*zij) + Cz*pow(Iz,2));
	    double f9 = (3*pow(B10,2) + B10*(6*Cy*yij + pow(yij,2) + 6*pow(Cy,2)) + pow(Cy,2)*pow(Iy,2));

	    I[0] += C[0]*W[a]*(Cz*f21*(xij*(xij + 2*Cx) + Px));
	    I[1] += C[0]*W[a]*(Cz*f10*(Px + Cx*xij));
	    I[2] += C[0]*W[a]*(Cz*Py*(xij*(xij + 2*Cx) + Px));
	    I[3] += C[0]*W[a]*(Iz*Py*(Px + Cx*xij));
	    I[4] += C[0]*W[a]*(Cy*f22*(xij*(xij + 2*Cx) + Px));
	    I[5] += C[0]*W[a]*(Cy*f13*(Px + Cx*xij));
	    I[6] += C[0]*W[a]*(Cy*Pz*(xij*(xij + 2*Cx) + Px));
	    I[7] += C[0]*W[a]*(Iy*Pz*(Px + Cx*xij));
	    I[8] += C[0]*W[a]*(Cx*Iy*Iz*f0);
	    I[9] += C[0]*W[a]*(Cx*Cy*f7);
	    I[10] += C[0]*W[a]*(Ix*Pz*f10);
	    I[11] += C[0]*W[a]*(Cy*Ix*Iz*f22);
	    I[12] += C[0]*W[a]*(Cy*Iz*f12);
	    I[13] += C[0]*W[a]*((3*pow(B10,2)*(5*Cz + 2*zij) + B10*Cz*(3*pow(zij,2) + 10*pow(Cz,2) + 12*Cz*zij) + pow(Cz,3)*pow(Iz,2)));
	    I[14] += C[0]*W[a]*(Cz*Ix*f4);
	    I[15] += C[0]*W[a]*(Ix*Py*f13);
	    I[16] += C[0]*W[a]*(Cy*Ix*f23);
	    I[17] += C[0]*W[a]*(Cy*Px*f16);
	    I[18] += C[0]*W[a]*(Cy*f16*f22);
	    I[19] += C[0]*W[a]*(Cx*Py*f16);
	    I[20] += C[0]*W[a]*(Cx*f0*f16);
	    I[21] += C[0]*W[a]*(Cx*Pz*f15);
	    I[22] += C[0]*W[a]*(Cx*f0*f15);
	    I[23] += C[0]*W[a]*(Cx*Cz*f11);
	    I[24] += C[0]*W[a]*(Cz*Iy*f12);
	    I[25] += C[0]*W[a]*(Cz*Ix*Iy*f21);
	    I[26] += C[0]*W[a]*(Cz*f15*f21);
	    I[27] += C[0]*W[a]*(Cz*Px*f15);
	    I[28] += C[0]*W[a]*(Iz*Px*f10);
	    I[29] += C[0]*W[a]*(Iy*Px*f13);
	    I[30] += C[0]*W[a]*(f23*(Px + Cx*xij));
	    I[31] += C[0]*W[a]*(f4*(Px + Cx*xij));
	    I[32] += C[0]*W[a]*(Cx*Iz*f4);
	    I[33] += C[0]*W[a]*(Cx*Iy*f23);
	    I[34] += C[0]*W[a]*(Cx*f10*f13);
	    I[35] += C[0]*W[a]*((B10*Cx*(12*Cx*xij + 10*pow(Cx,2) + 3*pow(xij,2)) + 3*pow(B10,2)*(5*Cx + 2*xij) + pow(Cx,3)*pow(Ix,2)));
	    I[36] += C[0]*W[a]*((3*pow(B10,2)*(5*Cy + 2*yij) + B10*Cy*(12*Cy*yij + 3*pow(yij,2) + 10*pow(Cy,2)) + pow(Cy,3)*pow(Iy,2)));
	    I[37] += C[0]*W[a]*(Cy*Cz*f2);
	    I[38] += C[0]*W[a]*(Py*f2);
	    I[39] += C[0]*W[a]*(Pz*f2);
	    I[40] += C[0]*W[a]*(Pz*f11);
	    I[41] += C[0]*W[a]*(Px*f11);
	    I[42] += C[0]*W[a]*(Px*f7);
	    I[43] += C[0]*W[a]*(Py*f7);
	    I[44] += C[0]*W[a]*(Cy*f5);
	    I[45] += C[0]*W[a]*(Cz*f5);
	    I[46] += C[0]*W[a]*(Cz*f9);
	    I[47] += C[0]*W[a]*(Cx*f9);
	    I[48] += C[0]*W[a]*(Cx*f3);
	    I[49] += C[0]*W[a]*(Cy*f3);
	    I[50] += C[0]*W[a]*(Ix*f19);
	    I[51] += C[0]*W[a]*(Ix*f26);
	    I[52] += C[0]*W[a]*(Iz*f26);
	    I[53] += C[0]*W[a]*(Iz*f14);
	    I[54] += C[0]*W[a]*(Iy*f14);
	    I[55] += C[0]*W[a]*(Iy*f19);
	    I[56] += C[0]*W[a]*(f13*f4);
	    I[57] += C[0]*W[a]*(f12*f13);
	    I[58] += C[0]*W[a]*(f10*f12);
	    I[59] += C[0]*W[a]*(f10*f23);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[60]) {
	double T[60];
	for (int i = 0; i < 60; ++i) {
	    T[i] = I[i];
	}
	I[2] = T[0];
	I[39] = T[1];
	I[6] = T[2];
	I[45] = T[3];
	I[1] = T[4];
	I[49] = T[5];
	I[8] = T[6];
	I[37] = T[7];
	I[50] = T[8];
	I[29] = T[9];
	I[38] = T[10];
	I[41] = T[11];
	I[43] = T[12];
	I[22] = T[13];
	I[36] = T[14];
	I[46] = T[15];
	I[48] = T[16];
	I[23] = T[17];
	I[21] = T[18];
	I[25] = T[19];
	I[20] = T[20];
	I[17] = T[21];
	I[10] = T[22];
	I[19] = T[23];
	I[34] = T[24];
	I[32] = T[25];
	I[12] = T[26];
	I[14] = T[27];
	I[53] = T[28];
	I[54] = T[29];
	I[47] = T[30];
	I[35] = T[31];
	I[55] = T[32];
	I[57] = T[33];
	I[59] = T[34];
	I[0] = T[35];
	I[11] = T[36];
	I[9] = T[37];
	I[5] = T[38];
	I[7] = T[39];
	I[18] = T[40];
	I[13] = T[41];
	I[24] = T[42];
	I[26] = T[43];
	I[3] = T[44];
	I[4] = T[45];
	I[16] = T[46];
	I[15] = T[47];
	I[27] = T[48];
	I[28] = T[49];
	I[42] = T[50];
	I[31] = T[51];
	I[51] = T[52];
	I[40] = T[53];
	I[30] = T[54];
	I[52] = T[55];
	I[56] = T[56];
	I[44] = T[57];
	I[33] = T[58];
	I[58] = T[59];
    }

};

template<>
struct impl<meta::braket<rysq::P, rysq::P, rysq::P, rysq::P> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[1], double (&I)[81]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];

	double xkl = rkl[0];
	double ykl = rkl[1];
	double zkl = rkl[2];

// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Kx = (xkl + Dx);
	    double Ky = (ykl + Dy);
	    double Kz = (Dz + zkl);
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f1 = (B00*(ykl + 2*Dy) + Iy*(B01 + Dy*Ky));
	    double f10 = (B00 + Dz*Iz);
	    double f12 = (Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx));
	    double f13 = (B01 + Dx*Kx);
	    double f14 = (Cy*Iy + B10);
	    double f15 = B01*B10;
	    double f16 = (B00*(yij + 2*Cy) + Ky*(Cy*Iy + B10));
	    double f18 = (B01*Cz + B00*(2*Dz + zkl) + Cz*Dz*Kz);
	    double f19 = (B00 + Iy*Ky);
	    double f2 = (Cy*Dy*Ky + B01*Cy + B00*(ykl + 2*Dy));
	    double f21 = (B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10));
	    double f23 = (Dz*Kz + B01);
	    double f24 = (B01 + Dy*Ky);
	    double f26 = (B00 + Cy*Ky);
	    double f27 = 2*pow(B00,2);
	    double f28 = (Cx*Ix + B10);
	    double f29 = (Dy*Iy + B00);
	    double f31 = (Dz*(B10 + Cz*Iz) + B00*(2*Cz + zij));
	    double f34 = (B10 + Cz*Iz);
	    double f4 = (B01*Cx + Cx*Dx*Kx + B00*(xkl + 2*Dx));
	    double f5 = (Kx*(Cx*Ix + B10) + B00*(xij + 2*Cx));
	    double f7 = (B00 + Ix*Kx);
	    double f9 = (Iz*Kz + B00);

	    I[0] += C[0]*W[a]*((B00*(xij + 2*Cx)*(xkl + 2*Dx) + f15 + f27 + B01*Cx*Ix + Dx*Kx*(Cx*Ix + B10)));
	    I[1] += C[0]*W[a]*((f15 + f27 + Dy*Ky*(Cy*Iy + B10) + B01*Cy*Iy + B00*(yij + 2*Cy)*(ykl + 2*Dy)));
	    I[2] += C[0]*W[a]*((Dz*Kz*(B10 + Cz*Iz) + f15 + f27 + B00*(2*Cz + zij)*(2*Dz + zkl) + B01*Cz*Iz));
	    I[3] += C[0]*W[a]*(Cx*(zij*(Dz*Kz + B01) + f18));
	    I[4] += C[0]*W[a]*(Cy*(zij*(Dz*Kz + B01) + f18));
	    I[5] += C[0]*W[a]*(Dy*(f31 + zkl*(B10 + Cz*Iz)));
	    I[6] += C[0]*W[a]*(Dx*(f31 + zkl*(B10 + Cz*Iz)));
	    I[7] += C[0]*W[a]*(Cz*(f4 + xij*(B01 + Dx*Kx)));
	    I[8] += C[0]*W[a]*(Cy*(f4 + xij*(B01 + Dx*Kx)));
	    I[9] += C[0]*W[a]*(Dy*Iz*(Cx*xkl + Qx));
	    I[10] += C[0]*W[a]*(Dz*Iy*(Cx*xkl + Qx));
	    I[11] += C[0]*W[a]*(f29*(Cx*xkl + Qx));
	    I[12] += C[0]*W[a]*(f10*(Cx*xkl + Qx));
	    I[13] += C[0]*W[a]*(Cy*Kz*(Dx*xij + Qx));
	    I[14] += C[0]*W[a]*(Cz*Ky*(Dx*xij + Qx));
	    I[15] += C[0]*W[a]*(f26*(Dx*xij + Qx));
	    I[16] += C[0]*W[a]*((Cz*zkl + Qz)*(Dx*xij + Qx));
	    I[17] += C[0]*W[a]*(Dy*Ix*(Cz*zkl + Qz));
	    I[18] += C[0]*W[a]*(Cz*Ix*f24);
	    I[19] += C[0]*W[a]*(Cx*Iz*f24);
	    I[20] += C[0]*W[a]*(Dx*Iz*f26);
	    I[21] += C[0]*W[a]*(Dx*Ky*f34);
	    I[22] += C[0]*W[a]*(Dz*Ky*f28);
	    I[23] += C[0]*W[a]*(Dz*Ix*f26);
	    I[24] += C[0]*W[a]*(Ix*Ky*Qz);
	    I[25] += C[0]*W[a]*(Cy*Kx*f10);
	    I[26] += C[0]*W[a]*(Dz*Kx*f14);
	    I[27] += C[0]*W[a]*(Cx*Dz*f19);
	    I[28] += C[0]*W[a]*(Cx*Ky*f10);
	    I[29] += C[0]*W[a]*(Iz*Ky*Qx);
	    I[30] += C[0]*W[a]*(Iy*Kz*Qx);
	    I[31] += C[0]*W[a]*(Dx*Kz*f14);
	    I[32] += C[0]*W[a]*(Cy*Dx*f9);
	    I[33] += C[0]*W[a]*(Cy*Dz*f7);
	    I[34] += C[0]*W[a]*(Cz*Dy*f7);
	    I[35] += C[0]*W[a]*(Dy*Kz*f28);
	    I[36] += C[0]*W[a]*(Ix*Kz*Qy);
	    I[37] += C[0]*W[a]*(Cx*Kz*f29);
	    I[38] += C[0]*W[a]*(Cx*Dy*f9);
	    I[39] += C[0]*W[a]*(Dy*Kx*f34);
	    I[40] += C[0]*W[a]*(Iy*Kx*Qz);
	    I[41] += C[0]*W[a]*(Dx*Iy*(Cz*zkl + Qz));
	    I[42] += C[0]*W[a]*(f29*(Cz*zkl + Qz));
	    I[43] += C[0]*W[a]*(Cz*Kx*f29);
	    I[44] += C[0]*W[a]*(Iz*Kx*Qy);
	    I[45] += C[0]*W[a]*(Cy*Iz*f13);
	    I[46] += C[0]*W[a]*(Cy*Ix*f23);
	    I[47] += C[0]*W[a]*(Cx*Iy*f23);
	    I[48] += C[0]*W[a]*(Cz*Iy*f13);
	    I[49] += C[0]*W[a]*(Cz*Dx*f19);
	    I[50] += C[0]*W[a]*(Dx*f16);
	    I[51] += C[0]*W[a]*(Cx*f1);
	    I[52] += C[0]*W[a]*(Cz*f1);
	    I[53] += C[0]*W[a]*(Dy*f5);
	    I[54] += C[0]*W[a]*(Dz*f5);
	    I[55] += C[0]*W[a]*(Dz*f16);
	    I[56] += C[0]*W[a]*(f10*f26);
	    I[57] += C[0]*W[a]*(Kx*f21);
	    I[58] += C[0]*W[a]*(Kx*f31);
	    I[59] += C[0]*W[a]*(Ky*f31);
	    I[60] += C[0]*W[a]*(Ky*f12);
	    I[61] += C[0]*W[a]*(Kz*f12);
	    I[62] += C[0]*W[a]*(Kz*f21);
	    I[63] += C[0]*W[a]*(Qx*f9);
	    I[64] += C[0]*W[a]*(Qx*f19);
	    I[65] += C[0]*W[a]*(Qz*f19);
	    I[66] += C[0]*W[a]*(Qz*f7);
	    I[67] += C[0]*W[a]*(Qy*f7);
	    I[68] += C[0]*W[a]*(Qy*f9);
	    I[69] += C[0]*W[a]*(f13*f14);
	    I[70] += C[0]*W[a]*(f13*f34);
	    I[71] += C[0]*W[a]*(f24*f34);
	    I[72] += C[0]*W[a]*(f24*f28);
	    I[73] += C[0]*W[a]*(f23*f28);
	    I[74] += C[0]*W[a]*(f14*f23);
	    I[75] += C[0]*W[a]*(Ix*f18);
	    I[76] += C[0]*W[a]*(Ix*f2);
	    I[77] += C[0]*W[a]*(Iz*f2);
	    I[78] += C[0]*W[a]*(Iz*f4);
	    I[79] += C[0]*W[a]*(Iy*f4);
	    I[80] += C[0]*W[a]*(Iy*f18);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[81]) {
	double T[81];
	for (int i = 0; i < 81; ++i) {
	    T[i] = I[i];
	}
	I[0] = T[0];
	I[40] = T[1];
	I[80] = T[2];
	I[78] = T[3];
	I[79] = T[4];
	I[71] = T[5];
	I[62] = T[6];
	I[2] = T[7];
	I[1] = T[8];
	I[15] = T[9];
	I[21] = T[10];
	I[12] = T[11];
	I[24] = T[12];
	I[55] = T[13];
	I[29] = T[14];
	I[28] = T[15];
	I[56] = T[16];
	I[65] = T[17];
	I[38] = T[18];
	I[42] = T[19];
	I[34] = T[20];
	I[35] = T[21];
	I[45] = T[22];
	I[46] = T[23];
	I[47] = T[24];
	I[25] = T[25];
	I[22] = T[26];
	I[48] = T[27];
	I[51] = T[28];
	I[33] = T[29];
	I[57] = T[30];
	I[58] = T[31];
	I[61] = T[32];
	I[19] = T[33];
	I[11] = T[34];
	I[63] = T[35];
	I[64] = T[36];
	I[66] = T[37];
	I[69] = T[38];
	I[17] = T[39];
	I[23] = T[40];
	I[59] = T[41];
	I[68] = T[42];
	I[14] = T[43];
	I[16] = T[44];
	I[7] = T[45];
	I[73] = T[46];
	I[75] = T[47];
	I[5] = T[48];
	I[32] = T[49];
	I[31] = T[50];
	I[39] = T[51];
	I[41] = T[52];
	I[9] = T[53];
	I[18] = T[54];
	I[49] = T[55];
	I[52] = T[56];
	I[13] = T[57];
	I[26] = T[58];
	I[53] = T[59];
	I[27] = T[60];
	I[54] = T[61];
	I[67] = T[62];
	I[60] = T[63];
	I[30] = T[64];
	I[50] = T[65];
	I[20] = T[66];
	I[10] = T[67];
	I[70] = T[68];
	I[4] = T[69];
	I[8] = T[70];
	I[44] = T[71];
	I[36] = T[72];
	I[72] = T[73];
	I[76] = T[74];
	I[74] = T[75];
	I[37] = T[76];
	I[43] = T[77];
	I[6] = T[78];
	I[3] = T[79];
	I[77] = T[80];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::P, rysq::D, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[1], double (&I)[108]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double Rx = (B01 + pow(Dx,2));
	    double Ry = (B01 + pow(Dy,2));
	    double Rz = (pow(Dz,2) + B01);
	    double f0 = (Iz*Rz + 2*B00*Dz);
	    double f1 = (Dz*(Iz*pow(Cz,2) + B10*(3*Cz + zij)) + 3*B00*Pz + 2*B00*Cz*zij);
	    double f12 = (Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx));
	    double f15 = (Cy*Iy + B10);
	    double f16 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f17 = (2*B00*Dx + Cx*Rx);
	    double f2 = (4*B00*Cy*Dy + Py*Ry + 2*pow(B00,2));
	    double f20 = (2*B00*Dx*(xij + 2*Cx) + Rx*(Cx*Ix + B10) + 2*pow(B00,2));
	    double f22 = (Px*Rx + 2*pow(B00,2) + 4*B00*Cx*Dx);
	    double f24 = (2*B00*Cy*yij + 3*B00*Py + Dy*(B10*(3*Cy + yij) + Iy*pow(Cy,2)));
	    double f26 = (Dy*Iy + B00);
	    double f29 = (Dx*Px + 2*B00*Cx);
	    double f3 = (Dy*Py + 2*B00*Cy);
	    double f30 = (2*B00*Dy + Iy*Ry);
	    double f31 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f32 = (Dx*(B10*(3*Cx + xij) + Ix*pow(Cx,2)) + 2*B00*Cx*xij + 3*B00*Px);
	    double f33 = (2*B00*Dy + Cy*Ry);
	    double f39 = (Rz*(B10 + Cz*Iz) + 2*B00*Dz*(2*Cz + zij) + 2*pow(B00,2));
	    double f40 = (2*pow(B00,2) + Pz*Rz + 4*B00*Cz*Dz);
	    double f5 = (2*B00*Dz + Cz*Rz);
	    double f6 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));
	    double f9 = (Dz*Pz + 2*B00*Cz);

	    I[0] += C[0]*W[a]*(Cy*Dx*(f9 + Qz*zij));
	    I[1] += C[0]*W[a]*((6*B00*Dz*Pz + 2*pow(B00,2)*(3*Cz + zij) + 4*B00*Cz*Dz*zij + Rz*(Iz*pow(Cz,2) + B10*(3*Cz + zij))));
	    I[2] += C[0]*W[a]*(Qx*(f9 + Qz*zij));
	    I[3] += C[0]*W[a]*(Qy*(f9 + Qz*zij));
	    I[4] += C[0]*W[a]*(Cx*Dy*(f9 + Qz*zij));
	    I[5] += C[0]*W[a]*(Cx*Dz*(Qy*yij + f3));
	    I[6] += C[0]*W[a]*(Qz*(Qy*yij + f3));
	    I[7] += C[0]*W[a]*(Qx*(Qy*yij + f3));
	    I[8] += C[0]*W[a]*(Cz*Dx*(Qy*yij + f3));
	    I[9] += C[0]*W[a]*(Cz*(yij*(2*B00*Dy + Cy*Ry) + f2));
	    I[10] += C[0]*W[a]*(Cx*(yij*(2*B00*Dy + Cy*Ry) + f2));
	    I[11] += C[0]*W[a]*((4*B00*Cy*Dy*yij + 2*pow(B00,2)*(3*Cy + yij) + Ry*(B10*(3*Cy + yij) + Iy*pow(Cy,2)) + 6*B00*Dy*Py));
	    I[12] += C[0]*W[a]*((6*B00*Dx*Px + Rx*(B10*(3*Cx + xij) + Ix*pow(Cx,2)) + 4*B00*Cx*Dx*xij + 2*pow(B00,2)*(3*Cx + xij)));
	    I[13] += C[0]*W[a]*(Dy*Pz*(Dx*xij + Qx));
	    I[14] += C[0]*W[a]*(Dx*Qy*(Cz*zij + Pz));
	    I[15] += C[0]*W[a]*(Cz*Qy*(Dx*xij + Qx));
	    I[16] += C[0]*W[a]*(Cy*Cz*(f17 + Rx*xij));
	    I[17] += C[0]*W[a]*(Pz*(f17 + Rx*xij));
	    I[18] += C[0]*W[a]*(Py*(f17 + Rx*xij));
	    I[19] += C[0]*W[a]*(Dz*Py*(Dx*xij + Qx));
	    I[20] += C[0]*W[a]*(Dx*Py*(Dz*zij + Qz));
	    I[21] += C[0]*W[a]*(Cy*Qx*(Dz*zij + Qz));
	    I[22] += C[0]*W[a]*(Cy*Qz*(Dx*xij + Qx));
	    I[23] += C[0]*W[a]*(Cy*Rz*(Px + Cx*xij));
	    I[24] += C[0]*W[a]*(Dz*Qy*(Px + Cx*xij));
	    I[25] += C[0]*W[a]*(Dy*Qz*(Px + Cx*xij));
	    I[26] += C[0]*W[a]*(Dy*Px*(Dz*zij + Qz));
	    I[27] += C[0]*W[a]*(Cx*Qy*(Dz*zij + Qz));
	    I[28] += C[0]*W[a]*(f3*(Dz*zij + Qz));
	    I[29] += C[0]*W[a]*(f29*(Dz*zij + Qz));
	    I[30] += C[0]*W[a]*(Dz*Iy*f29);
	    I[31] += C[0]*W[a]*(Dx*Dz*f6);
	    I[32] += C[0]*W[a]*(Dz*Qx*f15);
	    I[33] += C[0]*W[a]*(Iy*Qx*Qz);
	    I[34] += C[0]*W[a]*(Dx*Qz*f15);
	    I[35] += C[0]*W[a]*(Cx*Cy*f0);
	    I[36] += C[0]*W[a]*(Cy*Dz*f12);
	    I[37] += C[0]*W[a]*(Dy*Dz*f16);
	    I[38] += C[0]*W[a]*(Dz*Ix*f3);
	    I[39] += C[0]*W[a]*(Dy*Ix*f9);
	    I[40] += C[0]*W[a]*(Ix*Qy*Qz);
	    I[41] += C[0]*W[a]*(Iz*Qx*Qy);
	    I[42] += C[0]*W[a]*(Dy*Iz*f29);
	    I[43] += C[0]*W[a]*(Dx*Dy*f31);
	    I[44] += C[0]*W[a]*(Dx*Iy*f9);
	    I[45] += C[0]*W[a]*(f9*(Dx*xij + Qx));
	    I[46] += C[0]*W[a]*(f3*(Dx*xij + Qx));
	    I[47] += C[0]*W[a]*(Dx*Iz*f3);
	    I[48] += C[0]*W[a]*(Dx*Pz*f26);
	    I[49] += C[0]*W[a]*(Cz*Qx*f26);
	    I[50] += C[0]*W[a]*(Cz*Dy*f12);
	    I[51] += C[0]*W[a]*(Dy*Qx*(Cz*zij + Pz));
	    I[52] += C[0]*W[a]*(Cy*Rx*(Cz*zij + Pz));
	    I[53] += C[0]*W[a]*(Cx*Ry*(Cz*zij + Pz));
	    I[54] += C[0]*W[a]*(Cz*Ry*(Px + Cx*xij));
	    I[55] += C[0]*W[a]*(f5*(Px + Cx*xij));
	    I[56] += C[0]*W[a]*(f33*(Px + Cx*xij));
	    I[57] += C[0]*W[a]*(Cx*Iz*f33);
	    I[58] += C[0]*W[a]*(Cx*Cz*f30);
	    I[59] += C[0]*W[a]*(Cx*Rz*f15);
	    I[60] += C[0]*W[a]*(Iy*Px*Rz);
	    I[61] += C[0]*W[a]*(Dz*Px*f26);
	    I[62] += C[0]*W[a]*(Cx*Qz*f26);
	    I[63] += C[0]*W[a]*(Cx*Iy*f5);
	    I[64] += C[0]*W[a]*(Cy*Ix*f5);
	    I[65] += C[0]*W[a]*(Ix*Py*Rz);
	    I[66] += C[0]*W[a]*(Ix*Pz*Ry);
	    I[67] += C[0]*W[a]*(Iz*Px*Ry);
	    I[68] += C[0]*W[a]*(Iz*Py*Rx);
	    I[69] += C[0]*W[a]*(Iy*Pz*Rx);
	    I[70] += C[0]*W[a]*(Cz*Rx*f15);
	    I[71] += C[0]*W[a]*(Cz*Ix*f33);
	    I[72] += C[0]*W[a]*(f33*(Cz*zij + Pz));
	    I[73] += C[0]*W[a]*(f17*(Cz*zij + Pz));
	    I[74] += C[0]*W[a]*(Cz*Iy*f17);
	    I[75] += C[0]*W[a]*(Cy*Iz*f17);
	    I[76] += C[0]*W[a]*(Cy*f20);
	    I[77] += C[0]*W[a]*(Cz*f20);
	    I[78] += C[0]*W[a]*(Rx*f6);
	    I[79] += C[0]*W[a]*(Rx*f31);
	    I[80] += C[0]*W[a]*(Ry*f31);
	    I[81] += C[0]*W[a]*(Ry*f16);
	    I[82] += C[0]*W[a]*(Rz*f16);
	    I[83] += C[0]*W[a]*(Rz*f6);
	    I[84] += C[0]*W[a]*(Ix*f40);
	    I[85] += C[0]*W[a]*(Ix*f2);
	    I[86] += C[0]*W[a]*(Iz*f2);
	    I[87] += C[0]*W[a]*(Iz*f22);
	    I[88] += C[0]*W[a]*(Iy*f22);
	    I[89] += C[0]*W[a]*(Iy*f40);
	    I[90] += C[0]*W[a]*(f15*f17);
	    I[91] += C[0]*W[a]*(f15*f5);
	    I[92] += C[0]*W[a]*(Pz*f30);
	    I[93] += C[0]*W[a]*(Px*f30);
	    I[94] += C[0]*W[a]*(Px*f0);
	    I[95] += C[0]*W[a]*(Py*f0);
	    I[96] += C[0]*W[a]*(Cx*f39);
	    I[97] += C[0]*W[a]*(Cy*f39);
	    I[98] += C[0]*W[a]*(Qy*f12);
	    I[99] += C[0]*W[a]*(Qz*f12);
	    I[100] += C[0]*W[a]*(f26*f29);
	    I[101] += C[0]*W[a]*(f26*f9);
	    I[102] += C[0]*W[a]*(Dx*f1);
	    I[103] += C[0]*W[a]*(Dx*f24);
	    I[104] += C[0]*W[a]*(Dz*f24);
	    I[105] += C[0]*W[a]*(Dz*f32);
	    I[106] += C[0]*W[a]*(Dy*f32);
	    I[107] += C[0]*W[a]*(Dy*f1);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[108]) {
	double T[108];
	for (int i = 0; i < 108; ++i) {
	    T[i] = I[i];
	}
	I[89] = T[0];
	I[50] = T[1];
	I[88] = T[2];
	I[107] = T[3];
	I[106] = T[4];
	I[99] = T[5];
	I[101] = T[6];
	I[63] = T[7];
	I[65] = T[8];
	I[29] = T[9];
	I[27] = T[10];
	I[25] = T[11];
	I[0] = T[12];
	I[56] = T[13];
	I[71] = T[14];
	I[59] = T[15];
	I[5] = T[16];
	I[2] = T[17];
	I[1] = T[18];
	I[73] = T[19];
	I[85] = T[20];
	I[87] = T[21];
	I[77] = T[22];
	I[39] = T[23];
	I[93] = T[24];
	I[94] = T[25];
	I[102] = T[26];
	I[105] = T[27];
	I[103] = T[28];
	I[84] = T[29];
	I[78] = T[30];
	I[79] = T[31];
	I[81] = T[32];
	I[82] = T[33];
	I[83] = T[34];
	I[51] = T[35];
	I[75] = T[36];
	I[90] = T[37];
	I[91] = T[38];
	I[92] = T[39];
	I[95] = T[40];
	I[69] = T[41];
	I[66] = T[42];
	I[68] = T[43];
	I[80] = T[44];
	I[74] = T[45];
	I[55] = T[46];
	I[67] = T[47];
	I[62] = T[48];
	I[64] = T[49];
	I[58] = T[50];
	I[70] = T[51];
	I[17] = T[52];
	I[34] = T[53];
	I[22] = T[54];
	I[40] = T[55];
	I[21] = T[56];
	I[33] = T[57];
	I[28] = T[58];
	I[45] = T[59];
	I[42] = T[60];
	I[96] = T[61];
	I[100] = T[62];
	I[46] = T[63];
	I[41] = T[64];
	I[37] = T[65];
	I[20] = T[66];
	I[30] = T[67];
	I[13] = T[68];
	I[8] = T[69];
	I[11] = T[70];
	I[23] = T[71];
	I[35] = T[72];
	I[16] = T[73];
	I[10] = T[74];
	I[15] = T[75];
	I[3] = T[76];
	I[4] = T[77];
	I[7] = T[78];
	I[14] = T[79];
	I[32] = T[80];
	I[18] = T[81];
	I[36] = T[82];
	I[43] = T[83];
	I[38] = T[84];
	I[19] = T[85];
	I[31] = T[86];
	I[12] = T[87];
	I[6] = T[88];
	I[44] = T[89];
	I[9] = T[90];
	I[47] = T[91];
	I[26] = T[92];
	I[24] = T[93];
	I[48] = T[94];
	I[49] = T[95];
	I[52] = T[96];
	I[53] = T[97];
	I[57] = T[98];
	I[76] = T[99];
	I[60] = T[100];
	I[98] = T[101];
	I[86] = T[102];
	I[61] = T[103];
	I[97] = T[104];
	I[72] = T[105];
	I[54] = T[106];
	I[104] = T[107];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::S, rysq::F, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<4> &t2, const Vector<4> &W,
			    const double (&C)[1], double (&I)[100]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 4;



// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double Rx = (B01 + pow(Dx,2));
	    double Ry = (B01 + pow(Dy,2));
	    double Rz = (pow(Dz,2) + B01);
	    double f0 = (4*B00*Cy*Dy + Py*Ry + 2*pow(B00,2));
	    double f1 = (3*B10 + pow(Cx,2));
	    double f11 = (2*B00*Dz + Cz*Rz);
	    double f12 = (pow(Dz,2) + 3*B01);
	    double f13 = (6*Cy*pow(B00,2) + Cy*Ry*(3*B10 + pow(Cy,2)) + 6*B00*Dy*Py);
	    double f15 = (3*B01 + pow(Dy,2));
	    double f16 = 9*B00*B01*B10;
	    double f19 = 6*pow(B00,3);
	    double f2 = (6*Dz*pow(B00,2) + Dz*Pz*(pow(Dz,2) + 3*B01) + 6*B00*Cz*Rz);
	    double f21 = (3*B01 + pow(Dx,2));
	    double f22 = (3*B00*Pz + Cz*Dz*(3*B10 + pow(Cz,2)));
	    double f23 = (Dy*Py + 2*B00*Cy);
	    double f25 = (2*B00*Dx + Cx*Rx);
	    double f26 = (3*B00*Rx + Cx*Dx*(3*B01 + pow(Dx,2)));
	    double f27 = (Dx*Px + 2*B00*Cx);
	    double f28 = (Cy*Dy*(3*B01 + pow(Dy,2)) + 3*B00*Ry);
	    double f29 = (6*Cx*pow(B00,2) + Cx*Rx*(3*B10 + pow(Cx,2)) + 6*B00*Dx*Px);
	    double f3 = (Dz*Pz + 2*B00*Cz);
	    double f32 = (Cx*Dx*(3*B10 + pow(Cx,2)) + 3*B00*Px);
	    double f33 = (3*B10 + pow(Cz,2));
	    double f34 = (3*B10 + pow(Cy,2));
	    double f35 = (2*pow(B00,2) + Pz*Rz + 4*B00*Cz*Dz);
	    double f36 = (6*B00*Dz*Pz + Cz*Rz*(3*B10 + pow(Cz,2)) + 6*Cz*pow(B00,2));
	    double f4 = (2*B00*Dy + Cy*Ry);
	    double f5 = (Dx*Px*(3*B01 + pow(Dx,2)) + 6*Dx*pow(B00,2) + 6*B00*Cx*Rx);
	    double f6 = (Px*Rx + 2*pow(B00,2) + 4*B00*Cx*Dx);
	    double f7 = (6*B00*Cy*Ry + 6*Dy*pow(B00,2) + Dy*Py*(3*B01 + pow(Dy,2)));
	    double f8 = (Cy*Dy*(3*B10 + pow(Cy,2)) + 3*B00*Py);
	    double f9 = (Cz*Dz*(pow(Dz,2) + 3*B01) + 3*B00*Rz);

	    I[0] += C[0]*W[a]*((f16 + f19 + Cy*Dy*f15*f34 + 9*B00*(B01*pow(Cy,2) + Py*pow(Dy,2)) + 18*Cy*Dy*pow(B00,2)));
	    I[1] += C[0]*W[a]*((Cz*Dz*f12*f33 + f16 + f19 + 9*B00*(B01*pow(Cz,2) + Pz*pow(Dz,2)) + 18*Cz*Dz*pow(B00,2)));
	    I[2] += C[0]*W[a]*((18*Cx*Dx*pow(B00,2) + f16 + f19 + Cx*Dx*f1*f21 + 9*B00*(B01*pow(Cx,2) + Px*pow(Dx,2))));
	    I[3] += C[0]*W[a]*(Cx*Ry*f3);
	    I[4] += C[0]*W[a]*(Cx*Qz*f4);
	    I[5] += C[0]*W[a]*(Cz*Qx*f4);
	    I[6] += C[0]*W[a]*(Py*Qx*Rz);
	    I[7] += C[0]*W[a]*(Cy*Rz*f27);
	    I[8] += C[0]*W[a]*(Cy*Dx*Rz*f34);
	    I[9] += C[0]*W[a]*(Cy*Dx*f35);
	    I[10] += C[0]*W[a]*(Cy*Qx*f11);
	    I[11] += C[0]*W[a]*(Pz*Qx*Ry);
	    I[12] += C[0]*W[a]*(Cz*Dx*Ry*f33);
	    I[13] += C[0]*W[a]*(Cz*Ry*f27);
	    I[14] += C[0]*W[a]*(Px*Qz*Ry);
	    I[15] += C[0]*W[a]*(Px*Qy*Rz);
	    I[16] += C[0]*W[a]*(Cx*Qy*f11);
	    I[17] += C[0]*W[a]*(Cx*Dy*f35);
	    I[18] += C[0]*W[a]*(Dy*Dz*f32);
	    I[19] += C[0]*W[a]*(Dx*Dz*f8);
	    I[20] += C[0]*W[a]*(Dx*Dy*f22);
	    I[21] += C[0]*W[a]*(Dx*Pz*f4);
	    I[22] += C[0]*W[a]*(Cy*Dx*Pz*f21);
	    I[23] += C[0]*W[a]*(Cy*Dx*f21*f34);
	    I[24] += C[0]*W[a]*(Cy*Dz*Rx*f34);
	    I[25] += C[0]*W[a]*(Cy*Rx*f3);
	    I[26] += C[0]*W[a]*(Cy*Qz*f25);
	    I[27] += C[0]*W[a]*(Cz*Qy*f25);
	    I[28] += C[0]*W[a]*(Pz*Qy*Rx);
	    I[29] += C[0]*W[a]*(Dz*Qy*f27);
	    I[30] += C[0]*W[a]*(Dy*Qz*f27);
	    I[31] += C[0]*W[a]*(Cz*Dy*Rx*f33);
	    I[32] += C[0]*W[a]*(Cz*Rx*f23);
	    I[33] += C[0]*W[a]*(Dz*Qx*f23);
	    I[34] += C[0]*W[a]*(Dx*Qz*f23);
	    I[35] += C[0]*W[a]*(Cx*Rz*f23);
	    I[36] += C[0]*W[a]*(Cx*Dy*Rz*f1);
	    I[37] += C[0]*W[a]*(Cx*Dy*f1*f15);
	    I[38] += C[0]*W[a]*(Cx*Dy*Pz*f15);
	    I[39] += C[0]*W[a]*(Dy*Pz*f25);
	    I[40] += C[0]*W[a]*(Dy*Px*f11);
	    I[41] += C[0]*W[a]*(Cz*Dy*Px*f15);
	    I[42] += C[0]*W[a]*(Cz*Dy*f15*f33);
	    I[43] += C[0]*W[a]*(Cz*Dy*f6);
	    I[44] += C[0]*W[a]*(Dy*Qx*f3);
	    I[45] += C[0]*W[a]*(Dx*Qy*f3);
	    I[46] += C[0]*W[a]*(Qx*Qy*Qz);
	    I[47] += C[0]*W[a]*(Py*Qz*Rx);
	    I[48] += C[0]*W[a]*(Dx*Py*f11);
	    I[49] += C[0]*W[a]*(Cz*Dx*Py*f21);
	    I[50] += C[0]*W[a]*(Cz*Dx*f21*f33);
	    I[51] += C[0]*W[a]*(Cz*Dx*f0);
	    I[52] += C[0]*W[a]*(Cx*Dz*f0);
	    I[53] += C[0]*W[a]*(Cx*Dz*Ry*f1);
	    I[54] += C[0]*W[a]*(Cx*Dz*f1*f12);
	    I[55] += C[0]*W[a]*(Cx*Dz*Py*f12);
	    I[56] += C[0]*W[a]*(Dz*Py*f25);
	    I[57] += C[0]*W[a]*(Dz*Px*f4);
	    I[58] += C[0]*W[a]*(Cy*Dz*Px*f12);
	    I[59] += C[0]*W[a]*(Cy*Dz*f12*f34);
	    I[60] += C[0]*W[a]*(Cy*Dz*f6);
	    I[61] += C[0]*W[a]*(Cx*Cy*f9);
	    I[62] += C[0]*W[a]*(Cx*Cz*f28);
	    I[63] += C[0]*W[a]*(Cy*Cz*f26);
	    I[64] += C[0]*W[a]*(Py*f26);
	    I[65] += C[0]*W[a]*(Pz*f26);
	    I[66] += C[0]*W[a]*(Pz*f28);
	    I[67] += C[0]*W[a]*(Px*f28);
	    I[68] += C[0]*W[a]*(Px*f9);
	    I[69] += C[0]*W[a]*(Py*f9);
	    I[70] += C[0]*W[a]*(Cy*f5);
	    I[71] += C[0]*W[a]*(Cz*f5);
	    I[72] += C[0]*W[a]*(Cz*f7);
	    I[73] += C[0]*W[a]*(Cx*f7);
	    I[74] += C[0]*W[a]*(Cx*f2);
	    I[75] += C[0]*W[a]*(Cy*f2);
	    I[76] += C[0]*W[a]*(Rx*f8);
	    I[77] += C[0]*W[a]*(Rx*f22);
	    I[78] += C[0]*W[a]*(Ry*f22);
	    I[79] += C[0]*W[a]*(Ry*f32);
	    I[80] += C[0]*W[a]*(Rz*f32);
	    I[81] += C[0]*W[a]*(Rz*f8);
	    I[82] += C[0]*W[a]*(Dx*f36);
	    I[83] += C[0]*W[a]*(Dx*f13);
	    I[84] += C[0]*W[a]*(Dz*f13);
	    I[85] += C[0]*W[a]*(Dz*f29);
	    I[86] += C[0]*W[a]*(Dy*f29);
	    I[87] += C[0]*W[a]*(Dy*f36);
	    I[88] += C[0]*W[a]*(f11*f27);
	    I[89] += C[0]*W[a]*(f27*f4);
	    I[90] += C[0]*W[a]*(f3*f4);
	    I[91] += C[0]*W[a]*(f25*f3);
	    I[92] += C[0]*W[a]*(f23*f25);
	    I[93] += C[0]*W[a]*(f11*f23);
	    I[94] += C[0]*W[a]*(Qx*f35);
	    I[95] += C[0]*W[a]*(Qx*f0);
	    I[96] += C[0]*W[a]*(Qz*f0);
	    I[97] += C[0]*W[a]*(Qz*f6);
	    I[98] += C[0]*W[a]*(Qy*f6);
	    I[99] += C[0]*W[a]*(Qy*f35);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[100]) {
	double T[100];
	for (int i = 0; i < 100; ++i) {
	    T[i] = I[i];
	}
	I[11] = T[0];
	I[22] = T[1];
	I[0] = T[2];
	I[67] = T[3];
	I[69] = T[4];
	I[59] = T[5];
	I[75] = T[6];
	I[73] = T[7];
	I[71] = T[8];
	I[78] = T[9];
	I[79] = T[10];
	I[57] = T[11];
	I[52] = T[12];
	I[54] = T[13];
	I[64] = T[14];
	I[83] = T[15];
	I[89] = T[16];
	I[87] = T[17];
	I[90] = T[18];
	I[91] = T[19];
	I[92] = T[20];
	I[58] = T[21];
	I[8] = T[22];
	I[1] = T[23];
	I[41] = T[24];
	I[48] = T[25];
	I[49] = T[26];
	I[39] = T[27];
	I[38] = T[28];
	I[93] = T[29];
	I[94] = T[30];
	I[32] = T[31];
	I[36] = T[32];
	I[95] = T[33];
	I[96] = T[34];
	I[85] = T[35];
	I[80] = T[36];
	I[10] = T[37];
	I[17] = T[38];
	I[37] = T[39];
	I[84] = T[40];
	I[14] = T[41];
	I[12] = T[42];
	I[34] = T[43];
	I[97] = T[44];
	I[98] = T[45];
	I[99] = T[46];
	I[46] = T[47];
	I[76] = T[48];
	I[6] = T[49];
	I[2] = T[50];
	I[56] = T[51];
	I[65] = T[52];
	I[60] = T[53];
	I[20] = T[54];
	I[25] = T[55];
	I[45] = T[56];
	I[63] = T[57];
	I[23] = T[58];
	I[21] = T[59];
	I[43] = T[60];
	I[29] = T[61];
	I[19] = T[62];
	I[9] = T[63];
	I[5] = T[64];
	I[7] = T[65];
	I[18] = T[66];
	I[13] = T[67];
	I[24] = T[68];
	I[26] = T[69];
	I[3] = T[70];
	I[4] = T[71];
	I[16] = T[72];
	I[15] = T[73];
	I[27] = T[74];
	I[28] = T[75];
	I[31] = T[76];
	I[42] = T[77];
	I[62] = T[78];
	I[50] = T[79];
	I[70] = T[80];
	I[81] = T[81];
	I[72] = T[82];
	I[51] = T[83];
	I[61] = T[84];
	I[40] = T[85];
	I[30] = T[86];
	I[82] = T[87];
	I[74] = T[88];
	I[53] = T[89];
	I[68] = T[90];
	I[47] = T[91];
	I[35] = T[92];
	I[86] = T[93];
	I[77] = T[94];
	I[55] = T[95];
	I[66] = T[96];
	I[44] = T[97];
	I[33] = T[98];
	I[88] = T[99];
    }

};

template<>
struct impl<meta::braket<rysq::SP, rysq::P, rysq::P, rysq::P> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[2], double (&I)[108]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];

	double xkl = rkl[0];
	double ykl = rkl[1];
	double zkl = rkl[2];

// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB01[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 
// 	    vB01[a] = recurrence::coefficient(1.0/B, A, t2[a]);
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);
	    double B01 = recurrence::coefficient(1.0/B, A, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Kx = (xkl + Dx);
	    double Ky = (ykl + Dy);
	    double Kz = (Dz + zkl);
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f1 = B01*B10;
	    double f11 = (B00 + Iy*Ky);
	    double f12 = (B00 + Cy*Ky);
	    double f13 = (Dy*Iy + B00);
	    double f14 = (Dz*(B10 + Cz*Iz) + B00*(2*Cz + zij));
	    double f15 = (Cx*Ix + B10);
	    double f16 = (B00*(2*Dz + zkl) + Iz*(Dz*Kz + B01));
	    double f18 = (Kx*(Cx*Ix + B10) + B00*(xij + 2*Cx));
	    double f2 = (Dz*Kz + B01);
	    double f20 = (Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx));
	    double f21 = (B01 + Dx*Kx);
	    double f22 = (Cy*Iy + B10);
	    double f23 = (B00*(yij + 2*Cy) + Ky*(Cy*Iy + B10));
	    double f24 = (B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10));
	    double f26 = (B01 + Dy*Ky);
	    double f28 = (Cy*Dy*Ky + B01*Cy + B00*(ykl + 2*Dy));
	    double f29 = (B01*Cx + Cx*Dx*Kx + B00*(xkl + 2*Dx));
	    double f3 = (B00 + Cx*Kx);
	    double f30 = (B01*Cz + B00*(2*Dz + zkl) + Cz*Dz*Kz);
	    double f31 = 2*pow(B00,2);
	    double f36 = (B10 + Cz*Iz);
	    double f4 = (B00*(ykl + 2*Dy) + Iy*(B01 + Dy*Ky));
	    double f7 = (B00 + Dz*Iz);

	    I[0] += C[1]*W[a]*((B00*(xij + 2*Cx)*(xkl + 2*Dx) + f31 + B01*Cx*Ix + f1 + Dx*Kx*(Cx*Ix + B10)));
	    I[1] += C[1]*W[a]*((f31 + Dy*Ky*(Cy*Iy + B10) + B01*Cy*Iy + B00*(yij + 2*Cy)*(ykl + 2*Dy) + f1));
	    I[2] += C[1]*W[a]*((Dz*Kz*(B10 + Cz*Iz) + f31 + B00*(2*Cz + zij)*(2*Dz + zkl) + f1 + B01*Cz*Iz));
	    I[3] += C[1]*W[a]*(Dy*Ix*(Cz*zkl + Qz));
	    I[4] += C[1]*W[a]*(Dy*(f14 + zkl*(B10 + Cz*Iz)));
	    I[5] += C[1]*W[a]*(Dx*(f14 + zkl*(B10 + Cz*Iz)));
	    I[6] += C[1]*W[a]*(Cy*Dx*(Iz*zkl + f7));
	    I[7] += C[1]*W[a]*(Qx*(Iz*zkl + f7));
	    I[8] += C[1]*W[a]*(Cx*Dy*(Iz*zkl + f7));
	    I[9] += C[0]*W[a]*(Dy*(Iz*zkl + f7));
	    I[10] += C[1]*W[a]*(Qy*(Iz*zkl + f7));
	    I[11] += C[0]*W[a]*(Dx*(Iz*zkl + f7));
	    I[12] += C[1]*W[a]*(Cz*(f29 + xij*(B01 + Dx*Kx)));
	    I[13] += C[0]*W[a]*((f29 + xij*(B01 + Dx*Kx)));
	    I[14] += C[1]*W[a]*(Cy*(f29 + xij*(B01 + Dx*Kx)));
	    I[15] += C[1]*W[a]*(Cy*Dz*(Kx*xij + f3));
	    I[16] += C[0]*W[a]*(Dz*(Kx*xij + f3));
	    I[17] += C[1]*W[a]*(Qy*(Kx*xij + f3));
	    I[18] += C[1]*W[a]*(Cz*Dy*(Kx*xij + f3));
	    I[19] += C[0]*W[a]*(Dy*(Kx*xij + f3));
	    I[20] += C[1]*W[a]*(Qz*(Kx*xij + f3));
	    I[21] += C[1]*W[a]*((Cz*zkl + Qz)*(Dx*xij + Qx));
	    I[22] += C[1]*W[a]*(Dx*Iy*(Cz*zkl + Qz));
	    I[23] += C[1]*W[a]*(f13*(Cz*zkl + Qz));
	    I[24] += C[1]*W[a]*(Cx*Dz*f11);
	    I[25] += C[1]*W[a]*(Dx*Kz*f22);
	    I[26] += C[1]*W[a]*(Ix*Kz*Qy);
	    I[27] += C[1]*W[a]*(Cx*Kz*f13);
	    I[28] += C[1]*W[a]*(Dy*Kz*f15);
	    I[29] += C[0]*W[a]*(Dy*Ix*Kz);
	    I[30] += C[1]*W[a]*(Cy*Ix*f2);
	    I[31] += C[1]*W[a]*(Dz*Ix*f12);
	    I[32] += C[1]*W[a]*(Dx*Iz*f12);
	    I[33] += C[1]*W[a]*(Dx*Ky*f36);
	    I[34] += C[1]*W[a]*(Dz*Ky*f15);
	    I[35] += C[0]*W[a]*(Dz*Ix*Ky);
	    I[36] += C[1]*W[a]*(Ix*Ky*Qz);
	    I[37] += C[1]*W[a]*(Iy*Kx*Qz);
	    I[38] += C[1]*W[a]*(Dz*Kx*f22);
	    I[39] += C[1]*W[a]*(Dy*Kx*f36);
	    I[40] += C[1]*W[a]*(Iz*Kx*Qy);
	    I[41] += C[0]*W[a]*(Dy*Iz*Kx);
	    I[42] += C[1]*W[a]*(Dy*Iz*f3);
	    I[43] += C[1]*W[a]*(Dz*Iy*f3);
	    I[44] += C[0]*W[a]*(Dz*Iy*Kx);
	    I[45] += C[1]*W[a]*(Iy*Kz*Qx);
	    I[46] += C[0]*W[a]*(Dx*Iy*Kz);
	    I[47] += C[1]*W[a]*(Cy*Kz*(Dx*xij + Qx));
	    I[48] += C[0]*W[a]*(Kz*(Dx*xij + Qx));
	    I[49] += C[1]*W[a]*(f12*(Dx*xij + Qx));
	    I[50] += C[1]*W[a]*(Cz*Ky*(Dx*xij + Qx));
	    I[51] += C[0]*W[a]*(Ky*(Dx*xij + Qx));
	    I[52] += C[1]*W[a]*(Iz*Ky*Qx);
	    I[53] += C[0]*W[a]*(Dx*Iz*Ky);
	    I[54] += C[1]*W[a]*(Dx*f23);
	    I[55] += C[1]*W[a]*(Ix*f28);
	    I[56] += C[1]*W[a]*(Cx*Iz*f26);
	    I[57] += C[0]*W[a]*(Iz*f26);
	    I[58] += C[1]*W[a]*(Iz*f28);
	    I[59] += C[1]*W[a]*(f21*f36);
	    I[60] += C[1]*W[a]*(f26*f36);
	    I[61] += C[1]*W[a]*(Dy*f18);
	    I[62] += C[1]*W[a]*(Dz*f18);
	    I[63] += C[1]*W[a]*(Dz*f23);
	    I[64] += C[0]*W[a]*(Dz*f11);
	    I[65] += C[1]*W[a]*(Qx*f11);
	    I[66] += C[1]*W[a]*(Cz*Dx*f11);
	    I[67] += C[0]*W[a]*(Dx*f11);
	    I[68] += C[1]*W[a]*(Qz*f11);
	    I[69] += C[1]*W[a]*(Kx*f24);
	    I[70] += C[1]*W[a]*(Cy*Kx*f7);
	    I[71] += C[0]*W[a]*(Kx*f7);
	    I[72] += C[1]*W[a]*(Cx*Ky*f7);
	    I[73] += C[0]*W[a]*(Ky*f7);
	    I[74] += C[1]*W[a]*(f12*f7);
	    I[75] += C[1]*W[a]*(f3*f7);
	    I[76] += C[1]*W[a]*(f13*f3);
	    I[77] += C[1]*W[a]*(Cz*Kx*f13);
	    I[78] += C[0]*W[a]*(Kx*f13);
	    I[79] += C[1]*W[a]*(Kx*f14);
	    I[80] += C[1]*W[a]*(Ky*f14);
	    I[81] += C[1]*W[a]*(Ky*f20);
	    I[82] += C[1]*W[a]*(Kz*f20);
	    I[83] += C[0]*W[a]*(Kz*f13);
	    I[84] += C[1]*W[a]*(Kz*f24);
	    I[85] += C[1]*W[a]*(f15*f2);
	    I[86] += C[1]*W[a]*(f15*f26);
	    I[87] += C[1]*W[a]*(Cz*Ix*f26);
	    I[88] += C[0]*W[a]*(Ix*f26);
	    I[89] += C[1]*W[a]*(Ix*f30);
	    I[90] += C[0]*W[a]*(Ix*f2);
	    I[91] += C[1]*W[a]*(Cx*Iy*f2);
	    I[92] += C[0]*W[a]*(Iy*f2);
	    I[93] += C[1]*W[a]*(f2*f22);
	    I[94] += C[1]*W[a]*(f21*f22);
	    I[95] += C[1]*W[a]*(Cy*Iz*f21);
	    I[96] += C[0]*W[a]*(Iz*f21);
	    I[97] += C[1]*W[a]*(Iz*f29);
	    I[98] += C[1]*W[a]*(Iy*f29);
	    I[99] += C[1]*W[a]*(Cz*Iy*f21);
	    I[100] += C[0]*W[a]*(Iy*f21);
	    I[101] += C[1]*W[a]*(Iy*f30);
	    I[102] += C[1]*W[a]*(Cx*f16);
	    I[103] += C[1]*W[a]*(Cy*f16);
	    I[104] += C[1]*W[a]*(Cz*f4);
	    I[105] += C[1]*W[a]*(Cx*f4);
	    I[106] += C[0]*W[a]*(f4);
	    I[107] += C[0]*W[a]*(f16);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[108]) {
	double T[108];
	for (int i = 0; i < 108; ++i) {
	    T[i] = I[i];
	}
	I[1] = T[0];
	I[54] = T[1];
	I[107] = T[2];
	I[87] = T[3];
	I[95] = T[4];
	I[83] = T[5];
	I[82] = T[6];
	I[81] = T[7];
	I[93] = T[8];
	I[92] = T[9];
	I[94] = T[10];
	I[80] = T[11];
	I[3] = T[12];
	I[0] = T[13];
	I[2] = T[14];
	I[26] = T[15];
	I[24] = T[16];
	I[14] = T[17];
	I[15] = T[18];
	I[12] = T[19];
	I[27] = T[20];
	I[75] = T[21];
	I[79] = T[22];
	I[91] = T[23];
	I[65] = T[24];
	I[78] = T[25];
	I[86] = T[26];
	I[89] = T[27];
	I[85] = T[28];
	I[84] = T[29];
	I[98] = T[30];
	I[62] = T[31];
	I[46] = T[32];
	I[47] = T[33];
	I[61] = T[34];
	I[60] = T[35];
	I[63] = T[36];
	I[31] = T[37];
	I[30] = T[38];
	I[23] = T[39];
	I[22] = T[40];
	I[20] = T[41];
	I[21] = T[42];
	I[29] = T[43];
	I[28] = T[44];
	I[77] = T[45];
	I[76] = T[46];
	I[74] = T[47];
	I[72] = T[48];
	I[38] = T[49];
	I[39] = T[50];
	I[36] = T[51];
	I[45] = T[52];
	I[44] = T[53];
	I[42] = T[54];
	I[50] = T[55];
	I[57] = T[56];
	I[56] = T[57];
	I[58] = T[58];
	I[11] = T[59];
	I[59] = T[60];
	I[13] = T[61];
	I[25] = T[62];
	I[66] = T[63];
	I[64] = T[64];
	I[41] = T[65];
	I[43] = T[66];
	I[40] = T[67];
	I[67] = T[68];
	I[18] = T[69];
	I[34] = T[70];
	I[32] = T[71];
	I[69] = T[72];
	I[68] = T[73];
	I[70] = T[74];
	I[33] = T[75];
	I[17] = T[76];
	I[19] = T[77];
	I[16] = T[78];
	I[35] = T[79];
	I[71] = T[80];
	I[37] = T[81];
	I[73] = T[82];
	I[88] = T[83];
	I[90] = T[84];
	I[97] = T[85];
	I[49] = T[86];
	I[51] = T[87];
	I[48] = T[88];
	I[99] = T[89];
	I[96] = T[90];
	I[101] = T[91];
	I[100] = T[92];
	I[102] = T[93];
	I[6] = T[94];
	I[10] = T[95];
	I[8] = T[96];
	I[9] = T[97];
	I[5] = T[98];
	I[7] = T[99];
	I[4] = T[100];
	I[103] = T[101];
	I[105] = T[102];
	I[106] = T[103];
	I[55] = T[104];
	I[53] = T[105];
	I[52] = T[106];
	I[104] = T[107];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::D, rysq::S, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[1], double (&I)[36]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 

// 	for (int a = 0; a < N; ++a) {
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);


	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double f1 = (B10*(3*Cx + 2*xij) + Cx*pow(Ix,2));
	    double f10 = (B10 + pow(Iz,2));
	    double f11 = 3*pow(B10,2);
	    double f15 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f3 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));
	    double f6 = (B10*(3*Cz + 2*zij) + Cz*pow(Iz,2));
	    double f7 = (Cy*pow(Iy,2) + B10*(3*Cy + 2*yij));
	    double f8 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f9 = (B10 + pow(Iy,2));

	    I[0] += C[0]*W[a]*((Py + Cy*yij)*(Cz*zij + Pz));
	    I[1] += C[0]*W[a]*(Cz*Ix*(Py + Cy*yij));
	    I[2] += C[0]*W[a]*(Cx*Iz*(Py + Cy*yij));
	    I[3] += C[0]*W[a]*(Cy*Iz*(Px + Cx*xij));
	    I[4] += C[0]*W[a]*((Px + Cx*xij)*(Cz*zij + Pz));
	    I[5] += C[0]*W[a]*(Cx*Iy*(Cz*zij + Pz));
	    I[6] += C[0]*W[a]*(Cz*Iy*(Px + Cx*xij));
	    I[7] += C[0]*W[a]*(Cy*Cz*(xij*(xij + 2*Cx) + Px));
	    I[8] += C[0]*W[a]*((Px + Cx*xij)*(Py + Cy*yij));
	    I[9] += C[0]*W[a]*(Py*(xij*(xij + 2*Cx) + Px));
	    I[10] += C[0]*W[a]*(Pz*(xij*(xij + 2*Cx) + Px));
	    I[11] += C[0]*W[a]*((f11 + B10*(6*Cx*xij + pow(xij,2) + 6*pow(Cx,2)) + pow(Cx,2)*pow(Ix,2)));
	    I[12] += C[0]*W[a]*((B10*(6*Cy*yij + pow(yij,2) + 6*pow(Cy,2)) + f11 + pow(Cy,2)*pow(Iy,2)));
	    I[13] += C[0]*W[a]*((f11 + pow(Cz,2)*pow(Iz,2) + B10*(6*Cz*zij + pow(zij,2) + 6*pow(Cz,2))));
	    I[14] += C[0]*W[a]*(Cy*Ix*(Cz*zij + Pz));
	    I[15] += C[0]*W[a]*(Ix*Iy*Pz);
	    I[16] += C[0]*W[a]*(Ix*Iz*Py);
	    I[17] += C[0]*W[a]*(Iy*Iz*Px);
	    I[18] += C[0]*W[a]*(Cx*Cy*f10);
	    I[19] += C[0]*W[a]*(Cx*Cz*f9);
	    I[20] += C[0]*W[a]*(Pz*f9);
	    I[21] += C[0]*W[a]*(Px*f9);
	    I[22] += C[0]*W[a]*(Px*f10);
	    I[23] += C[0]*W[a]*(Py*f10);
	    I[24] += C[0]*W[a]*(Cy*f1);
	    I[25] += C[0]*W[a]*(Cz*f1);
	    I[26] += C[0]*W[a]*(Cz*f7);
	    I[27] += C[0]*W[a]*(Cx*f7);
	    I[28] += C[0]*W[a]*(Cx*f6);
	    I[29] += C[0]*W[a]*(Cy*f6);
	    I[30] += C[0]*W[a]*(Ix*f15);
	    I[31] += C[0]*W[a]*(Ix*f3);
	    I[32] += C[0]*W[a]*(Iz*f3);
	    I[33] += C[0]*W[a]*(Iz*f8);
	    I[34] += C[0]*W[a]*(Iy*f8);
	    I[35] += C[0]*W[a]*(Iy*f15);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[36]) {
	double T[36];
	for (int i = 0; i < 36; ++i) {
	    T[i] = I[i];
	}
	I[35] = T[0];
	I[23] = T[1];
	I[33] = T[2];
	I[27] = T[3];
	I[28] = T[4];
	I[34] = T[5];
	I[22] = T[6];
	I[5] = T[7];
	I[21] = T[8];
	I[1] = T[9];
	I[2] = T[10];
	I[0] = T[11];
	I[7] = T[12];
	I[14] = T[13];
	I[29] = T[14];
	I[20] = T[15];
	I[25] = T[16];
	I[30] = T[17];
	I[15] = T[18];
	I[10] = T[19];
	I[8] = T[20];
	I[6] = T[21];
	I[12] = T[22];
	I[13] = T[23];
	I[3] = T[24];
	I[4] = T[25];
	I[11] = T[26];
	I[9] = T[27];
	I[16] = T[28];
	I[17] = T[29];
	I[26] = T[30];
	I[19] = T[31];
	I[31] = T[32];
	I[24] = T[33];
	I[18] = T[34];
	I[32] = T[35];
    }

};

template<>
struct impl<meta::braket<rysq::F, rysq::P, rysq::SP, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<3> &t2, const Vector<3> &W,
			    const double (&C)[2], double (&I)[120]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 3;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f0 = (3*B10 + pow(Cx,2));
	    double f10 = (3*B00*Pz + Cz*Dz*(3*B10 + pow(Cz,2)));
	    double f15 = (Dy*Iy + B00);
	    double f16 = (3*pow(B10,2) + 3*B10*Cz*(2*Cz + zij) + Iz*pow(Cz,3));
	    double f17 = (Cy*Iy + B10);
	    double f18 = (3*pow(B10,2) + Iy*pow(Cy,3) + 3*B10*Cy*(yij + 2*Cy));
	    double f19 = (B10*(3*Cx + xij) + Ix*pow(Cx,2));
	    double f2 = (Dz*Pz + 2*B00*Cz);
	    double f21 = (Dx*(B10*(3*Cx + xij) + Ix*pow(Cx,2)) + 2*B00*Cx*xij + 3*B00*Px);
	    double f24 = 3*pow(B10,2);
	    double f28 = (Cx*Dx*(3*B10 + pow(Cx,2)) + 3*B00*Px);
	    double f3 = (Dz*(Iz*pow(Cz,2) + B10*(3*Cz + zij)) + 3*B00*Pz + 2*B00*Cz*zij);
	    double f32 = (3*B10 + pow(Cz,2));
	    double f33 = (3*B10 + pow(Cy,2));
	    double f34 = (Iz*pow(Cz,2) + B10*(3*Cz + zij));
	    double f4 = (B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10));
	    double f5 = (Dx*Px + 2*B00*Cx);
	    double f6 = (B00 + Dz*Iz);
	    double f7 = (Cy*Dy*(3*B10 + pow(Cy,2)) + 3*B00*Py);
	    double f8 = (Dy*Py + 2*B00*Cy);
	    double f9 = (B10*(3*Cy + yij) + Iy*pow(Cy,2));

	    I[0] += C[1]*W[a]*(Cx*Cy*(f2 + Qz*zij));
	    I[1] += C[1]*W[a]*(Cz*(f7 + yij*(Dy*Py + 2*B00*Cy)));
	    I[2] += C[1]*W[a]*(Cx*(f7 + yij*(Dy*Py + 2*B00*Cy)));
	    I[3] += C[1]*W[a]*((Dy*(f24 + Iy*pow(Cy,3) + 3*B10*Cy*(yij + 2*Cy)) + 3*B00*Py*yij + 4*B00*Cy*f33));
	    I[4] += C[1]*W[a]*((Dz*(f24 + 3*B10*Cz*(2*Cz + zij) + Iz*pow(Cz,3)) + 3*B00*Pz*zij + 4*B00*Cz*f32));
	    I[5] += C[1]*W[a]*(Px*(f2 + Qz*zij));
	    I[6] += C[1]*W[a]*(Py*(f2 + Qz*zij));
	    I[7] += C[1]*W[a]*(Cz*f32*(Dx*xij + Qx));
	    I[8] += C[1]*W[a]*(Cz*Py*(Dx*xij + Qx));
	    I[9] += C[1]*W[a]*(Cy*f33*(Dx*xij + Qx));
	    I[10] += C[1]*W[a]*(Cy*Pz*(Dx*xij + Qx));
	    I[11] += C[1]*W[a]*(Cy*Qx*(Cz*zij + Pz));
	    I[12] += C[0]*W[a]*(Cx*Cy*(Cz*zij + Pz));
	    I[13] += C[1]*W[a]*(Cx*Qy*(Cz*zij + Pz));
	    I[14] += C[1]*W[a]*(Cz*Qy*(Px + Cx*xij));
	    I[15] += C[1]*W[a]*(Cy*Qz*(Px + Cx*xij));
	    I[16] += C[0]*W[a]*(Cy*Cz*(Px + Cx*xij));
	    I[17] += C[1]*W[a]*(Cy*Cz*(Qx*xij + f5));
	    I[18] += C[1]*W[a]*(Pz*(Qx*xij + f5));
	    I[19] += C[1]*W[a]*(Py*(Qx*xij + f5));
	    I[20] += C[1]*W[a]*((Dx*(f24 + 3*B10*Cx*(xij + 2*Cx) + Ix*pow(Cx,3)) + 3*B00*Px*xij + 4*B00*Cx*f0));
	    I[21] += C[1]*W[a]*(Dz*(f24 + 3*B10*Cx*(xij + 2*Cx) + Ix*pow(Cx,3)));
	    I[22] += C[1]*W[a]*(Dy*(f24 + 3*B10*Cx*(xij + 2*Cx) + Ix*pow(Cx,3)));
	    I[23] += C[0]*W[a]*((f24 + 3*B10*Cx*(xij + 2*Cx) + Ix*pow(Cx,3)));
	    I[24] += C[1]*W[a]*(Cy*Ix*f2);
	    I[25] += C[1]*W[a]*(Ix*Py*Qz);
	    I[26] += C[1]*W[a]*(Cx*Qz*f17);
	    I[27] += C[1]*W[a]*(Iy*Px*Qz);
	    I[28] += C[1]*W[a]*(Iz*Px*Qy);
	    I[29] += C[1]*W[a]*(Iz*Py*Qx);
	    I[30] += C[1]*W[a]*(Cy*Iz*f5);
	    I[31] += C[1]*W[a]*(Cy*f33*f6);
	    I[32] += C[1]*W[a]*(Cy*Px*f6);
	    I[33] += C[0]*W[a]*(Cy*Iz*Px);
	    I[34] += C[1]*W[a]*(Cy*Dx*Iz*f33);
	    I[35] += C[0]*W[a]*(Cy*Iz*f33);
	    I[36] += C[1]*W[a]*(Cx*Iz*f8);
	    I[37] += C[1]*W[a]*(Cx*Cz*f4);
	    I[38] += C[1]*W[a]*(Cx*Iy*f2);
	    I[39] += C[1]*W[a]*(f2*(Px + Cx*xij));
	    I[40] += C[1]*W[a]*(f8*(Px + Cx*xij));
	    I[41] += C[1]*W[a]*(Dz*Py*(Px + Cx*xij));
	    I[42] += C[0]*W[a]*(Py*(Px + Cx*xij));
	    I[43] += C[1]*W[a]*(Dy*Pz*(Px + Cx*xij));
	    I[44] += C[0]*W[a]*(Pz*(Px + Cx*xij));
	    I[45] += C[1]*W[a]*(Cx*Pz*f15);
	    I[46] += C[1]*W[a]*(Cx*f0*f15);
	    I[47] += C[1]*W[a]*(Cx*f0*f6);
	    I[48] += C[1]*W[a]*(Cx*Py*f6);
	    I[49] += C[0]*W[a]*(Cx*Iz*Py);
	    I[50] += C[1]*W[a]*(Cx*Dy*Iz*f0);
	    I[51] += C[0]*W[a]*(Cx*Iz*f0);
	    I[52] += C[1]*W[a]*(Cx*Dz*Iy*f0);
	    I[53] += C[0]*W[a]*(Cx*Iy*f0);
	    I[54] += C[0]*W[a]*(Cx*Iy*Pz);
	    I[55] += C[1]*W[a]*(Iy*Pz*Qx);
	    I[56] += C[1]*W[a]*(Cz*Qx*f17);
	    I[57] += C[0]*W[a]*(Cx*Cz*f17);
	    I[58] += C[1]*W[a]*(Cz*f15*f32);
	    I[59] += C[1]*W[a]*(Cz*Px*f15);
	    I[60] += C[0]*W[a]*(Cz*Iy*Px);
	    I[61] += C[1]*W[a]*(Cz*Dx*Iy*f32);
	    I[62] += C[0]*W[a]*(Cz*Iy*f32);
	    I[63] += C[1]*W[a]*(Cz*Iy*f5);
	    I[64] += C[1]*W[a]*(f5*(Cz*zij + Pz));
	    I[65] += C[1]*W[a]*(Dx*Py*(Cz*zij + Pz));
	    I[66] += C[0]*W[a]*(Py*(Cz*zij + Pz));
	    I[67] += C[1]*W[a]*(Dy*Px*(Cz*zij + Pz));
	    I[68] += C[0]*W[a]*(Px*(Cz*zij + Pz));
	    I[69] += C[1]*W[a]*(f8*(Cz*zij + Pz));
	    I[70] += C[1]*W[a]*(Cz*Ix*f8);
	    I[71] += C[0]*W[a]*(Cz*Ix*Py);
	    I[72] += C[1]*W[a]*(Cz*Dy*Ix*f32);
	    I[73] += C[0]*W[a]*(Cz*Ix*f32);
	    I[74] += C[1]*W[a]*(Ix*Pz*Qy);
	    I[75] += C[0]*W[a]*(Cy*Ix*Pz);
	    I[76] += C[1]*W[a]*(Cy*Dz*Ix*f33);
	    I[77] += C[0]*W[a]*(Cy*Ix*f33);
	    I[78] += C[1]*W[a]*(Cy*f21);
	    I[79] += C[1]*W[a]*(Cz*f21);
	    I[80] += C[1]*W[a]*(Px*f4);
	    I[81] += C[1]*W[a]*(Pz*f4);
	    I[82] += C[1]*W[a]*(Dz*f18);
	    I[83] += C[1]*W[a]*(Ix*f10);
	    I[84] += C[1]*W[a]*(Ix*f7);
	    I[85] += C[1]*W[a]*(Iz*f7);
	    I[86] += C[1]*W[a]*(Iz*f28);
	    I[87] += C[1]*W[a]*(Iy*f28);
	    I[88] += C[1]*W[a]*(Iy*f10);
	    I[89] += C[1]*W[a]*(f17*f5);
	    I[90] += C[1]*W[a]*(Dx*Pz*f17);
	    I[91] += C[0]*W[a]*(Pz*f17);
	    I[92] += C[1]*W[a]*(Dz*Px*f17);
	    I[93] += C[0]*W[a]*(Px*f17);
	    I[94] += C[1]*W[a]*(f17*f2);
	    I[95] += C[1]*W[a]*(Qy*f34);
	    I[96] += C[1]*W[a]*(Qy*f19);
	    I[97] += C[1]*W[a]*(Qz*f19);
	    I[98] += C[1]*W[a]*(Qz*f9);
	    I[99] += C[1]*W[a]*(Qx*f9);
	    I[100] += C[1]*W[a]*(Qx*f34);
	    I[101] += C[1]*W[a]*(Cx*Dy*f34);
	    I[102] += C[0]*W[a]*(Cx*f34);
	    I[103] += C[1]*W[a]*(Cy*Dx*f34);
	    I[104] += C[0]*W[a]*(Cy*f34);
	    I[105] += C[1]*W[a]*(Cy*Dz*f19);
	    I[106] += C[0]*W[a]*(Cy*f19);
	    I[107] += C[1]*W[a]*(Cz*Dy*f19);
	    I[108] += C[0]*W[a]*(Cz*f19);
	    I[109] += C[1]*W[a]*(Cz*Dx*f9);
	    I[110] += C[0]*W[a]*(Cz*f9);
	    I[111] += C[1]*W[a]*(Cx*Dz*f9);
	    I[112] += C[0]*W[a]*(Cx*f9);
	    I[113] += C[1]*W[a]*(Cx*f3);
	    I[114] += C[1]*W[a]*(Cy*f3);
	    I[115] += C[1]*W[a]*(Dy*f16);
	    I[116] += C[1]*W[a]*(Dx*f16);
	    I[117] += C[1]*W[a]*(Dx*f18);
	    I[118] += C[0]*W[a]*(f18);
	    I[119] += C[0]*W[a]*(f16);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[120]) {
	double T[120];
	for (int i = 0; i < 120; ++i) {
	    T[i] = I[i];
	}
	I[119] = T[0];
	I[76] = T[1];
	I[75] = T[2];
	I[71] = T[3];
	I[112] = T[4];
	I[114] = T[5];
	I[116] = T[6];
	I[32] = T[7];
	I[36] = T[8];
	I[31] = T[9];
	I[38] = T[10];
	I[59] = T[11];
	I[29] = T[12];
	I[89] = T[13];
	I[69] = T[14];
	I[99] = T[15];
	I[9] = T[16];
	I[39] = T[17];
	I[37] = T[18];
	I[35] = T[19];
	I[30] = T[20];
	I[90] = T[21];
	I[60] = T[22];
	I[0] = T[23];
	I[98] = T[24];
	I[96] = T[25];
	I[109] = T[26];
	I[104] = T[27];
	I[83] = T[28];
	I[55] = T[29];
	I[53] = T[30];
	I[111] = T[31];
	I[113] = T[32];
	I[23] = T[33];
	I[51] = T[34];
	I[21] = T[35];
	I[85] = T[36];
	I[79] = T[37];
	I[107] = T[38];
	I[97] = T[39];
	I[65] = T[40];
	I[95] = T[41];
	I[5] = T[42];
	I[67] = T[43];
	I[7] = T[44];
	I[77] = T[45];
	I[70] = T[46];
	I[110] = T[47];
	I[115] = T[48];
	I[25] = T[49];
	I[80] = T[50];
	I[20] = T[51];
	I[100] = T[52];
	I[10] = T[53];
	I[17] = T[54];
	I[47] = T[55];
	I[49] = T[56];
	I[19] = T[57];
	I[72] = T[58];
	I[74] = T[59];
	I[14] = T[60];
	I[42] = T[61];
	I[12] = T[62];
	I[44] = T[63];
	I[54] = T[64];
	I[56] = T[65];
	I[26] = T[66];
	I[84] = T[67];
	I[24] = T[68];
	I[86] = T[69];
	I[66] = T[70];
	I[6] = T[71];
	I[62] = T[72];
	I[2] = T[73];
	I[68] = T[74];
	I[8] = T[75];
	I[91] = T[76];
	I[1] = T[77];
	I[33] = T[78];
	I[34] = T[79];
	I[73] = T[80];
	I[78] = T[81];
	I[101] = T[82];
	I[92] = T[83];
	I[61] = T[84];
	I[81] = T[85];
	I[50] = T[86];
	I[40] = T[87];
	I[102] = T[88];
	I[43] = T[89];
	I[48] = T[90];
	I[18] = T[91];
	I[103] = T[92];
	I[13] = T[93];
	I[108] = T[94];
	I[88] = T[95];
	I[63] = T[96];
	I[94] = T[97];
	I[106] = T[98];
	I[45] = T[99];
	I[57] = T[100];
	I[87] = T[101];
	I[27] = T[102];
	I[58] = T[103];
	I[28] = T[104];
	I[93] = T[105];
	I[3] = T[106];
	I[64] = T[107];
	I[4] = T[108];
	I[46] = T[109];
	I[16] = T[110];
	I[105] = T[111];
	I[15] = T[112];
	I[117] = T[113];
	I[118] = T[114];
	I[82] = T[115];
	I[52] = T[116];
	I[41] = T[117];
	I[11] = T[118];
	I[22] = T[119];
    }

};

template<>
struct impl<meta::braket<rysq::SP, rysq::P, rysq::P, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[2], double (&I)[36]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 	double vB00[N] __attribute__ ((aligned(16)));
// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 
// 	double vDx[N] __attribute__ ((aligned(16)));
// 	double vDy[N] __attribute__ ((aligned(16)));
// 	double vDz[N] __attribute__ ((aligned(16)));
// 

// 	for (int a = 0; a < N; ++a) {
// 
// 	    vB00[a] = 0.5*t2[a];
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	    vDx[a] = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
// 	    vDy[a] = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
// 	    vDz[a] = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);
// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B00 = 0.5*t2[a];
	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);

	    double Dx = recurrence::coefficient<0>(rBk, -A, rAB, t2[a]);
	    double Dy = recurrence::coefficient<1>(rBk, -A, rAB, t2[a]);
	    double Dz = recurrence::coefficient<2>(rBk, -A, rAB, t2[a]);

	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Qx = (Cx*Dx + B00);
	    double Qy = (Cy*Dy + B00);
	    double Qz = (Cz*Dz + B00);
	    double f1 = (B00 + Dz*Iz);
	    double f11 = (B10 + Cz*Iz);
	    double f12 = (Cy*Iy + B10);
	    double f6 = (Cx*Ix + B10);

	    I[0] += C[1]*W[a]*((B00*(yij + 2*Cy) + Dy*(Cy*Iy + B10)));
	    I[1] += C[1]*W[a]*((Dz*(B10 + Cz*Iz) + B00*(2*Cz + zij)));
	    I[2] += C[1]*W[a]*((Dx*(Cx*Ix + B10) + B00*(xij + 2*Cx)));
	    I[3] += C[1]*W[a]*(Cy*(Dx*xij + Qx));
	    I[4] += C[1]*W[a]*(Cz*(Dx*xij + Qx));
	    I[5] += C[0]*W[a]*((Dx*xij + Qx));
	    I[6] += C[1]*W[a]*(Cz*(Dy*yij + Qy));
	    I[7] += C[1]*W[a]*(Cx*(Dy*yij + Qy));
	    I[8] += C[0]*W[a]*((Dy*yij + Qy));
	    I[9] += C[1]*W[a]*(Cy*Dz*Ix);
	    I[10] += C[1]*W[a]*(Ix*Qy);
	    I[11] += C[1]*W[a]*(Cx*Dy*Iz);
	    I[12] += C[0]*W[a]*(Dy*Iz);
	    I[13] += C[1]*W[a]*(Iz*Qy);
	    I[14] += C[1]*W[a]*(Dx*f11);
	    I[15] += C[1]*W[a]*(Dy*f11);
	    I[16] += C[1]*W[a]*(Dz*f6);
	    I[17] += C[1]*W[a]*(Dy*f6);
	    I[18] += C[1]*W[a]*(Cz*Dy*Ix);
	    I[19] += C[0]*W[a]*(Dy*Ix);
	    I[20] += C[1]*W[a]*(Ix*Qz);
	    I[21] += C[0]*W[a]*(Dz*Ix);
	    I[22] += C[1]*W[a]*(Cx*Dz*Iy);
	    I[23] += C[0]*W[a]*(Dz*Iy);
	    I[24] += C[1]*W[a]*(Dz*f12);
	    I[25] += C[1]*W[a]*(Dx*f12);
	    I[26] += C[1]*W[a]*(Cy*Dx*Iz);
	    I[27] += C[0]*W[a]*(Dx*Iz);
	    I[28] += C[1]*W[a]*(Iz*Qx);
	    I[29] += C[1]*W[a]*(Iy*Qx);
	    I[30] += C[1]*W[a]*(Cz*Dx*Iy);
	    I[31] += C[0]*W[a]*(Dx*Iy);
	    I[32] += C[1]*W[a]*(Iy*Qz);
	    I[33] += C[1]*W[a]*(Cx*f1);
	    I[34] += C[1]*W[a]*(Cy*f1);
	    I[35] += C[0]*W[a]*(f1);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[36]) {
	double T[36];
	for (int i = 0; i < 36; ++i) {
	    T[i] = I[i];
	}
	I[18] = T[0];
	I[35] = T[1];
	I[1] = T[2];
	I[2] = T[3];
	I[3] = T[4];
	I[0] = T[5];
	I[19] = T[6];
	I[17] = T[7];
	I[16] = T[8];
	I[26] = T[9];
	I[14] = T[10];
	I[21] = T[11];
	I[20] = T[12];
	I[22] = T[13];
	I[11] = T[14];
	I[23] = T[15];
	I[25] = T[16];
	I[13] = T[17];
	I[15] = T[18];
	I[12] = T[19];
	I[27] = T[20];
	I[24] = T[21];
	I[29] = T[22];
	I[28] = T[23];
	I[30] = T[24];
	I[6] = T[25];
	I[10] = T[26];
	I[8] = T[27];
	I[9] = T[28];
	I[5] = T[29];
	I[7] = T[30];
	I[4] = T[31];
	I[31] = T[32];
	I[33] = T[33];
	I[34] = T[34];
	I[32] = T[35];
    }

};

template<>
struct impl<meta::braket<rysq::D, rysq::P, rysq::S, rysq::S> > {
    typedef void enable;
    static const bool value = true; 

    static inline void eval(double A, double B,
			    const Vector<3> &rAi, const Vector<3> &rBk,
			    const Vector<3> &rAB,
			    const Vector<3> &rij, const Vector<3> &rkl,
			    const Vector<2> &t2, const Vector<2> &W,
			    const double (&C)[1], double (&I)[18]) {

#define pow(x,y) util::pow<(y)>::eval((x))

	const int N = 2;

	double xij = rij[0];
	double yij = rij[1];
	double zij = rij[2];


// 
// 
// 	double vB10[N] __attribute__ ((aligned(16)));
// 
// 

// 
// 	double vCx[N] __attribute__ ((aligned(16)));
// 	double vCy[N] __attribute__ ((aligned(16)));
// 	double vCz[N] __attribute__ ((aligned(16)));
// 

// 

// 	for (int a = 0; a < N; ++a) {
// 
// 
// 	    vB10[a] = recurrence::coefficient(1.0/A, B, t2[a]);
// 
// 

// 
// 	    vCx[a] = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
// 	    vCy[a] = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
// 	    vCz[a] = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);
// 

// 
// 	}

// #define B00 vB00[a]
// #define B10 vB10[a]
// #define B01 vB01[a]
// #define Cx vCx[a]
// #define Cy vCy[a]
// #define Cz vCz[a]
// #define Dx vDx[a]
// #define Dy vDy[a]
// #define Dz vDz[a]

#if defined (__INTEL_COMPILER) 
#pragma ivdep
#pragma vector aligned
#endif
	for (int a = 0; a < N; ++a) {

	    double B10 = recurrence::coefficient(1.0/A, B, t2[a]);

	    double Cx = recurrence::coefficient<0>(rAi, B, rAB, t2[a]);
	    double Cy = recurrence::coefficient<1>(rAi, B, rAB, t2[a]);
	    double Cz = recurrence::coefficient<2>(rAi, B, rAB, t2[a]);


	    double Ix = (Cx + xij);
	    double Iy = (Cy + yij);
	    double Iz = (Cz + zij);
	    double Px = (B10 + pow(Cx,2));
	    double Py = (B10 + pow(Cy,2));
	    double Pz = (B10 + pow(Cz,2));
	    double f2 = (B10 + Cz*Iz);

	    I[0] += C[0]*W[a]*(Cx*(Py + Cy*yij));
	    I[1] += C[0]*W[a]*(Cz*(Py + Cy*yij));
	    I[2] += C[0]*W[a]*((B10*(3*Cy + yij) + Iy*pow(Cy,2)));
	    I[3] += C[0]*W[a]*(Cx*Cz*Iy);
	    I[4] += C[0]*W[a]*((Iz*pow(Cz,2) + B10*(3*Cz + zij)));
	    I[5] += C[0]*W[a]*((B10*(3*Cx + xij) + Ix*pow(Cx,2)));
	    I[6] += C[0]*W[a]*(Cz*(Px + Cx*xij));
	    I[7] += C[0]*W[a]*(Cy*(Px + Cx*xij));
	    I[8] += C[0]*W[a]*(Cx*Cy*Iz);
	    I[9] += C[0]*W[a]*(Cy*Cz*Ix);
	    I[10] += C[0]*W[a]*(Ix*Py);
	    I[11] += C[0]*W[a]*(Ix*Pz);
	    I[12] += C[0]*W[a]*(Iy*Pz);
	    I[13] += C[0]*W[a]*(Iy*Px);
	    I[14] += C[0]*W[a]*(Iz*Px);
	    I[15] += C[0]*W[a]*(Iz*Py);
	    I[16] += C[0]*W[a]*(Cx*f2);
	    I[17] += C[0]*W[a]*(Cy*f2);
	}

#undef B00
#undef B10
#undef B01
#undef Cx
#undef Cy
#undef Cz
#undef Dx
#undef Dy
#undef Dz

#undef pow
    }

    static inline void reorder(double (&I)[18]) {
	double T[18];
	for (int i = 0; i < 18; ++i) {
	    T[i] = I[i];
	}
	I[9] = T[0];
	I[11] = T[1];
	I[7] = T[2];
	I[10] = T[3];
	I[14] = T[4];
	I[0] = T[5];
	I[4] = T[6];
	I[3] = T[7];
	I[15] = T[8];
	I[5] = T[9];
	I[1] = T[10];
	I[2] = T[11];
	I[8] = T[12];
	I[6] = T[13];
	I[12] = T[14];
	I[13] = T[15];
	I[16] = T[16];
	I[17] = T[17];
    }

};

END_NAMESPACE(rysq, kernel, quadrature)

#endif /* _RYSQ_KERNEL_QUADRATURE2_IMPL_HPP_ */

