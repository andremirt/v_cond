#ifndef _CURYSQ_JK_H
#define _CURYSQ_JK_H

#include "rysq.h"

/** 
    @file 
    @brief function prototypes for 2, 3, 4  to thread mappings
*/

/** 
    @brief Evaluate the <AB|CD> ERI with <AB| mapped explicitly to cuda threads
    @param flags Flags
    @param tol Tolerance
    @param a A shell
    @param b B shell
    @param ni Number of A centers
    @param nj Number of B centers
    @param dev_rij Packed <AB| center pairs in device memory
    @param c C shell
    @param d D shell
    @param nk Number of C centers
    @param nl Number of D centers
    @param dev_rkl Packed |CD> center pairs in device memory
    @param dev_Ix 2-D x integrals in device memory
    @param dev_Iy 2-D y integrals in device memory
    @param dev_Iz 2-D z integrals in device memory
    @param Dij ij density in device memory
    @param Dkl kl density in device memory
    @param Dil il density in device memory
    @param Djk jk density in device memory
    @param alpha Constant scale factor
    @param[out] Fij ij Fock matrix in device memory
    @param[out] Fkl kl Fock matrix in device memory
    @param[out] Fil il Fock matrix in device memory
    @param[out] Fjk jk Fock matrix in device memory
*/
int cuRysq_jk2(int flags, double tol, 
	       Rysq_shell_t a, Rysq_shell_t b, int ni, int nj, double *dev_rij, 
	       Rysq_shell_t c, Rysq_shell_t d, int nk, int nl, double* dev_rkl,
	       double *dev_Ix, double *dev_Iy, double *dev_Iz,
	       double *Dij, double *Dkl, double *Dil, double *Djk,
	       double alpha, double *Fij, double *Fkl, double *Fil, double *Fjk);

/** 
    @brief Evaluate the <AB|CD> ERI with <AB|C> mapped explicitly to cuda threads
    @param flags Flags
    @param tol Tolerance
    @param a A shell
    @param b B shell
    @param ni Number of A centers
    @param nj Number of B centers
    @param dev_rij Packed <AB| center pairs in device memory
    @param c C shell
    @param d D shell
    @param nk Number of C centers
    @param nl Number of D centers
    @param dev_rkl Packed |CD> center pairs in device memory
    @param dev_Ix 2-D x integrals in device memory
    @param dev_Iy 2-D y integrals in device memory
    @param dev_Iz 2-D z integrals in device memory
    @param Dij ij density in device memory
    @param Dkl kl density in device memory
    @param Dil il density in device memory
    @param Djk jk density in device memory
    @param alpha Constant scale factor
    @param[out] Fij ij Fock matrix in device memory
    @param[out] Fkl kl Fock matrix in device memory
    @param[out] Fil il Fock matrix in device memory
    @param[out] Fjk jk Fock matrix in device memory
*/
int cuRysq_jk3(int flags, double tol, 
	       Rysq_shell_t a, Rysq_shell_t b, int ni, int nj, double *dev_rij, 
	       Rysq_shell_t c, Rysq_shell_t d, int nk, int nl, double* dev_rkl,
	       double *dev_Ix, double *dev_Iy, double *dev_Iz,
	       double *Dij, double *Dkl, double *Dil, double *Djk,
	       double alpha, double *Fij, double *Fkl, double *Fil, double *Fjk);

/** 
    @brief Evaluate the <AB|CD> ERI with <AB|CD> mapped explicitly to cuda threads
    @param flags Flags
    @param tol Tolerance
    @param a A shell
    @param b B shell
    @param ni Number of A centers
    @param nj Number of B centers
    @param dev_rij Packed <AB| center pairs in device memory
    @param c C shell
    @param d D shell
    @param nk Number of C centers
    @param nl Number of D centers
    @param dev_rkl Packed |CD> center pairs in device memory
    @param dev_Ix 2-D x integrals in device memory
    @param dev_Iy 2-D y integrals in device memory
    @param dev_Iz 2-D z integrals in device memory
    @param Dij ij density in device memory
    @param Dkl kl density in device memory
    @param Dil il density in device memory
    @param Djk jk density in device memory
    @param alpha Constant scale factor
    @param[out] Fij ij Fock matrix in device memory
    @param[out] Fkl kl Fock matrix in device memory
    @param[out] Fil il Fock matrix in device memory
    @param[out] Fjk jk Fock matrix in device memory
*/
int cuRysq_jk4(int flags, double tol, 
	       Rysq_shell_t a, Rysq_shell_t b, int ni, int nj, double *dev_rij, 
	       Rysq_shell_t c, Rysq_shell_t d, int nk, int nl, double* dev_rkl,
	       double *dev_Ix, double *dev_Iy, double *dev_Iz,
	       double *Dij, double *Dkl, double *Dil, double *Djk,
	       double alpha, double *Fij, double *Fkl, double *Fil, double *Fjk);

#endif // _CURYSQ_JK_H
