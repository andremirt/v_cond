#ifndef _CURYSQ_ERI_H
#define _CURYSQ_ERI_H

/** 
    @file 
    @brief function prototypes for 2, 3, 4 ERI to thread mappings
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
    @param alpha Constant scale factor
    @param[out] dev_I Final integral in device memory
*/
int cuRysq_eri2(int flags, double tol, 
		Rysq_shell_t a, Rysq_shell_t b, 
		Rysq_shell_t c, Rysq_shell_t d, int n,
		double* dev_Ix, double* dev_Iy, double *dev_Iz, double scale, 
		double* dev_I);

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
    @param alpha Constant scale factor
    @param[out] dev_I Final integral in device memory
*/
int cuRysq_eri3(int flags, double tol, 
		Rysq_shell_t a, Rysq_shell_t b, 
		Rysq_shell_t c, Rysq_shell_t d, int n,
		double* dev_Ix, double* dev_Iy, double *dev_Iz, double scale, 
		double* dev_I);

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
    @param alpha Constant scale factor
    @param[out] dev_I Final integral in device memory
*/
int cuRysq_eri4(int flags, double tol, 
		Rysq_shell_t a, Rysq_shell_t b, 
		Rysq_shell_t c, Rysq_shell_t d, int n,
		double* dev_Ix, double* dev_Iy, double *dev_Iz, double scale,
		double* dev_I);

#endif // _CURYSQ_ERI_H
