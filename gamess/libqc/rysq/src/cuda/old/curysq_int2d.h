#ifndef _CURYSQ_INT2D_H
#define _CURYSQ_INT2D_H

/**
   @file
 */

#include <rysq.h>


/**
   @brief
   @param a
 */
template <class T> int cuRysq_int2d(Rysq_shell_t a, Rysq_shell_t b, Rysq_shell_t c, Rysq_shell_t d,
				    double *dev_rij, double *dev_rkl, int n, int dev_index4[],
				    T *dev_Ix, T *dev_Iy, T *dev_Iz);







#endif // _CURYSQ_INT2D_H


