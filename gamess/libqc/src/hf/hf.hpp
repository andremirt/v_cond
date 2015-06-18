#ifndef QC_HF_H
#define QC_HF_H

/**
 * @file 
 * @brief  HF functions
 */

#include "core/basis.hpp"
#include "matrix/matrix.hpp"

namespace hf {


    template<class M1, class M2>
    void Dmax(const Basis &basis, const M1 &D, M2 &Dmax);

    /** 
     * Construct Fock JK matrix operator
     * 
     * @param basis basis set
     * @param D density matrix
     * @param F fock matrix
     * @param tolerance integral tolerance threshold
     * @param screen integral screening matrix
     */
    template<class M>
    void fock(const Basis &basis,
	      const matrix::meta_matrix<BlockMatrix> &D,
	      matrix::meta_matrix<BlockMatrix> &F,
	      const M *Kmax = NULL, const M *Dmax = NULL,
	      double tolerance=1.0e-10);

}

#include "hf/screen.hpp"

#endif // QC_HF_H
