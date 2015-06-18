/**
   @file
   @brief libqc  Fortran bindings
   @details The fortran binding provides lower case subroutines with and without
   trailing underscore for compatibility with different fortran compilers.
*/

#ifdef HAVE_CONFIG_H 
#include "config.h"
#endif

#include "bindings/fortran.h" 
#include "core/molecule.hpp"
#include "core/basis.hpp"
#include "core/normalize.hpp"
#include "hf/hf.hpp"

#include "util/util.hpp"
#include "matrix/matrix.hpp"

namespace ublas = boost::numeric::ublas;

//adapter to use allocated memory for ublas containers
typedef ublas::shallow_array_adaptor<double> ShallowAdapter;

//symmetric matrix using adapter
typedef ublas::symmetric_matrix<double, ublas::upper, ublas::column_major,
				ShallowAdapter> SymmetricAdapter;



extern "C" {

    //     /**
    //        @brief creates new molecule
    //        @return new molecule handle
    //     */
    //     Integer molecule_new (...) {
    // 	return Molecule_new (...);
    //     }

    /**
       @brief deletes molecule
       @param molecule molecule handle
    */
    void molecule_delete(Integer *molecule) {
	Molecule_delete(*molecule);
    }

    //     /**
    //        @brief creates new basis
    //        @param molecule molecule handle
    //        @return new basis handle
    //     */
    //     Integer basis_new (Integer *molecule,...) {
    // 	return Basis_new (*molecule,...);
    //     }

    /**
       @brief deletes basis
       @param basis basis handle
    */
    void basis_delete(Integer *basis) {
	Basis_delete(*basis);
    }

    void basis_sort(Integer *basis) {
        Basis *basisObj = Basis_find_object(*basis);
        basisObj->sort();
    }


    /**
       @brief computes jk fock matrix operator
       @param basis basis handle
       @param D density matrix
       @param F fock matrix
    */
    void hf_fockjk(const Integer *basis_handle, const double *D, double *F,
		   const double *screen, const double *tolerance) {
	const Basis &basis = *Basis_find_object(*basis_handle);

	int N = basis.numShells();
	//int N2 = binomial2 (N);
	int n = basis.numFunctions();
	int n2 = binomial2 (n);

	//copy packed symmetric matrix to square matrix
	//const Matrix tmp = SymmetricAdapter(n, ShallowAdapter(n2, const_cast<double*>(D)));
	//Matrix screen2 = SymmetricAdapter(N, ShallowAdapter(N2, const_cast<double*>(screen)));

	typedef matrix::symmetric_adapter<double> symmetric;
	typedef matrix::const_symmetric_adapter<double> const_symmetric;
	typedef matrix::block_meta_matrix<BlockMatrix> MetaMatrix;
	typedef MetaMatrix::size_vector size_vector;

	size_vector block_sizes, shell_sizes, matrix_sizes;
	foreach (const Basis::Block& block, basis.blocks()) {
	    shell_sizes.push_back(block.shell().size());
	    block_sizes.push_back(block.size());
	    matrix_sizes.push_back(shell_sizes.back()*block_sizes.back());
	}
	const size_vector matrix_dims[] = { matrix_sizes, matrix_sizes };
	const size_vector matrix_block_dims[] = { shell_sizes, shell_sizes };
	MetaMatrix D_(matrix_dims, matrix_block_dims);
	MetaMatrix F_(matrix_dims, matrix_block_dims);

	// matrix::meta_matrix<Matrix> Dmax(block_sizes);//, Kmax(block_sizes);
	Matrix Dmax(N, N);

	// Matrix K_(N,N);
	// typedef matrix::meta_matrix<Matrix> meta_matrix;
	// // matrix::meta_matrix_decorator<Matrix, Matrix> Kmax(block_sizes, K_);
	// matrix::meta_matrix<Matrix> Kmax(block_sizes);
	// Kmax = K_;
	Matrix Kmax(N, N);

	typedef matrix::permutation<int> P;
	matrix::const_symmetric_adapter<double>  screen_adapter(screen,N);
	matrix::assign(Kmax, P(basis.shell_permutations())(screen_adapter));

	P Pf(basis.function_permutations());

	const matrix::const_symmetric_adapter<double> D_adapter(D,n);
	D_ = Pf(D_adapter);
	hf::Dmax(basis, D_, Dmax);
	basis.normalize(D_);

	//matrix::zero(F_);
	F_ = 0;

	//matrix::meta_matrix_decorator<BlockMatrix, MetaMatrix&> d(dims, dims, D_, 6);

	//double cutoff = 1.0e-15;
	double cutoff = *tolerance;
	hf::fock(basis, D_, F_, &Kmax, &Dmax, cutoff);
	//hf::fock(basis, d, F_, &Kmax, &Dmax, cutoff);
	basis.normalize(F_);

	//copy square symmetric matrix to packed upper
	SymmetricAdapter Fsymm(n, ShallowAdapter(n2, F));

	matrix::symmetric_adapter<double> F_adapter(F,n);
	Pf(F_adapter) = F_;

     }

}


extern "C" {

    void basis_delete_(Integer *basis);
    void basis_sort_(Integer *basis);
    void molecule_delete_(Integer *molecule);
    void hf_fockjk_(const Integer *basis, const double *D, double *F,
		    const double *screen, const double *tolerance);

#pragma weak molecule_delete_ = molecule_delete
#pragma weak basis_delete_ = basis_delete
#pragma weak basis_sort_ = basis_sort
#pragma weak hf_fockjk_ = hf_fockjk

}
