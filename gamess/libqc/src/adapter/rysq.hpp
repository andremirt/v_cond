#ifndef _ADAPTER_RYSQ_HPP_
#define _ADAPTER_RYSQ_HPP_

#include <rysq.hpp>
#include "core/basis.hpp"
#include "matrix/matrix.hpp"

namespace adapter {
    namespace rysq {

	struct Shell : ::rysq::Shell {
	    typedef ::rysq::Shell super;
	    Shell(const Basis::Shell &shell)
		: super(::rysq::type(shell.type()), shell.K(),
			shell.exps(), &shell.coeffs().front()) {
	    }
	};

	template<typename T>
	struct block_matrix : ::rysq::block_matrix_adapter<T> {
	    typedef ::rysq::block_matrix_adapter<T> base;
	    block_matrix() {}
	    template<class Matrix>
	    block_matrix(Matrix &A)
		: base(A.size1(), A.size2(), A.data(), A.block1(), A.block2()) {}
	    };

	struct fock_matrix : block_matrix<double> {
	    typedef block_matrix<double> base;
	    fock_matrix() {}
	    template<class Matrix>
	    fock_matrix(Matrix &A) : base(A) {}
	};

	struct density_matrix : block_matrix<const double> {
	    typedef block_matrix<const double> base;
	    density_matrix() {}
	    template<class Matrix>
	    density_matrix(const Matrix &A) : base(A){}
	};

    }
}

#endif /* _ADAPTER_RYSQ_HPP_ */
