#ifndef _DISTRIBUTED_OPS_HPP_
#define _DISTRIBUTED_OPS_HPP_

//#include "matrix/meta.hpp"

namespace distributed {


    template<class Array, class E>
    void assign(Array &D, const matrix::meta_matrix<E> &e) {
	typedef typename Array::index_array index_array;
	index_array lo(0), hi = D.dims() - 1;
	index_array ld = hi - lo;
	typedef typename matrix::meta_matrix<E>::index_iterator::indexed indexed;
	foreach (indexed &ab, e.matrix_index()) {
	    unpack((int a, b, c), ab);
	    //unpack<indexed>(a,b) = ab;
	    typedef typename E::value_type* pointer;
	    const E &Mab = e.m(i);
	    pointer src = const_cast<pointer>(Mab.data().begin());
	    D.put(lo, hi, src, ld);
	}
    }

    template<class Array, class E>
    void assign(Array &D, const E &e) {
    }

}

#endif /* _DISTRIBUTED_OPS_HPP_ */
