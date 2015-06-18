#ifndef _RYSQ_CUDA_MATRIX_HPP_
#define _RYSQ_CUDA_MATRIX_HPP_

#include <boost/smart_ptr.hpp>

#include "rysq/cuda-memory.hpp"

namespace rysq {
    namespace cuda {

	// struct zero_matrix {
	//     zero_matrix(size_t size1, size_t  size2)
	// 	: size1_(size1), size2_(size2) {}
	//     size_t size1() const { return size1_; }
	//     size_t size2() const { return size2_; }
	// private:
	//     size_t size1_, size2_;
	// };

	template<typename T, class A_ = array<T> >
	struct block_matrix :  block_matrix_base <T> {
	    typedef A_ array_type;
	    typedef block_matrix_base<T> base;
	    typedef rysq::block_matrix<T,boost::shared_array<T> > host_matrix;
	    using base::size1;
	    using base::size2;
	    using base::size;
	    using base::block1;
	    using base::block2;
	    block_matrix() {}
	    T* block(int a,int b) { return data() +  base::layout_.block_at(a,b); }
	    const T* block(int a,int b) const {
		return data() + base::layout_.block_at(a,b);
	    }
	    void clear(){//(const zero_matrix &A) {
		// base::check_size(A);
		memset(device_ptr<void>(data_.data()),  char(0), size()*sizeof(T));
	    }
	    template<class M>
	    void assign(const M &A) {
		base::check_size(A);
		host_matrix B(size1(), size2(), block1(), block2());
		B.assign(A);
		memcpy<T>(B.data(), data_, B.size());
	    }
	    host_matrix copy() const {
		host_matrix B(size1(), size2(), block1(), block2());
		memcpy<T>(data_, B.data(), B.size());
		return B;
	    }
	    void resize(size_t size1, size_t size2, size_t block1, size_t block2) {

		base::layout_ = block_matrix_layout(size1, size2, block1, block2);
		base::size1_ = size1;  base::size2_ = size2;
		data_.resize(size());
		// std:: cout << this->size1() <<this-> size2() <<this-> block1() <<this-> block2() << std::endl;
	    }
	    T* data() { return data_.data(); }
	    const T* data() const { return data_.data(); }
	private:
	    array_type data_;
	};

}
}


#endif /* __RYSQ_CUDA_MATRIX_HPP_ */
