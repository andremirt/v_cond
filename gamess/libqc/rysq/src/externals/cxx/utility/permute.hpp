#ifndef _CXX_UTILITY_PERMUTE_HPP
#define _CXX_UTILITY_PERMUTE_HPP

#include <boost/array.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include "cxx/array.hpp"

#ifndef __host__
#define __host__
#endif

#ifndef __device__
#define __device__
#endif


namespace cxx {
    namespace utility {

	namespace permutation_ {

	    template<class Q>
	    struct value_{
		typedef typename Q::value_type type;
	    };

	    template< typename T>
	    struct value_<T[4]> { typedef T type; };

	    __host__ __device__
	    static int index(int index, int mask) {
		return (mask>>(index*2) & 0x03);
	    }

	}

	template<class Q>
	struct permutation {
	    typedef typename  permutation_::value_<Q>::type value_type;
	    __host__ __device__
	    permutation(const Q &quartet, int mask) : data_(quartet), mask_(mask) {}
	    __host__ __device__
	    const value_type& operator[](size_t i) const {
		return data_[permutation_::index(i, mask_)];
	    }
	private:
	    const Q &data_;
	    const int mask_;
	};

	template<class Q>
	__host__ __device__
	permutation<Q> make_permutation(const Q &quartet, int mask) {
	    return permutation<Q>(quartet, mask);
	}

	template <typename T>
	__host__ __device__
	void permute(T &i, T &j, T &k, T &l, int mask) {
	    T q[4] = {i, j, k, l};
	    permutation<T[4]> p(q, mask);
	    i = p[0];
	    j = p[1];
	    k = p[2];
	    l = p[3];
	}

	template<typename T>
	__host__ __device__
	T permute(int i, const T *array, int mask) {
	    return array[permutation_::index(i, mask)];
	}

	template <typename T>
	__host__ __device__
	boost::array<T,4> permute(const boost::array<T,4> &array, int mask) {
	    permutation<boost::array<T,4> > p(array, mask);
	    boost::array<T,4> q = {{ p[0], p[1], p[2], p[3] }};
	    // std::cout << mask << q << array << std::endl;
	    return q;
	}

	template< class Q>
	__host__ __device__
	 Q permute(const Q &array, int mask) {
	    Q q;
	    permutation<Q> p(array, mask);
	    for (int i = 0; i < 4; ++i) {
		q[i] = p[i];
	    }
	    return q;
	}



    }
}

#endif /* _CXX_UTILITY_PERMUTE_HPP */
