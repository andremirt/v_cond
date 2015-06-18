#ifndef _MP2_HPP_
#define _MP2_HPP_

#include <boost/mpl/range_c.hpp>
#include <boost/mpl/back_inserter.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/transform.hpp>


#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/as_vector.hpp>

namespace detail {

    namespace mpl = boost::mpl;
    namespace fusion = boost::fusion;

    template<size_t N>
    struct index {
	struct value_ {
	    typedef const int& type;
	    template<class> struct apply { typedef typename value_::type type; };
	};
	typedef typename value_::type value_type;
	typedef typename mpl::copy<mpl::range_c<int,0,N>,
				   mpl::back_inserter<mpl::vector_c<value_type> >
				   >::type v;
	typedef typename fusion::result_of::as_vector<
	    typename mpl::transform <v, value_>::type>::type type;
    };

    template<size_t N, typename T>
    struct  layout {
	typedef typename index<N>::type index_vector;
	size_t element_at(const index_vector &index) const;
    };

}

template<class T, size_t N, typename value_type>
struct tensor_operator;

template< class T, typename value_type>
struct tensor_operator<T, 4, value_type> {
    value_type& operator()(int i, int j, int k, int l) {
	return tensor_[typename T::index_vector(i, j, k, l)];
    }
    const value_type& operator()(int i, int j, int k, int l) const {
	return tensor_[typename T::index_vector(i, j, k, l)];
    }
protected:
    tensor_operator(T &tensor) : tensor_(tensor) {}
    T &tensor_;
};

template<size_t N, typename T = double>
struct tensor : tensor_operator<tensor<N,T>, N, T> {
    typedef T value_type;
    static const size_t rank = N;
    typedef detail::layout<N,T> layout_type;
    typedef typename layout_type::index_vector index_vector;
    value_type& operator[](const index_vector &index);
    const value_type& operator[](const index_vector &index) const;
private:
    layout_type layout_;
};

struct mp2 {
    tensor<4> T;
    void apply() {
	T(0,0,0,0) = 0;
    }
};

#endif /* _MP2_HPP_ */
