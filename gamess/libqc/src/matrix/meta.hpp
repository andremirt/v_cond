#ifndef _MATRIX_META_MATRIX_H_
#define _MATRIX_META_MATRIX_H_


#include <vector>
#include <algorithm>
#include <numeric>
#include <deque>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/bind.hpp>

#include "matrix/matrix.hpp"
#include "matrix/matrix_base.hpp"
#include "matrix/slice.hpp"

#include "iterator/iterator.hpp"
#include "util/util.hpp"
#include "util/array.hpp"
#include "util/range.hpp"


    template<typename T, size_t N, typename U>
    boost::array<T,N> make_array(const U *in) {
	boost::array<T,N> out;
	std::copy(in, in+N, out.begin());
	return out;
    }

    template<typename T, typename U, size_t N>
    boost::array<T,N> make_array(const boost::array<U,N> &in) {
	return make_array<T,N>(in.elems);
    }

namespace matrix {

    using namespace ::iterator;

    template<class M>
    struct submatrix_view {
	static const size_t dimensionality = 2;
	typedef boost::fortran_storage_order storage_order_type;
	typedef typename M::matrix reference;
	typedef typename M::const_matrix const_reference;
	typedef typename M::size_type size_type;
	typedef array<size_type,2> size_array;
	submatrix_view(M &A) : A_(A) {}
	size_array shape() const {
	    return make_array<size_type,2>(A_.submatrix_array_.shape());
	}
	size_type num_elements() const { return A_.submatrix_array_.num_elements(); }
	reference operator()(array<size_type,2>) { return A_.m(index); }
	const_reference operator()(array<size_type,2>) const {
	    return (const_cast<submatrix_view&>(*this))(index);
	}
    private: M &A_;
    };


    template<typename T>
    size_t element_bin(const std::vector<T> &bins, T e) {
	typename std::vector<T>::const_iterator it;
	it = std::upper_bound(bins.begin(), bins.end(), e) - 1;
	assert(it < bins.end()-1);
	//std::cout << size_t(it - bins.begin()) << "\n";	    
	return size_t(it - bins.begin());
    }


    template<class M>
    struct meta_matrix_base : matrix_semantics<typename M::value_type>
    {
    public:
	typedef meta_matrix_base self;
	typedef typename M::value_type value_type;

	typedef int index_type;
	typedef boost::array<index_type,2> index_array;
	typedef size_t size_type;
	typedef size_type shape_type[2];

	typedef M matrix;
	typedef const M const_matrix;

	typedef boost::multi_array<matrix*,2> matrix_array_type;
	typedef indexed_iterator<2,submatrix_view<self> > matrix_iterator;

	meta_matrix_base(size_type size1, size_type size2)
	    : submatrix_array_(boost::extents[size1][size2]),
	      submatrix_view_(*this)
	{ }

	virtual ~meta_matrix_base() {
	    for (RANGE(it, matrix_range())) {
		free_(make_array<index_type>(it.index()));
	    }
	}

	const matrix_array_type& array() const {return submatrix_array_; }

	iterator_range<matrix_iterator> matrix_range() {
	    return indexed_range(submatrix_view_);
	}

	virtual matrix& m(int a, int b) { return *matrix_(a,b); }

	virtual const_matrix& m(int a, int b) const { return *matrix_(a,b); }

// 	template<class Index>
// 	matrix& m(const Index &ab) { return m(ab[0], ab[1]); }

// 	template<class Index>
// 	const_matrix& m(const Index &ab) const { return m(ab[0], ab[1]); }

    protected:
	friend class submatrix_view<self>;

	mutable matrix_array_type submatrix_array_;
	submatrix_view<meta_matrix_base> submatrix_view_;

	void free_(const boost::array<index_type,2> &ab) {
	    matrix* &m = submatrix_array_(ab);
	    if (m)  { delete m; m = NULL; }
	}

	matrix* matrix_(int a, int b) const { return submatrix_array_[a][b]; }

	matrix& matrix_(index_type &i, index_type &j,
			const std::vector<index_type> (&index)[2]) {
	    int a = element_bin(index[0], i);
	    int b = element_bin(index[1], j);
	    i -= index[0].at(a);
	    j -= index[1].at(b);
	    //printf("%i %i\n" ,i,j);
	    return m(a,b);
	}	    

	const_matrix& matrix_(index_type &i, index_type &j,
			      const std::vector<index_type> (&index)[2]) const {
	    return const_cast<self*>(this)->matrix_(i, j, index);
	}	    
    };

    template<class M>
    class meta_matrix : public meta_matrix_base<M> {
    public:
	typedef meta_matrix self;
	typedef meta_matrix_base<M> base;

	typedef typename base::reference reference;
	typedef typename base::const_reference const_reference;

	typedef typename base::matrix matrix;
	typedef typename base::const_matrix const_matrix;



	typedef util::Array<index_type,2> index_array;
	typedef util::Array<size_type,2> size_array;

	typedef std::vector<index_type> index_vector;
	typedef std::vector<size_type> size_vector;


	meta_matrix(const size_vector (&matrix_dim)[2], bool allocate = true)
	    : base(matrix_dim[0].size(), matrix_dim[1].size())
	{
	    initialize(matrix_dim, allocate);
	}

	meta_matrix(const size_vector &matrix_dim, bool allocate = true)
	    : base(matrix_dim.size(), matrix_dim.size())
	{
	    size_vector dims[] = { matrix_dim, matrix_dim };
	    initialize(dims, allocate);
	}


    protected:

	std::vector<int> elems_index[2];
	size_type max_matrix_size_;

	void allocate_(const index_array ab) const {
	    size_array size = matrix_size(ab[0], ab[1]);
	    base::submatrix_array_(ab) = new matrix(size[0], size[1]);
	}


	    

	void initialize(const size_vector (&matrix_dim)[2], bool allocate) {
	    size_array extents( matrix_dim[0].size(), matrix_dim[1].size() );

 	    set_index(matrix_dim[0], elems_index[0]);
 	    set_index(matrix_dim[1], elems_index[1]);

	    max_matrix_size_ = 1;
	    foreach (const size_vector &d, matrix_dim) {
		max_matrix_size_ *= std::max_element(d.begin(), d.end())[0];
	    }

	    for ( RANGE( it, this->matrix_range())){
		if (allocate) allocate_(make_array<index_type>(it.index()));
	    }

	}




	static void set_index(const size_vector &dims,
			      index_vector &index) {
	    index.push_back(0);
	    foreach (size_type size, dims) {
		index.push_back(index.back() + size);
		//std::cout << num_blocks << " " <<  block_index.back() << "\n";
	    }
	}




	typedef matrix* matrix_pointer;

	//const matrix_pointer& matrix_(int a, int b) const { return matrix_array_[a][b]; }

    public:
	size_t size1() const { return elems_index[0].back(); }
	size_t size2() const { return elems_index[1].back(); }

	const index_array matrix_index(index_type a, index_type b) const {
	    return index_array(elems_index[0].at(a), elems_index[1].at(b));
	}

	const size_array matrix_size(int a, int b) const {
	    size_type m = elems_index[0].at(a+1) - elems_index[0].at(a);
	    size_type n = elems_index[1].at(b+1) - elems_index[1].at(b);
	    return size_array(m, n);
	}

	reference operator()(int i, int j) {
 	    matrix &m = base::matrix_(i, j, elems_index);
 	    return m(i,j);
	}

	const_reference operator()(int i, int j) const {
	    return const_cast<self*>(this)->operator()(i, j);
	}

	template<class E>
	void operator=(const E &e) {
	    assign(*this, e);
	}

    };


    template<class M>
    class block_meta_matrix : public meta_matrix<M>,
			      public block_matrix< typename M::value_type  >{

    public:
	typedef block_meta_matrix self;
	typedef meta_matrix<M> base;
	typedef block_matrix<typename M::value_type> block_base;

	typedef typename base::value_type value_type;
	//typedef typename  base::element element;
	typedef typename base::reference reference;
	typedef typename base::const_reference const_reference;
	//typedef typename base::pointer pointer;

	typedef typename base::size_vector size_vector;
	typedef typename base::index_vector index_vector;
	typedef typename base::matrix_iterator matrix_iterator;

	reference operator()(index_type i, index_type j) {
	    return base::operator()(i,j);
	}

	const_reference operator()(index_type i, index_type j) const {
	    return base::operator()(i,j);
	}

	size_type size1() const { return base::size1(); }
	size_type size2() const { return base::size2(); }

	block_meta_matrix(const size_vector (&matrix_dim)[2],
			  const size_vector (&block_dim)[2],
			  bool allocate = true)
	    : base(matrix_dim, allocate)
	{
	    for (RANGE(it, base::matrix_range())) {
		const index_type a = it.index()[0], b = it.index()[1];
		//printf("%i, %i\n", a,b);
		shape_type shape = { block_dim[0].at(a), block_dim[1].at(b) };
		base::matrix_(a,b)->set_block(shape);
	    }
 	    set_block_index(matrix_dim[0],  block_dim[0], block_index[0]);
 	    set_block_index(matrix_dim[1],  block_dim[1], block_index[1]);
	    block_base::set(this->shape());
	}


	static void set_block_index(const size_vector &matrix_dims,
				    const size_vector &block_dims,
				    index_vector &block_index) {
	    block_index.push_back(0);
	    for (uint i = 0; i < block_dims.size(); ++i) {
		size_type num_blocks =  matrix_dims[i]/block_dims[i];
		block_index.push_back(block_index.back() + num_blocks);
		//std::cout << num_blocks << " " <<  block_index.back() << "\n";
	    }
	}

	value_type* b(int i, int j) {
	    M &A = matrix_(i, j, block_index);
	    return A.b(i,j);
	}

	const value_type* b(int i, int j) const {
	    return const_cast<self*>(this)->b(i, j);
	}

	template<class E>
	void operator=(const E &e) { base::operator=(e); }

    protected:
 	index_vector block_index[2];

    };

//     template<class M>
//     struct data {
// 	typedef const M const_matrix;
// 	typedef M::array_type array_type;
// 	data(const_matrix &A) : A_(A) {}
//     private:
// 	const_matrix &A_;
//     };

    template<class M, class Adaptable>
    class meta_matrix_decorator : public meta_matrix<M> {
    public:
	typedef meta_matrix<M> base;
	typedef typename base::value_type value_type;
	//typedef typename base::pointer pointer;
	typedef typename base::index_array index_array;
	typedef typename base::size_vector size_vector;

	typedef M matrix;
	typedef const M const_matrix;

	meta_matrix_decorator(const size_vector &matrix_dim,
			      Adaptable& A, size_t cache_size = 1)
	    : base(matrix_dim, false), A(A)
	{
	}

	meta_matrix_decorator(const size_vector (&matrix_dim)[2],
			      Adaptable& A, size_t cache_size = 1)
	    : base(matrix_dim, false), A(A)
	{
	}

	~meta_matrix_decorator() { }

    protected:

	Adaptable& A;

	
	typedef std::deque<index_array> fifo;
	typedef typename fifo::iterator fifo_iterator;
	mutable fifo fifo_;

	void push_matrix(index_type a, index_type b) const {
	    index_array index(a,b);
	    fifo_iterator it = std::find(fifo_.begin(), fifo_.end(), index);
	    if (it == fifo_.end()) {
		fifo_.push_back(index);
		base::allocate_(index);

		index_array from = this->matrix_index(a,b);
 		index_array to = from + this->matrix_size(a,b) - 1;
 		assign(*base::matrix_(a,b), slice<>(from, to)(A));

	    }
	    else if (it != fifo_.end() - 1) {
		fifo_.erase(it);
		fifo_.push_back(index);
	    }
	    
	}

	virtual const_matrix& m(int a, int b) const {
	    push_matrix(a,b);
	    return *base::matrix_(a,b);
	}		

	virtual matrix& m(int a, int b) {
	    push_matrix(a,b);
	    return *base::matrix_(a,b);
	}		

    };

}


#endif /* _MATRIX_META_MATRIX_H_ */
