#ifndef _CUDA_HOST_ARRAY_HPP_
#define _CUDA_HOST_ARRAY_HPP_

/**
 * @file 
 * @brief  Basic CUDA C++ template types and operations
 */

#include <iostream>
#include <vector>
#include <string>
#include <boost/array.hpp>

#include "cuda/host/core.hpp"

namespace cuda {
    namespace host {

	template<typename T, size_t N>
	struct array_base;

	template<typename T, size_t N>
	struct array_base {
	    size_t size() const { return size_; }
	    void memset(char value) { host::memset(this->data(), N*size_, value); }
	    device_ptr<T> data() { return data_; }
	    device_ptr<const T> data() const { return data_; }
	    device_ptr<T> operator+(size_t offset) { return data_ + offset*N; }
	    device_ptr<const T> operator+(size_t offset) const {
		return data_ + offset*N;
	    }
	    T operator[](size_t i) const { return *((data_ + i*N).data()); }
	protected:
	    array_base() : data_(wrap<T>(NULL)), size_(0) {}
	    array_base(size_t size) { allocate(size); }
	    array_base(device_ptr<T> data, size_t size) : data_(data), size_(size) {}
	    void allocate(size_t size) { 
		data_ = malloc<T>(N*size);
		size_ = size;
	    }
	    device_ptr<T> data_;
	    size_t size_;
	};	

	template<class T>
	void copy(const std::vector<T> &from, array_base<T,1> &to) {
	    if (!from.empty())
		detail::copy(&from.front(), to.data(), from.size()*sizeof(T));
	}

	template<class T, size_t N>
	void copy(const std::vector<boost:: array<T,N> > &from, array_base<T,N> &to) {
	    if (!from.empty())
		detail::copy(&from.front(), to.data(), N*from.size()*sizeof(T));
	}

	template<typename T, size_t N>
	void copy(const array_base<T,N> &from,  const std::string &symbol) {
	    // std::cout << N <<  " " << from.size() << std::endl;
	    copy(from.data(), symbol, from.size()*N);
	}

	template<typename T,size_t N = 1>
	struct array : array_base<T,N> {
	    typedef array_base<T,N> base;
	    array() {}
	    explicit array(size_t size ) : base(size) {}
	    array(const std::vector<boost::array<T,N> > &v) : base(v.size()) {
		copy(v, *this);
	    }
	    ~array() { if (base::data_) free(base::data_); }
	    void reserve(size_t size) {
		if (base::data_) throw;
		base::allocate(size);
	    }
	    array& operator=(char value) { base::fill(value); return *this; }
	private:
	    array(const array&) {}
	    void operator=(const array&) {}
	};

	template<typename T,size_t N = 1>
	struct array_ref : array_base<T,N> {
	    typedef array_base<T,N> base;
	    array_ref() {}
	    array_ref& operator+(size_t offset) const;
	};

	template<typename T,size_t N>
	struct array_ref<const T, N> : array_base<const T,N> {
	    typedef array_base<const T,N> base;
	    array_ref() {}
	    array_ref(device_ptr<const T> data, size_t size) : base(data, size) {}
	    array_ref(const array<T,N> &a) : base(a.data(), a.size()) {}
	    array_ref& operator+(size_t offset) const;
	};


    }
}

#endif /* _CUDA_HOST_ARRAY_HPP_ */

