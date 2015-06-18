#ifndef _CUDAXX_CUDAXX_HPP_
#define _CUDAXX_CUDAXX_HPP_

/**
 * @file 
 * @brief  Basic CUDA C++ template types and operations
 */

#include <iostream>
#include <vector>
#include <string>
#include <boost/array.hpp>

namespace cudaxx {

    int initialize(size_t device = 0);
 
    template<typename T>
    struct device_ptr;

    template<typename T>
    device_ptr<T> wrap(T *data);

    template<typename T>
    device_ptr<T> cast(device_ptr<const T> ptr);

    template<typename T, typename U>
    device_ptr<T> reinterpret(device_ptr<U> ptr);

    /**
     * @brief device pointer wrapper
     * @tparam T datatype
     */
    template<typename T>
    struct device_ptr_base {
	T* data() { return data_; }
	const T* data() const { return data_; }
	device_ptr<T> operator+(size_t offset) {
	    return wrap(data_ + offset);
	}
	device_ptr<const T> operator+(size_t offset) const {
	    return wrap(data_ + offset);
	}
	operator bool() const { return data_; }
    protected:
	friend device_ptr<T> wrap<>(T*);
	T* data_;		/**< device pointer */
	device_ptr_base() {}
    };


    /**
     * @brief device pointer wrapper
     * @tparam void datatype
     */
    template<>
    struct device_ptr<const void> : device_ptr_base<const void> {
    };

    /**
     * @brief device pointer wrapper
     * @tparam void datatype
     */
    template<>
    struct device_ptr<void> : device_ptr_base<void> {
	operator device_ptr<const void>() const {
	    return reinterpret<const void>(*this);
	}
    };

    /**
     * @brief device pointer wrapper
     * @tparam T datatype
     */
    template<typename T>
    struct device_ptr<const T> : device_ptr_base<const T> {
	operator device_ptr<const void>() const {
	    return reinterpret<const void>(*this);
	}
    };

    /**
     * @brief device pointer wrapper
     * @tparam T datatype
     */
    template<typename T>
    struct device_ptr : device_ptr_base<T> {
	operator device_ptr<const T>() const { return reinterpret<const T>(*this); }
	operator device_ptr<const void>() const {
	    return reinterpret<const void>(*this);
	}
	operator device_ptr<void>() const { return reinterpret<void>(*this); }
    };

    template<typename T>
    device_ptr<T> wrap(T *data) {
	device_ptr<T> ptr;
	return (ptr.data_ = data, ptr);
    }


    template<typename T, typename U>
    device_ptr<T> reinterpret(device_ptr<U> ptr) {
	return wrap(reinterpret_cast<T*>(ptr.data()));
    }

    template<typename T>
    device_ptr<T> cast(device_ptr<const T> ptr) {
	return wrap(const_cast<T*>(ptr.data()));
    }

    template<typename T, size_t N>
    struct array_base;

    namespace detail {

	/**
	 * void type traits
	 * non-void type traits are partially specialized.
	 * @tparam T datatype
	 * @tparam R reference type
	 */
	template<typename T>
	struct type_traits {
	    typedef T& reference;
	    static const size_t size = sizeof(T);
	};

	/**
	 * non-void type traits.
	 * @tparam T datatype
	 */
	template<>
	struct type_traits<const void> {
	    static const size_t size = 1; /**< type size */
	};

	/**
	 * non-void type traits.
	 * @tparam T datatype
	 */
	template<>
	struct type_traits<void> {
	    static const size_t size = 1; /**< type size */
	};

	template<typename T>
	size_t sizeof_() { return type_traits<T>::size; }

	template<typename T>
	size_t sizeof_(size_t size) { return sizeof_<T>()*size; }
	
	/** 
	 * allocate device memory
	 * @param size  memory size in bytes
	 * @return device pointer
	 */
	void* malloc(size_t size);

	/** 
	 * free device memory
	 * @param ptr device pointer
	 */
	void free(void *ptr);

	/** 
	 * copy from device to host
	 * @param from source device pointer
	 * @param to destination host pointer
	 * @param size size in bytes
	 */
	void copy(device_ptr<const void> from, void *to, size_t size);

	/** 
	 * copy from host to device
	 * @param from  source host pointer
	 * @param to   destination device pointer
	 * @param size  size in bytes
	 */
	void copy(const void *from, device_ptr<void> to, size_t size);

	/** 
	 * copy from device to device symbol
	 * @param from source device pointer
	 * @param symbol destination device symbol
	 * @param size  size in bytes
	 */
	void copy(device_ptr<const void> from, const std::string &symbol, size_t size);

	/** 
	 * set device memory to byte value
	 * @param ptr device pointer
	 * @param size  size in bytes
	 * @param value  byte value
	 */
	void set(device_ptr<void> ptr, size_t size, char value);

    }

    template<typename T>
    device_ptr<T> malloc(size_t size) {
	size_t bytes = detail::sizeof_<T>(size);
	return wrap(reinterpret_cast<T*>(detail::malloc(bytes)));
    }

    template<typename T>
    void free(device_ptr<T> ptr) { detail::free(ptr.data()); }

    inline void copy(const void* from, device_ptr<void> to, size_t size) {
	detail::copy(from, to, size);
    }

    inline void copy(const device_ptr<const void> from, void *to, size_t size) {
	detail::copy(from, to, size);
    }

    template<typename T>
    void copy(const T *from, device_ptr<T> to, size_t size) {
	detail::copy(from, to, detail::sizeof_<T>(size));
    }

    template<typename T>
    void copy(const device_ptr<T> from, T *to, size_t size) {
	detail::copy(from, to, detail::sizeof_<T>(size));
    }

    template<typename T>
    void copy(const device_ptr<T> from, const std::string &symbol, size_t size) {
	size_t bytes = detail::sizeof_<T>(size);
	detail::copy(device_ptr<const void>(from), symbol, bytes);
    }

    template<typename T>
    void fill(device_ptr<T> ptr, size_t size, char value) {
	fill(device_ptr<char>(ptr), size* sizeof(T), value);
    } 



    template<typename T, size_t N>
    struct array_base {
	size_t size() const { return size_; }
	void fill(char value) { cudaxx::fill(this->data(), N*size_, value); }
	device_ptr<T> data() { return data_; }
	device_ptr<const T> data() const { return data_; }
	device_ptr<T> operator+(size_t offset) { return data_ + offset*N; }
	device_ptr<const T> operator+(size_t offset) const { return data_ + offset*N; }
	T  operator[](size_t i) const { return *((data_ + i*N).data()); }
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

#endif /* _CUDAXX_HPP_ */

#ifndef _RYSQ_CUDA_HPP_
#define _RYSQ_CUDA_HPP_

#include "rysq.hpp"


#ifndef _CUDAXX_CUDAXX_HPP_
#include "cuda/cudaxx.hpp"
#endif

namespace rysq {
    namespace cuda {

	using namespace ::cudaxx;


	static const int DEFAULT_DEVICE = -1;

	int initialize(int device = DEFAULT_DEVICE);

	void zero(void *gPtr, size_t bytes);

	void add(double *dest, const double *source, size_t size);

	template< typename T>
	void zero(T *gPtr, size_t size) { zero((void*)gPtr, size*sizeof(T)); }

	class Eri {
	public:
	    class Kernel;
	    Eri(const rysq::Quartet<rysq::Shell> &shells);
	    ~Eri();
	    operator bool() const { return kernel_ != NULL; }
	    void operator()(const cuda::array<double,3> &centers,
			    const cuda::array<int,4> &quartets,
			    double *Eri, const Parameters &parameters);
	private:
	    Kernel *kernel_;
	};

	class Fock {
	public:
	    struct Kernel;
	    Fock(const rysq::Quartet<rysq::Shell> &shells);
	    ~Fock();
	    operator bool() const { return eri_ != NULL; }
	    // operator bool() const { return  kernel_ != NULL; }
	    void operator()(const cuda::array<double,3> &centers,
			    const std::vector<Int4> &quartets,
			    const matrix::Adapter &D,
			    matrix::Adapter &F,
			    const Parameters &parameters);
	private:
	    Eri::Kernel *eri_;
	    Kernel *kernel_;
	};
	    

    }
}

#endif /* _RYSQ_CUDA_HPP_ */
