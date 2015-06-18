#ifndef _RYSQ_CUDA_MEMORY_HPP_
#define _RYSQ_CUDA_MEMORY_HPP_

#include <stdlib.h>

namespace rysq {
    namespace cuda {

	template<typename T, class = void>
	struct device_ptr {
	    T *data_;
	    device_ptr() : data_(NULL) {}
	    explicit device_ptr(T *data) : data_(data) {}
	    operator device_ptr<const T>()  const {
		return device_ptr<const T>(data_);
	    }
	    operator bool()const { return data_; }
	    T* data() { return data_; }
	    const T* data() const { return data_; }
	};

	template<>
	struct device_ptr<const void>
	    : device_ptr <const void, device_ptr<const void> > {
	    typedef device_ptr<const void,device_ptr> base;
	    device_ptr() :  base(NULL) {}
	    explicit device_ptr(void *data) : base(data) {}
	    template<typename U>
	    device_ptr(device_ptr<U> ptr) : base(ptr.data()) {}
	};

	template<>
	struct device_ptr<void> : device_ptr <void, device_ptr<void> > {
	    typedef device_ptr<void,device_ptr> base;
	    device_ptr() :  base(NULL) {}
	    explicit device_ptr(void *data) : base(data) {}
	    template<typename U>
	    device_ptr(device_ptr<U> ptr) : base(ptr.data()) {}
	};

	device_ptr<void> malloc(size_t size);
	void free(device_ptr<void> ptr);

	void memcpy(device_ptr<const void> from, void *to, size_t size);
	void memcpy(const void *from, device_ptr<void> to, size_t size);
	void memset(device_ptr<void> from, char byte, size_t size);

	template<typename T>
	device_ptr<T> malloc( size_t size) {
	    void *data = cuda::malloc(size*sizeof(T)).data();
	    return device_ptr<T>(reinterpret_cast<T*>(data));
	}		

	// template<typename T>
	// void free(device_ptr<T> ptr) {
	//     cuda::free(ptr);
	// }		

	template<typename T>
	void memcpy(device_ptr<const T> from, T *to, size_t size) {
	    memcpy(device_ptr<const void>(from), to, size*sizeof(T));
	}

	template<typename T>
	void memcpy(const T *from, device_ptr<T> to, size_t size) {
	    memcpy(from, device_ptr<void>(to), size*sizeof(T));
	}

	template<typename T>
	struct array {
	    array() : size_(0), capacity_(0), data_(NULL) {}
	    explicit array(const std::vector<T> &centers)
		: size_(0), capacity_(0), data_(NULL) {
		assign(centers);
	    }
	    explicit array(size_t size){
		resize(size);
	    }
	    array(const array & rhs) { operator=(rhs); }
	    array& operator=( const  array &rhs) {
		if (data_) throw;
		size_ = 0;
		capacity_ = 0;
		return *this;
	    }		
	    ~array() { if ( data_)  cuda::free(data_); }
	    size_t size() const { return size_; }
	    void reserve(size_t size){
		if (capacity_ < size) {
		    free(data_);
		    data_ = device_ptr<T>(NULL);
		}
		if (!data()) {
		    capacity_ = size;
		    data_ = malloc<T>(size);
		}
	    }
	    size_t capacity() const { return capacity_; }
	    void resize(size_t size) {
		reserve(size);
		size_ = size;
	    }
	    void assign(const std::vector<T> &v) {
		resize(v.size());
		if (!v.empty()) memcpy( &v.front(),  data_, v.size());
	    }
	    void clear() { size_ = 0 ; }	    
	    T* data() { return data_.data(); }
	    const T* data() const { return data_.data(); }
	    operator device_ptr<T>() { return data_; }
	    operator device_ptr<const T>() const { return data_; }
	    operator device_ptr<void>() { return data_; }
	    operator device_ptr< const void>()  const { return data_; }
	private:
	    size_t size_;
	    size_t capacity_;
	    device_ptr<T> data_;
	};

    }
}

#endif /* _RYSQ_CUDA_MEMORY_HPP_ */
