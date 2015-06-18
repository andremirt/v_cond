#ifndef _CUDA_CORE_HPP_
#define _CUDA_CORE_HPP_

#include <string>

namespace cuda {
    namespace host {

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
	    void copy(device_ptr<const void> from,
		      const std::string &symbol, size_t size);

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
	void set(device_ptr<T> element, const T &value) {
	    copy(value, element, 1);
	}

	template<typename T>
	T get(const device_ptr<T> element) {
	    T value;
	    copy(element, &value, 1);
	    return value;
	}

	template<typename T>
	void copy(const device_ptr<T> from, const std::string &symbol, size_t size) {
	    size_t bytes = detail::sizeof_<T>(size);
	    detail::copy(device_ptr<const void>(from), symbol, bytes);
	}

	template<typename T>
	void memset(device_ptr<T> ptr, size_t size, char value) {
	    memset(device_ptr<char>(ptr), size* sizeof(T), value);
	} 


	template<typename T>
	struct device_reference {
	    static device_reference bind(device_ptr< T> data) {
		return device_reference(data);
	    }
	    operator T() const {
		return get(data_);
	    } 
	    device_reference& operator=(const T&value) {
		set(data_, value);
		return *this;
	    } 
	private:
	    device_ptr<T> data_;
	    device_reference(device_ptr< T> data) : data_(data) {}
	};


    }
}

#endif /* _CUDA_CORE_HPP_ */
