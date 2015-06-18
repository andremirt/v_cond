#ifndef _VECTOR_HPP_
#define _VECTOR_HPP_

#include <boost/array.hpp>

template<size_t n, typename T = double>
class Vector  {
    boost::array<T,n> elems __attribute__ ((aligned(16)));
public:
    Vector() {}

    Vector(const T *v) {
	for (uint i = 0; i < n; ++i) {
	    this->elems[i] = v[i];
	}
    }

    template<class C>
    Vector(const C &v) {
	for (uint i = 0; i < n; ++i) {
	    this->elems[i] = v[i];
	}
    }
    
    typedef T* pointer;
    typedef const T* const_pointer;
    operator pointer() { return elems.c_array(); }
    operator const_pointer() const { return elems.c_array(); }

    T* data() { return elems.data(); }
    const T* data() const { return elems.data(); }

    T& operator[](uint i)  { return elems[i]; }
    const T& operator[](uint i) const { return elems[i]; }

    T inner() const {
	T q = 0.0;
	for (uint i = 0; i < n; ++i) {
	    q += (*this)[i]*(*this)[i];
	}
	return q;
    }

    Vector operator-(const Vector &rhs) const {
	Vector v;
	for (uint i = 0; i < n; ++i) {
	    v[i] = (*this)[i] - rhs[i];
	}
	return v;
    }	

    Vector& operator/=(T a) { divide(*this, a); return *this; }

    static Vector center(T w1, const Vector &v1,
			 T w2, const Vector &v2) {
	Vector v;
	T w = 1.0/(w1 + w2);
	for (uint i = 0; i < n; ++i) {
	    v[i] = (w1*v1[i] + w2*v2[i])*w;
	}
	return v;
    }

};


namespace {

    template<size_t N, typename T>
    void divide(Vector<N,T> &v, T a) {
	for (uint i = 0; i < N; ++i) v[i] /= a;
    }

    template<typename T>
    void divide(Vector<0,T> &v, T a) { }
}

#endif /* _VECTOR_HPP_ */
