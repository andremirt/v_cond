#ifndef _RYSQ_ROOTS_EVALUATE_H_
#define _RYSQ_ROOTS_EVALUATE_H_


namespace {

    struct serial_tag {};

    template<typename T>
    __device__
    bool master(const T &thread){ return thread == 0; }

    __device__
    bool master(const serial_tag&) { return true; }

    template<typename T>
    __device__
    T index(const T &thread){ return thread; }

    __device__
    size_t index(const serial_tag&) { return 0; }

    template<int L> __device__
    static void evaluate1(const double X, double &t2, const double (&C)[L]) {
	double q = X;
	t2 = C[0] + C[1]*q;
	for (int i = 2; i < L; ++i) {
	    q *= X;
	    t2 += C[i]*q;
	}
    }

    template<int N, int L, typename T> __device__
    static void evaluate1(const double X, double *t2, const double (&C)[N][L],
				 const T& thread) {
	double q = X;
	t2[thread] = C[thread][0] + C[thread][1]*q;
	for (int i = 2; i < L; ++i) {
	    q *= X;
	    t2[thread] += C[thread][i]*q;
	}
    }

    template<int N>
    __device__
    static inline double evaluate_polynomial(const double (&C)[N], double x) {
	// evaluate1(x, t2, C, a);
	double p = 0;
	for (int i = N-1; i > 0; --i) {
	    p = (p + C[i])*x;
	}
	return (p + C[0]);
    }

    template<int N, int L, typename T>
    __device__
    static inline void evaluate_polynomial(double *t2,
					    const double (&C)[N][L], double x,
					    const T &thread) {
	evaluate1(x, t2, C, thread);
    }

    template<int N, int L>
    __device__
    static inline void evaluate_polynomial(double *t2,
					    const double (&C)[N][L], double x,
					   const serial_tag& = serial_tag()) {
	for (int a = 0; a < N; ++a) evaluate_polynomial(t2, C, x, a);
    }

    template<int N, int L> __device__
    static inline void evaluate1(const double X, double *t2, const double (&C)[N][L],
				 const serial_tag&) {
	for (int a = 0; a < N; ++a) evaluate1(X, t2, C, a);
    }

    template<int N, int L, int M, typename T>
    __device__
    static void evaluate1(const double &X, double *R, double *W,
			  const double (&C)[N][L], const double (&D)[N][M],
			  const T &thread) {
	double q = X;
	R[thread] = C[thread][0] + C[thread][1]*q;
	W[thread] = D[thread][0] + D[thread][1]*q;
	for (int i = 2; i < L; ++i) {
	    q *= X;
	    R[thread] += C[thread][i]*q;
	    W[thread] += D[thread][i]*q;
	}
	for (int i = L; i < M; ++i) {
	    q *= X;
	    W[thread] += D[thread][i]*q;
	}
    }

    template<int N, int L, int M> __device__
    static void evaluate1(const double &X, double *R, double *W,
				 const double (&C)[N][L], const double (&D)[N][M],
				 const serial_tag& = serial_tag()) {
	for (int a = 0; a < N; ++a) evaluate1(X, R, W, C, D, a);
    }

    template<int L, int M> __device__
    static inline void evaluate1(const double X, double &R, double &W,
				 const double (&C)[L], const double (&D)[M]) {
	double q = X;
	R = C[0] + C[1]*q;
	W = D[0] + D[1]*q;
	for (int i = 2; i < L; ++i) {
	    q *= X;
	    R += C[i]*q;
	    W += D[i]*q;
	}
	for (int i = L; i < M; ++i) {
	    q *= X;
	    W += D[i]*q;
	}
    }

    template<int N, typename T>
    __device__
    void scale_add(const double e, double *R, double *W,
		   const double X, const double (&RN)[N],
		   const double w, const double (&WN)[N],
		   const T &thread) {
	R[thread] = e*R[thread] + RN[thread]/(X - RN[thread]);
	W[thread] = e*W[thread] + w*WN[thread];
    }

    template<int N> __device__
    void scale_add(const double e, double *R, double *W,
		   const double X, const double (&RN)[N],
		   const double w, const double (&WN)[N],
		   const serial_tag& = serial_tag()) {
	for (int a = 0; a < N; ++a) scale_add(e, R, W, X, RN, w, WN, a);
    }

    template<int N> __device__
    void scale_add(const double e1, double *R, const double e2, double *W,
		   const double X, const double (&RN)[N],
		   const double w, const double (&WN)[N],
		   const serial_tag& = serial_tag()) {
	for (int a = 0; a < N; ++a) {
	    R[a] = e1*R[a] + RN[a]/(X - RN[a]);
	    W[a] = e2*W[a] + w*WN[a];
	}
    }

    template<int N, typename T>
    __device__
    void scale_add(const double e, double *t2, const double X, const double (&RN)[N],
		   const T &thread) {
	t2[thread] = e*t2[thread] + RN[thread]/(X - RN[thread]);
    }

    template<int N> __device__
    void scale_add(const double e, double *t2, const double X, const double (&RN)[N],
		   const serial_tag& = serial_tag()) {
	for (int a = 0; a < N; ++a) scale_add(e, t2, X, RN, a);
    }
 
    // template<int N> __device__
    // void scale_add(const double e, double *t2, double *W,
    // 		   const double X, const double (&RN)[N],
    // 		   const double w, const double (&WN)[N]) {
    // 	for (int a = 0; a < N; ++a) {
    // 	    t2[a] = e*t2[a] + RN[a]/(X - RN[a]);
    // 	    W[a] = e*W[a] + w*WN[a];
    // 	}
    // }

    __device__
    static inline void scale_add(const double e, double &R, double &W,
				 const double X, const double RN,
				 const double w, const double WN) {
	R = e*R + RN/(X - RN);
	W = e*W + w*WN;
    }

    __device__
    void scale_add(const double e1, double &R, const double e2, double &W,
		   const double X, const double (&RN),
		   const double w, const double (&WN)) {
	    R = e1*R + RN/(X - RN);
	    W = e2*W + w*WN;
    }

    template<int N> __device__
    double sum(const double *W) {
	double q = W[N-1];
	for (int i = N-2; i >= 0; --i) {
	    q += W[i];
	}
	return q;
    }

    // template<int N> __device__
    // static inline void scale(double *p, const double x, const double (&r)[N]){
    // 	for (int a = 0; a < N; ++a) {
    // 	    p[a] = x*r[a];
    // 	}
    // }

    // template<int N> __device__
    // static inline void divide(double *p, const double x, const double (&r)[N]){
    // 	for (int a = 0; a < N; ++a) {
    // 	    p[a] = r[a]/x;
    // 	}
    // }

    template<int N, typename T>
    __device__
    static inline void evaluate_inf(double *R, double r, const double (&CR)[N],
				    double *W, double w, const double (&CW)[N],
				    const T& thread, bool simt = true) {
	R[thread] = CR[thread]/(r - CR[thread]);
	W[thread] = w*CW[thread];
	if (simt && master(thread)) W[0] -= sum<N-1>(W+1);
    }

    template<int N>
    __device__
    static inline void evaluate_inf(double *R, double r, const double (&CR)[N],
				    double *W, double w, const double (&CW)[N],
				    const serial_tag& = serial_tag()) {
	for (int i = 0; i < N; ++i) evaluate_inf(R, r, CR, W, w, CW, i, false);
	W[0] -= sum<N-1>(W+1);
    }

    __device__
    static inline void change_variable(double &R) {
	    R = R/(1.0 + R);
    }

    template<int N, typename T> __device__
    static void change_variable(double *t2, const T &thread) {
	t2[thread] = t2[thread]/(1.0 + t2[thread]);
    }

    template<int N> __device__
    static void change_variable(double *t2, const serial_tag& = serial_tag()) {
	for(int i = 0; i < N; ++i) change_variable< N>(t2,i);
    }


    __device__
    static inline void scale(double &p, const double &x, const double &r){
	p = x*r;
    }

    __device__
    static inline void divide(double &p, const double &x, const double &r){
	p = r/x;
    }

    template<size_t N, typename T>
    __device__
    static inline void divide(double *p, const double &x, const double (&r)[N],
    			      const T & thread) {
    	p[thread] = r[thread]/x;
    }

    template<size_t N, typename T>
    __device__
    static inline void scale(double *p, const double &x, const double (&r)[N],
    			      const T & thread) {
    	p[thread] = r[thread]*x;
    }

    template<size_t N>
    __device__
    static void scale(double *p, const double &x, const double (&r)[N],
    		      const serial_tag &tag = serial_tag()) {
    	for (size_t i = 0; i < N; ++i) scale(p, x, r, i);
    }

    template<size_t N>
    __device__
    static void divide(double *p, const double &x, const double (&r)[N],
    		      const serial_tag &tag = serial_tag()) {
    	for (size_t i = 0; i < N; ++i) divide(p, x, r, i);
    }

    template<size_t N>
    __device__
    static double evaluate_inverse_polynomial(const double (&C)[N], double x) {
	double p = 0.0;
	for (size_t i = N - 1; i > 0; --i) {
	    p = (p + C[i])/x;
	}
	return (C[0] + p);
    }


    template<size_t M, size_t N, typename T>
    __device__
    static void evaluate_inverse_polynomial(double *R,
					    const double (&C)[M][N], double x,
					    const T & thread) {
	for (size_t i = 0; i < M; ++i)
	    R[i] = evaluate_inverse_polynomial(C[i], x);
    }


    template<size_t N, typename T>
    __device__
    static void scale(double s, double *R, const T &thread) {
	    R[thread] *= s;
    }

    template<size_t N>
    __device__
    static void scale(double s, double *R, const serial_tag& = serial_tag()) {
	for (size_t i = 0; i < N; ++i) scale<N>(s, R, i);
    }

    template<typename T> struct disable_if_serial { typedef void type; };
    template<> struct disable_if_serial<serial_tag> { };

    template<int M, int N, typename T>
    __device__
    static void //typename disable_if_serial<T>::type
    add_evaluate_polynomial1(double *R,
			     const double (&C)[M][N], double x,
			     const T &thread) {
	double p = 0.0;
	for (int i = N - 1; i >= 0; --i) {
	    p = (p + C[thread][i])*x;
	}
	R[thread] += p;
    }

    template<int M, int N>
    __device__
    static void add_evaluate_polynomial1(double *R,
					 const double (&C)[M][N], double x,
					 const serial_tag&) {
	for (int i = 0; i < M; ++i)
	    add_evaluate_polynomial1(R, C, x, i);
    }


    template<int N, typename T>
    __device__
    static inline void add_evaluate_inf(double *R, const double (&CR)[N], double r,
					const T&thread) {
	R[thread] += CR[thread]/(r - CR[thread]);
    }

    template<int N>
    __device__
    static inline void add_evaluate_inf(double *R, const double (&CR)[N], double r,
					 const serial_tag& = serial_tag()) {
	for (int i = 0; i < N; ++i)
	    add_evaluate_inf(R, CR, r, i);
    }


}

#endif /* _RYSQ_ROOTS_EVALUATE_H_ */
