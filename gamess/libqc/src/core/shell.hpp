#ifndef _BASIS_SHELL_HPP_
#define _BASIS_SHELL_HPP_

#include <vector>
#include <iostream>

namespace basis {

    /**
       @brief
    */
    class Shell {

	friend class Basis;

	/**
	   @brief
	   @param output output stream
	   @param shell shell
	   @return output shell stream
	*/
	

    public:
	class Function {
	public:
	    static double computeC(int i, int j, int k);
	    
	    Function(int l, int m, int n):
		mask(l + (m << 8) + (n << 16)),
		C(computeC(l, m, n)) {}
	    
	    int mask;
	    double C;

	    int operator() (int i) const { return (0xFF & (mask >> (i*8))); }
	};

	const Function& function(int i) const { return functions.at(i); }

	/**
	   @brief Shell Constructormocked in
	   @param K Number of primitives
	   @param exps Primitive exponents, exps[K]
	   @param nc Number of contractions
	   @param coeff Contraction coefficients, coeff[nc][K]
	   @param L Shell angular momentums, L[nc]
	*/
	Shell(int K, double *exps, double *C, int L);

	/**
	   @brief SP shell constructor
	   @param exps exponents
	   @param Cs S shell coefficients
	   @param Cp P shell coefficients
	*/
	Shell(int K, double *exps, double *Cs, double *Cp);

	/** @brief Get number of primitives */
	int K() const { return _K; }

	int L() const { return _Lmax; }

	/** @brief Get max shell angular momentum */
	int Lmax() const { return _Lmax; }

	/**
	   @brief maximum angular momentum
	   @return maximum angular momentum
	*/
	int Lmin() const { return _Lmin; }

	/** @brief Get primitive exponents */
	const double* exps() const { return &(_exponents.at(0)); }

	std::vector<const double*> coeffs() const {
	    std::vector<const double*> C;
	    for (int i = 0; i < _nc; ++i) C.push_back(coeffs(i));
	    return C;
	}

	/** @brief Get contraction coefficients */
	const double* coeffs(int C) const { return &(_coefficients.at(C).at(0)); }

	/** shell size number of functions in shell */
	int size() const { return _size; }

	/**
	   @brief Get a primitive exponent
	   @param exp index of the primitive
	   @return primitive exponent
	*/
	double operator()(int prim) const { return _exponents[prim]; }

	int nc() const { return _nc; };

	/**
	   @brief Get contraction coefficient of a shell for a particular primitive
	   @param exp index of the primitive
	   @param cn
	*/
	double operator()(int prim, int cn) const { return _coefficients[cn][prim]; }

	/** @brief Get type of shell*/
	int type() const { return _type[0]; }

	bool operator==(const Shell &shell) const;

    private:
	/** number of primitives */
	int _K;
	/** exponents */
	std::vector<double> _exponents;
	/** number of contractions */
	int _nc;
	/** contraction coefficients */
	std::vector<std::vector<double> > _coefficients;
	/** angular momentums */
	std::vector<int> _type;
	/** Max shell angular momentum */
	int _Lmax;
	int _Lmin;
	int _size;

	std::vector<Function> functions;
    };

    struct shell {
	static size_t size(const Shell &shell) { return shell.size(); }
    };

}

std::ostream& operator<<(std::ostream &output, const basis::Shell &shell);

#endif /* _BASIS_SHELL_HPP_ */
