#include "core/shell.hpp"

#include "util/numeric.hpp"
#include <math.h>
#include <algorithm>
#include <iostream>

double basis::Shell::Function::computeC(int i, int j, int k) {
    int n = std::max(i + j + k, 1);
    int nx = std::max(i, 1);
    int ny = std::max(j, 1);
    int nz = std::max(k, 1);
    double denom = factorial2(2*nx -1 )*factorial2(2*ny -1 )*factorial2(2*nz -1 );
    return sqrt(factorial2(2*n - 1)/denom);
}

std::vector<basis::Shell::Function> generateFunctions(int L) {
    std::vector<basis::Shell::Function> result;

    for (int i = L; i >= 0; --i) {
	for (int k = 0; k <= (L-i)/2; ++k) {
	    int j = L - i - k;
	    if (i < j || j < k) continue;
	    permute(i, j, k, result);
	}
    }

    return result;
}

/** @see basis.h */
basis::Shell::Shell(int K, double *exps, double *C, int L) {
    _K = K;
    _exponents.assign(exps, exps + K);

    _nc = 1;
    _coefficients.assign(1, std::vector<double>(K));
    _coefficients[0].assign(C, C + K);

    _type.assign(1, L);
    _Lmax = L;
    _Lmin = L;

    _size = (((L+1)*(L+1) + (L+1))/2);

    functions = generateFunctions(L);
}

/** @see basis.h */
basis::Shell::Shell(int K, double *exps, double *Cs, double *Cp) {
    _K = K;

    _exponents.assign(exps, exps+ K);

    _nc = 2;
    _coefficients.assign(2, std::vector<double>(K));
    _coefficients[0].assign(Cs, Cs + K);
    _coefficients[1].assign(Cp, Cp + K);

    _type.assign(1, -1);
    _Lmax = 1;
    _Lmin = 0;
    _size = 4;

    for (int l = _Lmin; l <= _Lmax; ++l) {
	std::vector<basis::Shell::Function> list = generateFunctions(l);
	functions.insert(functions.end(), list.begin(), list.end());
    }

}

bool basis::Shell::operator==(const basis::Shell &shell) const {
    if (_K != shell._K) return false;
    else if (_type != shell._type) return false;
    else if (_nc != shell._nc) return false;
    else if (_exponents != shell._exponents) return false;
    else if (_coefficients != shell._coefficients) return false;

    return true;
}

std::ostream& operator<<(std::ostream &output, const basis::Shell &s) {
    output << "ptr = " << &s << "\t K = " << s.K() << "\t L = " << s.type() << std::endl;

//     for (unsigned int i = 0; i < s.functions.size(); ++i) {
// 	output << "("
// 	       << (s.functions.at(i))(0) << ", "
// 	       << (s.functions.at(i))(1) << ", "
// 	       << (s.functions.at(i))(2) << ")"
// 	       << "\t" << s.functions.at(i).C
// 	       << std::endl;
//     }

    for (int i = 0; i < s.K(); ++i) {
	output.width(13);
        output.precision(8);
	output << s(i);

	for (int j = 0; j < s.nc(); ++j) {
	    output.width(13);
	    output.precision(8);
	    output << '\t' << s(i, j);
	}

	output << std::endl;
    }

    return output;
}

