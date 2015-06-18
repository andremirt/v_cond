#include "cuda/quadrature.hpp"
#include "cuda/kernels/fock2.hpp"

using namespace rysq;

typedef rysq::cuda::Quadrature::Fock2 Fock2;

Fock2::Fock2(const ShellQuartet &shells, const Parameters &parameters) {
    kernel = cuda::kernels::fock2;
}

void Fock2::evaluate(const HFMatrix &gD, HFMatrix &gF,
		     const Array<int,4> &gQuartets,
		     const Parameters &parameters) {

    //kernel(gD, gF, gQuartets, gEri, parameters);

}


void Fock2::evaluate(const HFMatrix &gD, HFMatrix &gF,
		     const Array<int,4> &gQuartets, const double *gEri,
		     const Parameters &parameters) {
    kernel(gD, gF, gQuartets, gEri, parameters);
}
