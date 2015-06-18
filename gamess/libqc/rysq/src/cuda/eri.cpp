#include "externals/cuda/host.hpp"

#include <rysq/cuda.hpp>
#include "cuda/detail.hpp"

namespace host = cuda::host;

rysq::cuda::Eri::Eri(const rysq::Quartet<rysq::Shell> &quartet) {
    Transpose transpose;
    try { this->kernel_ = new cuda::detail::Eri(quartet, transpose); }
    catch (std::exception&) { this->kernel_ = NULL; }
}

rysq::cuda::Eri::~Eri() {
    delete kernel_;
}

void rysq::cuda::Eri::operator()(const  Centers &centers,
				  const  Quartets &quartets,
				 rysq::cuda::array <double> &eri,
				 const rysq::Parameters &parameters) {

	// for (int i = 0; i < quartets.size(); ++i)
	//     std::cout << quartets[i] << "\n";

    // host::device_ptr<double> eri =
    // 	host::malloc<double>(quartets.size()* kernel_->quartet().size());
    (*kernel_)(centers.data(), detail::Quartets(quartets.data(), quartets.size()),
	       eri.data(), parameters);

}
