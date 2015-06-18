#include "cuda/kernel/quadrature.hpp"
#include <exception>


// using namespace ::cuda;
using namespace rysq::cuda;

#define IMPL_TYPE kernel::Eri<kernel::eri::Transform>

rysq::cuda::detail::Eri::Eri(const rysq::Quartet<rysq::Shell> &quartet,
			     const rysq::Transpose &transpose) {
    this->impl_ = IMPL_TYPE::new_(quartet, transpose);
    if (!this->impl_) throw std::exception();
}

rysq::cuda::detail::Eri::~Eri() {
    delete static_cast<IMPL_TYPE*>(this->impl_);
}

