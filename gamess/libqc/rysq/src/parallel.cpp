// #ifdef HAVE_CONFIG_H
// #include "config.h"
// #endif


// #include <assert.h>
// #include <stdlib.h>
// #include <algorithm>
// #include <vector>
// #include <boost/tuple/tuple.hpp>
// #include <boost/tuple/tuple_comparison.hpp>

// #ifdef _OPENMP
// #include <omp.h>
// #endif

// #include <rysq.hpp>
// #include "cxx/sugar.hpp"
// #include "quadrature.h"
// #include "memory.h"

// using namespace rysq;


// parallel::Fock::Fock(const Quartet<Shell> &shells) : shells(shells) {} 

// parallel::Fock::~Fock() {}

// void parallel::Fock::operator()(const std::vector<Center> &centers,
// 				const std::vector<Int4> &quartets,
// 				const matrix::Adapter &D,
// 				matrix::Adapter &F,
// 				const Parameters &parameters) {

//     if (quartets.empty()) return;

// #ifdef _OPENMP
// #pragma omp parallel 
// #endif
//     {


// #ifdef _OPENMP
// 	int num_threads = omp_get_num_threads();
// 	int thread_num = omp_get_thread_num();
// #else
// 	int num_threads = 1;
// 	int thread_num = 0;
// #endif

// 	rysq::Fock fock(this->shells);

// 	int nij = shells.size(0,1);
// 	int nkl = shells.size(2,3);
// 	int nik = shells.size(0,2);
// 	int nil = shells.size(0,3);
// 	int njk = shells.size(1,2);
// 	int njl = shells.size(1,3);

// 	memory::Map<double> myMemory(32*1024);

// #ifdef _OPENMP
// #pragma omp for schedule(dynamic,1) nowait
// #endif
// 	for (uint q = 0; q < quartets.size(); ++q) {

// 	    int i, j, k, l;
// 	    util::unpack(quartets[q], i, j, k, l);

// 	    double *Fij = myMemory.get(i, j, nij);
// 	    double *Fkl = myMemory.get(k, l, nkl);
// 	    double *Fik = myMemory.get(i, k, nik);
// 	    double *Fil = myMemory.get(i, l, nil);
// 	    double *Fjk = myMemory.get(j, k, njk);
// 	    double *Fjl = myMemory.get(j, l, njl);

// 	    Quartet<Center> C(centers[i], centers[j], centers[k], centers[l]);

// 	    matrix::Density<6> D6 = matrix::index(D, i, j, k, l);
// 	    double *f[] = { Fij, Fkl, Fik, Fil, Fjk, Fjl };
// 	    matrix::Fock<6> F6(f);
		
// 	    Parameters p = parameters;
// 	    p.scale /= Quartet<Shell>::symmetry(i, j, k, l);

// 	    // fock(C, D6, F6, p);
    
// 	}

// 	typedef boost::tuple<int,int,double*,size_t> Block;
// 	std::vector<Block> blocks;
// 	typedef std::vector<Block>::iterator block_it;

// 	while (!myMemory.empty()) {
// 	    memory::Map<double>::Entry kv = myMemory.pop();
// 	    blocks.push_back(Block(kv.j, kv.i, kv.ptr, kv.n));
// 	}
// 	std::sort(blocks.begin(), blocks.end());

// 	std::vector<block_it> columns;
// 	typedef std::vector<block_it>::iterator column_it;
// 	block_it it = blocks.begin();

// 	while (it < blocks.end()) {
// 	    int j = it->get<0>();
// 	    columns.push_back(it);
// 	    while (it < blocks.end() && it->get<0>() == j) ++it;
// 	}
// 	columns.push_back(it);	    

// 	for (int t = 0; t < num_threads; ++t) {
// 	    int k = (thread_num + t)%num_threads;

// 	    for (column_it col = columns.begin(); col < columns.end() - 1; ++col) {
// 		block_it block = *col;
// 		int j = block->get<0>(); 
// 		if (j%num_threads != k) continue;

// 		while (block < *(col+1)) {
// 		    double *dest = F(block->get<1>(), block->get<0>());
// 		    double *source = block->get<2>();
// 		    size_t size = block->get<3>();
// 		    util::add(size, dest, source);
// 		    ++block;
// 		}
// 	    }
// #ifdef _OPENMP
// #pragma omp barrier
// #endif
// 	}

//     }

// }
