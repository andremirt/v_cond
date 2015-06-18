/**
   @file
   @brief Librysq GAMESS bindings
   @details The default integer size is 4 bytes. If RYSQ_WITH_INTEGER8 is defined, the integer size is 8 bytes.
   The size of arrays in common blocks is controlled by MXATM, MXGTOT, and MXSH parameters. If the corresponding GAMESS parameters are different, the parameters in this file must be changed to match
   those of GAMESS.
*/

#ifdef HAVE_CONFIG_H 
#include "config.h"
#endif

#include "gamess.h"

#include "core/basis.hpp"



extern "C" {

    /** \brief creates a new molecule from games, common block
	\return new molecule handle
    */
    Integer molecule_new_gamess() {
	int ian[INFOA.nat];

	std::copy(INFOA.ian, INFOA.ian + INFOA.nat, ian);

	//Create a molecule

	int mol = Molecule_new(INFOA.nat, ian, (double*) INFOA.c, INFOA.zan);

	//std::cout << Molecule_find_object(mol);

	return mol;
    }

    Integer molecule_new_gamess_();
#pragma weak molecule_new_gamess_ = molecule_new_gamess

    /** \brief creates a new basis from games, common block end a given molecule
	\param molecule molecule handle
	\return new basis handle
    */
    Integer basis_new_gamess(Integer *molecule) {
	int basis = Basis_new(*molecule);
	Basis* basisObj = Basis_find_object(basis);

	assert(basisObj && "Basis not found");

	double *C[7] = { NSHEL.cs, NSHEL.cp, NSHEL.cd, NSHEL.cf, NSHEL.cg, NSHEL.ch, NSHEL.ci };

	for (int i = 0; i < NSHEL.nshell; ++i) {
	    int kstart = NSHEL.kstart[i] - 1;
	    int K = NSHEL.kng[i];
	    int L = NSHEL.ktype[i] - 1;
	    int center = NSHEL.katom[i] - 1;

	    if (NSHEL.kmin[i] == 1 && NSHEL.kmax[i] == 4) {
		basisObj->add(K, &NSHEL.ex[kstart], &C[0][kstart], &C[1][kstart], center);
	    }
	    else {
		basisObj->add(K, &NSHEL.ex[kstart], &C[L][kstart], L, center);
	    }
	}

	//std::cout << *basisObj;
	return basis;
    }

    Integer basis_new_gamess_(Integer *molecule);
#pragma weak basis_new_gamess_ = basis_new_gamess

}
