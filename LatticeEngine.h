#pragma once
#include "Engine.h"

/**
No lattice site is occupied by the center of more than one particle.

Consequentially, the size of the lattice sites gk must not be larger than sqrt(2)*Rmin. -> should be slightly smaller due to not considering deformation.

Each lattice site is assigned a number which is equal to the index of the particle whose center resides at this site or -1 if there is no such particle.

The diameter of the largest particle of the system covers gm = int(2*rmax/gx)+1 lattice sites in each dimension.

The centers of possible collision partners of a particle at site (i, j) are
hence located at the sites (i + k, j + l), with k, l = +-1, +-2, ..., +-gm.

To find the possible collision partneers of a particle we have to loop through all those sites (i + k, j + l) and
check whether any of the site values (pindex) is non-negative.

If such sites are found, the corresponding particles are inserted into the list of particles interacting with the one at hand.

This interaction list is constructed at the beginning of each time step.

This is order (N) complexity.

**/

class LatticeEngine :
	public Engine
{
public:
	explicit LatticeEngine(const char* fname, int save_interval) : Engine(fname, save_interval) {
		init_algorithm();
	}
	void step();

private:
	void make_ilist();
	void clear_pindex();
	void init_algorithm();
	void make_forces() override;


	int nx, ny, gm;
	double rmin, rmax, gk;
	std::vector<std::vector<int>> partners, pindex;
};

