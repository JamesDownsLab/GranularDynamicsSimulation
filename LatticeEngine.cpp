#include "LatticeEngine.h"
#include <algorithm>

void LatticeEngine::step()
{
	make_ilist();
	integrate();
	check_dump();
}

void LatticeEngine::make_ilist()
{
	for (unsigned int i{ 0 }; i < particle.size(); i++) {
		double x{ particle[i].x() }, y{ particle[i].y() };
		if ((x >= x_0) && (x < x_0 + lx) && (y >= y_0) && (y < y_0 + ly)) {
			int ix = int((x - x_0) / gk);
			int iy = int((y - y_0) / gk);
			pindex[ix][iy] = i;
			partners[i].clear();
		}
	}

	// for each particle i located at (ix, iy), the adjacent sites up to the distance
	// +- gm are scanned for possible interaction partners and these particles are 
	// recorded in the vector partners[i].
	for (unsigned int i{ 0 }; i < particle.size(); i++) {
		double x{ particle[i].x() }, y{ particle[i].y() };
		if ((x >= x_0) && (x < x_0 + lx) && (y >= y_0) && (y < y_0 + ly)) {
			int ix = int((x - x_0) / gk);
			int iy = int((y - y_0) / gk);
			for (int dx{ -gm }; dx < gm; dx++) {
				for (int dy{ -gm }; dy < gm; dy++) {
					int i0 = (ix + dx + nx) % nx;
					int i1 = (iy + dy + ny) % ny;
					int k = pindex[i0][i1];
					if (k > (int)i) {
						partners[i].push_back(k);
					}
				}
			}
		}
	}
	clear_pindex();
}

void LatticeEngine::clear_pindex()
{
	std::for_each(pindex.begin(), pindex.end(), [](auto& vec) {std::fill(vec.begin(), vec.end(), -1); });
}

void LatticeEngine::init_algorithm()
{
	rmin = particle[0].r();
	rmax = particle[0].r();
	for (unsigned int i = 1; i < particle.size(); i++) {
		if (particle[i].r() < rmin) {
			rmin = particle[i].r();
		}
		if (particle[i].r() > rmax) {
			rmax = particle[i].r();
		}
	}
	gk = sqrt(2) * rmin;
	gm = int(2 * rmax / gk) + 1;
	nx = int(lx / gk) + 1;
	ny = int(ly / gk) + 1;

	partners.resize(no_of_particles);
	pindex.resize(nx);
	for (unsigned int ix{ 0 }; ix < pindex.size(); ix++) {
		pindex[ix].resize(ny);
	}
	clear_pindex();
	make_ilist();
}

void LatticeEngine::make_forces()
{
	for (unsigned int i{ 0 }; i < no_of_particles; i++) {
		for (unsigned int k{ 0 }; k < partners[i].size(); k++) {
			int pk = partners[i][k];
			force(particle[i], particle[pk], lx, ly);
		}
	}
}

