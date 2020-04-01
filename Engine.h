#pragma once

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdio>

#include "Vector.h"
#include "Sphere.h"

class Engine
{
public:
	explicit Engine(const char* fname, int save_int) : save_interval{ save_int } {
		init_system(fname); 
	};


	void step();
	void phase_plot(std::ostream& os);
	double total_kinetic_energy();
	std::vector<Vector> get_particle_positions();

protected:
	Vector G{ null };
	double Time{ 0 };
	double timestep;

	int nstep;
	int nprint;
	int nenergy;
	int save{ 0 };
	int save_interval;

	double lx;
	double ly;
	double x_0;
	double y_0;

	unsigned int no_of_particles;
	//std::ofstream os{ "C:/Users/james/Data/output.dump"};
	FILE* f1 = std::fopen("C:/Users/james/Data/output.dump", "w");

	

	std::vector<Sphere> particle;

	void integrate();
	void init_system(const char* fname);
	virtual void make_forces();
	void check_dump();
	void dump();
	
private:

	int bufferLimit{ 10 };
	int buffer{ 0 };



};

