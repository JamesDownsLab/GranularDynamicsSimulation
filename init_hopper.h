#pragma once
#pragma once
#include <iostream>

void dump_particle(std::ostream& os,
	double x, double y, double vx, double vy, double radius, double mass, int type);

int make_file();
