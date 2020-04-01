#include <fstream>
#include <stdlib.h>
#include <random>
#include <chrono>

double youngs{ 1e10 };
double mu{ 0.5 };
double A{ 0.01 };
double gamma_t{ 10 };
double rho{ 8 };


void dump_particle(std::ofstream& os,
	double x, double y, double vx, double vy, double radius, double mass, int type) {
	os.precision(11);
	os << x << "\t" << y << "\t0\t" << vx << "\t" << vy << "\t0\t"
		<< radius << "\t" << mass << "\t" << type << "\t"
		<< youngs << "\t" << A << "\t" << mu << "\t" << gamma_t
		<< "\t0\t0\t0"
		<< "\t0\t0\t0" << "\t0\t0\t0" << "\t0\t0\t0\n";
}

int make_file() {
	std::ofstream fout("closed_hopper.random");
	fout << "#gravity: 0 -9.81 0\n";
	fout << "#timestep: 1e-6\n";
	fout << "#nstep: 500000\n";
	fout << "#nprint: 1000\n";
	fout << "#nenergy: 1000\n";
	fout << "#Time: 0\n";
	fout << "#lx: 1\n";
	fout << "#ly: 1\n";
	fout << "#x_0: 0\n";
	fout << "#y_0: 0\n";
	for (int i = 0; i < 11; i++) {
		dump_particle(fout, 0.45 + i * 0.01, 0.19, 0, 0, 0.005, 1, 1);
	}
	for (int i = 0; i < 50; i++) {
		dump_particle(fout, 0.55 + (i + 0.5) * 0.005, i * 0.01 + 0.2, 0, 0, 0.005, 1, 1);
		dump_particle(fout, 0.45 - (i + 0.5) * 0.005, i * 0.01 + 0.2, 0, 0, 0.005, 1, 1);
	}
	const double Rmax = 0.006, Rmin = 0.004;

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine engine(seed);
	std::uniform_real_distribution<double> z_selector(0.0, 1.0);


	for (int i = 0; i < 30; i++) {
		for (int k = 0; k < 30; k++) {
			double centerx = 0.3115 + 0.013 * i;
			double centery = 0.6 + 0.013 * k;
			double z = z_selector(engine);
			double r = Rmin * Rmax / (Rmax - z * (Rmax - Rmin));
			dump_particle(fout, centerx, centery, 0, 0, r, r * r / (Rmax * Rmax), 0);
		}
	}
	return 0;
}
