#pragma once

#include <iostream>
#include "Vector.h"

inline double normalize(double dx, double L) {
	while (dx < -L / 2) dx += L;
	while (dx >= L / 2) dx -= L;
	return dx;
}

class Sphere
{
	friend std::istream& operator >> (std::istream& is, Sphere& p);
	friend std::ostream& operator << (std::ostream& os, const Sphere& p);
	friend double Distance(const Sphere& p1, const Sphere& p2, double lx, double ly) {
		double dx = normalize(p1.rtd0.x() - p2.rtd0.x(), lx);
		double dy = normalize(p1.rtd0.y() - p2.rtd0.y(), ly);
		return std::sqrt(dx * dx + dy * dy);
	}

	friend void force(Sphere& p1, Sphere& p2, double lx, double ly);

public:
	Sphere() : rtd0(null), rtd1(null), rtd2(null), rtd3(null), rtd4(null) {}
	Vector& pos() { return rtd0; }
	Vector pos() const { return rtd0; }
	double& x() { return rtd0.x(); }
	double x() const { return rtd0.x(); }
	double& y() { return rtd0.y(); }
	double y() const { return rtd0.y(); }
	double& phi() { return rtd0.phi(); }
	double phi() const { return rtd0.phi(); }
	double& vx() { return rtd1.x(); }
	double vx() const { return rtd1.x(); }
	double& vy() { return rtd1.y(); }
	double vy() const { return rtd1.y(); }
	double& omega() { return rtd1.phi(); }
	double omega() const { return rtd1.phi(); }

	const Vector& velocity() const { return rtd1; }
	double& r() { return _r; }
	double r() const { return _r; }
	double m() const { return _m; }
	int ptype() const { return _ptype; }
	void predict(double dt);
	void add_force(const Vector& f) { _force += f; }
	void correct(double dt, const Vector& G);
	double kinetic_energy() const;
	void set_force_to_zero() { _force = null; }
	void boundary_conditions(double timestep, double Time);
	void periodic_bc(double x_0, double y_0, double lx, double ly);

private:
	double _r; // radius
	double _m; // mass
	double _J; // moment of inertia
	int _ptype;
	double Y; // Young's modulus
	double A; // dissipative constant
	double mu; // Poisson ratio
	double gamma;

	// Position, velocity and higher order time derivatives
	Vector rtd0, rtd1, rtd2, rtd3, rtd4;
	Vector _force;
};

