#include "Sphere.h"
#include <assert.h>


// Implements Gear integration. 
// Computes the positions, velocities and higher-order time derivatives as Taylor expansions of the present values.
void Sphere::predict(double dt)
{
	double a1 = dt;
	double a2 = a1 * dt / 2;
	double a3 = a2 * dt / 3;
	double a4 = a3 * dt / 4;

	rtd0 += a1 * rtd1 + a2 * rtd2 + a3 * rtd3 + a4 * rtd4;
	rtd1 += a1 * rtd2 + a2 * rtd3 + a3 * rtd4;
	rtd2 += a1 * rtd3 + a2 * rtd4;
	rtd3 += a1 * rtd4;
}

void Sphere::correct(double dt, const Vector& G)
{
	static Vector accel, corr;
	double dtrez = 1 / dt;
	const double coeff0 = double(19) / double(90) * (dt * dt / double(2));
	const double coeff1 = double(3) / double(4) * (dt / double(2));
	const double coeff3 = double(1) / double(2) * (double(3) * dtrez);
	const double coeff4 = double(1) / double(12) * (double(12) * (dtrez * dtrez));

	accel = Vector((1 / _m) * _force.x() + G.x(),
		(1 / _m) * _force.y() + G.y(),
		(1 / _J) * _force.phi() + G.phi());

	corr = accel - rtd2;
	rtd0 += coeff0 * corr;
	rtd1 += coeff1 * corr;
	rtd2 = accel;
	rtd3 += coeff3 * corr;
	rtd4 += coeff4 * corr;

}

double Sphere::kinetic_energy() const
{
	return _m * (rtd1.x() * rtd1.x() / 2 + rtd1.y() * rtd1.y() / 2)
		+ _J * rtd1.phi() * rtd1.phi() / 2;
}

// Container walls are built of particles that behave differently
void Sphere::boundary_conditions(double timestep, double Time)
{
	switch (ptype())
	{
		// Normal particle that can move
	case(0): break;

	case(1): break; // Particle cannot move
	case(2): {
		x() = 0.5 - 0.4 * cos(10 * Time);
		y() = 0.1;
		vx() = 10 * 0.4 * sin(Time);
		vy() = 0;
	} break;
	case(3): {
		double xx = x() - 0.5;
		double yy = y() - 0.5;
		double xp = xx * cos(timestep) - yy * sin(timestep);
		double yp = xx * sin(timestep) + yy * cos(timestep);

		x() = 0.5 + xp;
		y() = 0.5 + yp;
		vx() = -yp;
		vy() = xp;
		omega() = 1;
	} break;
		//case(4): {
		//	x() = 0.5 + 0.1 * cos(Time) + 0.4 * cos(Time + 2 * n * 3.141 / 128);
		//	y() = 0.5 + 0.1 * sin(Time) + 0.4 * sin(Time + 2 * n * 3.141 / 128);
		//	vx() = -0.1 * sin(Time) - 0.4 * sin(Time + 2 * n * 3.141 / 128);
		//	vy() = 0.1 * cos(Time) - 0.4 * cos(Time + 2 * n * 3.141 / 128);
		//	omega() = 1;
		//} break;
	case(5): {
		y() = 0.1 + 0.02 * sin(30 * Time);
		vx() = 0;
		vy() = 0.02 * 30 * cos(30 * Time);
	} break;
		//case(6): {
		//	int i = n / 2;
		//	y() = i * 0.02 + 0.1 + 0.02 * sin(30 * Time);
		//	vx() = 0;
		//	vy() = 0.02 * 30 * cos(30 * Time);
		//} break;
	default: {
		std::cout << "\n" << *this << "\n";
		std::cerr << "ptype: " << ptype() << " not implemented\n";
		abort();
	}
	}
}

// Enforces the periodic boundary confitions.
// Places the particle at the position of its periodic image if it crossed the boundary
void Sphere::periodic_bc(double x_0, double y_0, double lx, double ly)
{
	while (rtd0.x() < x_0) rtd0.x() += lx;
	while (rtd0.x() > x_0 + lx) rtd0.x() -= lx;
	while (rtd0.y() < y_0) rtd0.y() += ly;
	while (rtd0.y() > y_0 + ly) rtd0.y() -= ly;
}

std::istream& operator>>(std::istream& is, Sphere& p)
{
	is >> p.rtd0 >> p.rtd1
		>> p._r >> p._m >> p._ptype
		>> p.Y >> p.A >> p.mu >> p.gamma
		>> p._force
		>> p.rtd2 >> p.rtd3 >> p.rtd4;
	p._J = p._m * p._r * p._r / 2;
	return is;
}

std::ostream& operator<<(std::ostream& os, const Sphere& p)
{
	os << p.rtd0 << " " << p.rtd1 << " ";
	os << p._r << " " << p._m << " " << p._ptype << " ";
	os << p.Y << " " << p.A << " " << p.mu << " " << p.gamma << " ";
	os << p._force << " ";
	os << p.rtd2 << " " << p.rtd3 << " " << p.rtd4 << "\n" << std::flush;
	return os;
}

// Computes the force which is exerted by two spheres
void force(Sphere& p1, Sphere& p2, double lx, double ly)
{
	double dx = normalize(p1.x() - p2.x(), lx);
	double dy = normalize(p1.y() - p2.y(), ly);
	double rr = sqrt(dx * dx + dy * dy);
	double r1 = p1.r();
	double r2 = p2.r();
	double xi = r1 + r2 - rr;

	if (xi > 0) {
		double Y = p1.Y * p2.Y / (p1.Y + p2.Y);
		double A = 0.5 * (p1.A + p2.A);
		double mu = (p1.mu < p2.mu ? p1.mu : p2.mu);
		double gamma = (p1.gamma < p2.gamma ? p1.gamma : p2.gamma);
		double reff = (r1 * r2) / (r1 + r2);
		double dvx = p1.vx() - p2.vx();
		double dvy = p1.vy() - p2.vy();
		double rr_rez = 1 / rr;
		double ex = dx * rr_rez;
		double ey = dy * rr_rez;
		double xidot = -(ex * dvx + ey * dvy);
		double vtrel = -dvx * ey + dvy * ex + p1.omega() * p1.r() - p2.omega() * p2.r();
		double fn = std::sqrt(xi) * Y * sqrt(reff) * (xi + A * xidot);
		double ft = -gamma * vtrel;

		if (fn < 0) fn = 0;
		if (ft < -mu * fn) ft = -mu * fn;
		if (ft > mu * fn) ft = mu * fn;
		if (p1.ptype() == 0) {
			p1.add_force(Vector(fn * ex - ft * ey, fn * ey + ft * ex, r1 * ft));
		}
		if (p2.ptype() == 0) {
			p2.add_force(Vector(-fn * ex + ft * ey, -fn * ey - ft * ex, -r2 * ft));
		}
	}
}
