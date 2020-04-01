#pragma once

#include <math.h>
#include <iostream>

class Vector
{
	friend std::istream& operator >> (std::istream& is, Vector& v) {
		is >> v._x >> v._y >> v._phi;
		return is;
	}
	friend std::ostream& operator << (std::ostream& os, const Vector& v) {
		os << v._x << " " << v._y << " " << v._phi;
		return os;
	}
	friend Vector operator + (const Vector& v1, const Vector& v2) {
		Vector res(v1);
		res += v2;
		return res;
	}
	friend Vector operator - (const Vector& v1, const Vector& v2) {
		Vector res(v1);
		res -= v2;
		return res;
	}
	friend Vector operator * (double c, const Vector& p) {
		Vector res = p;
		res *= c;
		return res;
	}
	friend Vector operator * (const Vector& p, double c) {
		Vector res = p;
		res *= c;
		return res;
	}
	friend double norm2d(const Vector& v) {
		return std::sqrt(v._x * v._x + v._y * v._y);
	}
	friend double scalprod2d(const Vector& v1, const Vector& v2) {
		return v1._x * v2._x + v1._y * v2._y;
	}
	friend double vecprod2d(const Vector& v1, const Vector& v2) {
		return v1._x * v2._y - v1._y * v2._x;
	}

public:
	explicit Vector(double x = 0, double y = 0, double phi = 0) : _x(x), _y(y), _phi(phi) {};

	double& x() { return _x; }
	double x() const { return _x; }
	double& y() { return _y; }
	double y() const { return _y; }
	double& phi() { return _phi; }
	double phi() const { return _phi; }

	const Vector& operator += (const Vector& p) {
		_x += p._x;
		_y += p._y;
		_phi += p._phi;
		return *this;
	}

	const Vector& operator -= (const Vector& p) {
		_x -= p._x;
		_y -= p._y;
		_phi -= p._phi;
		return *this;
	}

	const Vector& operator *= (double c) {
		_x *= c;
		_y *= c;
		_phi += c;
		return *this;
	}

private:
	double _x, _y, _phi;
};

const Vector null(0, 0, 0);