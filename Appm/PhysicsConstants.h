#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
/**
* Class for physical constants. 
*
* Implementation note: Singleton pattern.
*/
class PhysicsConstants {

private:
	/** Boltzmann constant (J/K) */
	const double _kB = 1.38e-23;

	/** Electrical permittivity (F/m) */
	const double _eps0 = 8.85e-12;

	/** Elementary charge (C) */
	const double _q = 1.609e-19;

	/** Speed of light in vacuum (m/s) */
	const double _c0 = 299792458;

	/** Planck constant (Js) */
	const double _h = 6.626e-34; 

	/** Atomic mass (kg) */
	const double _u = 1.66053906660e-27;

	PhysicsConstants();
	PhysicsConstants(const PhysicsConstants &);
	PhysicsConstants & operator=(const PhysicsConstants &);

public:
	
	static PhysicsConstants & instance();
	~PhysicsConstants();

	const double kB();
	const double eps0();
	const double q();
	const double c0();
	const double planckConstant();
	const double pi();
	const double atomicMass();
};
