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
	const double _kB = 1.380649e-23;

	/** Electrical permittivity (F/m) */
	const double _eps0 = 8.8541878128e-12;

	/** Elementary charge (C) */
	const double _q = 1.602176634e-19;

	/** Speed of light in vacuum (m/s) */
	const double _c0 = 299792458;

	/** Planck constant (Js) */
	const double _h = 6.62607015e-34;

	/** Atomic mass (kg) */
	const double _u = 1.66053906660e-27;

	PhysicsConstants();
	PhysicsConstants(const PhysicsConstants &);
	PhysicsConstants & operator=(const PhysicsConstants &);

public:
	
	static PhysicsConstants & instance();
	~PhysicsConstants();

	/**
	* @return Boltzmann constant.
	*/
	const double kB();

	/**
	* @return electrical permittivity in vacuum.
	*/
	const double eps0();

	/**
	* @return elementary charge q = 1.609e-19 C.
	*/
	const double q();

	/**
	* @return vacuum speed of light.
	*/
	const double c0();

	/**
	* @return Planck constant.
	*/	
	const double planckConstant();

	/**
	* @return mathematical constant pi, ratio of circumference of a circle to its diameter.
	*/
	const double pi();

	/**
	* @return atomic mass constant, defined as 1/12 of carbon atomic mass.
	*/
	const double atomicMass();
};
