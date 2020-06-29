#pragma once
/**
* Class for physical constants. 
*
* Implementation note: Singleton pattern.
*/
class PhysicsConstants {

private:
	/** Boltzmann constant */
	const double _kB = 1.38e-23;

	/** Electrical permittivity */
	const double _eps0 = 8.85e-12;

	/** Elementary charge */
	const double _q = 1.609e-19;


	PhysicsConstants();
	PhysicsConstants(const PhysicsConstants &);
	PhysicsConstants & operator=(const PhysicsConstants &);

public:
	
	static PhysicsConstants & instance();
	~PhysicsConstants();

	const double kB();
	const double eps0();
	const double q();
};
