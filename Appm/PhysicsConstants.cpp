#include "PhysicsConstants.h"

PhysicsConstants::PhysicsConstants()
{}

PhysicsConstants::~PhysicsConstants() 
{}


PhysicsConstants & PhysicsConstants::instance()
{
	static PhysicsConstants _instance;
	return _instance;
}

const double PhysicsConstants::kB()
{
	return _kB;
}

const double PhysicsConstants::eps0()
{
	return _eps0;
}

const double PhysicsConstants::q()
{
	return _q;
}

const double PhysicsConstants::c0()
{
	return _c0;
}

const double PhysicsConstants::planckConstant()
{
	return _h;
}

const double PhysicsConstants::pi()
{
	return M_PI;
}
