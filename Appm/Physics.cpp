#include "Physics.h"

/** Ratio of heat capacities */
const double Physics::gamma = 1.4;

Physics::Physics()
{
}


Physics::~Physics()
{
}

const Eigen::VectorXd Physics::primitive2state(const double n, const double p, const Eigen::Vector3d u)
{
	assert(n > 0);
	assert(p > 0);
	const double e = 1. / (gamma - 1) * p / n; // ideal gas law
	assert(e > 0);
	const double etot = e + 0.5 * u.squaredNorm();
	assert(etot > 0);
	Eigen::VectorXd state = Eigen::VectorXd::Zero(5);
	state(0) = n;
	state.segment(1, 3) = n * u;
	state(4) = n * etot;
	return state;
}
