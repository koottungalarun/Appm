#include "FluidState.h"

const double FluidState::gamma = 1.4;

FluidState::FluidState()
{
}

FluidState::FluidState(const double rho, const Eigen::Vector3d rhoU, const double rhoEtot)
{
	this->massDensity = rho;
	this->momentumDensity = rhoU;
	this->totalEnergyDensity = rhoEtot;
}

FluidState::FluidState(const FluidPrimitiveState & primitiveState) 
{
	const double rho = primitiveState.rho;
	const double p = primitiveState.p;
	const Eigen::Vector3d u = primitiveState.u;
	const double e = p / rho * (1. / (gamma - 1.));
	const double etot = 0.5 * u.squaredNorm() + e;

	this->massDensity = rho;
	this->momentumDensity = rho * u;
	this->totalEnergyDensity = rho * etot;
}

FluidState::FluidState(const Eigen::VectorXd & q) 
{
	assert(q.size() == 5);
	this->massDensity = q(0);
	this->momentumDensity = q.segment(1, 3);
	this->totalEnergyDensity = q(4);
}


FluidState::~FluidState()
{
}

const FluidPrimitiveState FluidState::getPrimitiveState() const
{
	const double rho_times_e = totalEnergyDensity - 0.5 * momentumDensity.squaredNorm() / massDensity;
	const double p = (gamma - 1.) * rho_times_e;

	FluidPrimitiveState primitiveState;
	primitiveState.rho = massDensity;
	primitiveState.u = momentumDensity / massDensity;
	primitiveState.p = p;
	return primitiveState;
}

const Eigen::VectorXd FluidState::getVecState() const
{
	Eigen::VectorXd vecState(5);
	vecState(0) = massDensity;
	vecState.segment(1, 3) = momentumDensity;
	vecState(4) = totalEnergyDensity;
	return vecState;
}

FluidPrimitiveState::FluidPrimitiveState()
{
	this->p = 0.0;
	this->rho = 1.0;
	this->u = Eigen::Vector3d::Zero();
}

FluidPrimitiveState::FluidPrimitiveState(const double rho, const double p, const Eigen::Vector3d & u)
{
	this->rho = rho;
	this->p = p;
	this->u = u;
}
