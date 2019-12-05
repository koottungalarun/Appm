#pragma once

#include <Eigen/Dense>

struct FluidPrimitiveState {
	double rho;
	double p;
	Eigen::Vector3d u;

	FluidPrimitiveState();
	FluidPrimitiveState(const double rho, const double p, const Eigen::Vector3d & u);
};

class FluidState
{
public:
	FluidState();
	FluidState(const double rho, const Eigen::Vector3d rhoU, const double rhoEtot);
	FluidState(const FluidPrimitiveState & primitiveState);
	FluidState(const Eigen::VectorXd & q);
	~FluidState();

	static const double gamma;

	const FluidPrimitiveState getPrimitiveState() const;
	const Eigen::VectorXd getVecState() const;

private:
	double massDensity = 0;
	Eigen::Vector3d momentumDensity = Eigen::Vector3d::Zero();
	double totalEnergyDensity = 0;

};

