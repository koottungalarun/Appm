#pragma once

#include <Eigen/Dense>
#include "FluidState.h"

class Numerics
{
public:
	Numerics();
	~Numerics();

	static Eigen::VectorXd fluidFlux_rusanov(const Eigen::VectorXd & qL_3d, const Eigen::VectorXd & qR_3d, const Eigen::Vector3d & fn, const double dx, double & dt_loc);
	static Eigen::VectorXd fluidFlux_rusanov(const Eigen::Vector3d & qL, const Eigen::Vector3d & qR);
	static Eigen::Vector3d fluidFlux(const Eigen::Vector3d & q);

	static double maxWaveSpeed(const Eigen::Vector3d & q);
};

