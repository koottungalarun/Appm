#pragma once

#include <Eigen/Dense>

class Physics
{
public:
	Physics();
	~Physics();

	/** Ratio of heat capacities */
	static const double gamma;

	static const Eigen::VectorXd primitive2state(const double n, const double p, const Eigen::Vector3d u);

};

