#pragma once

#include <Eigen/Dense>
#include <iostream>

class Physics
{
public:
	Physics();
	~Physics();

	/** Ratio of heat capacities */
	static const double gamma;

	static const Eigen::VectorXd primitive2state(const double n, const double p, const Eigen::Vector3d u);
	
	static double getMaxWavespeed(const Eigen::Vector3d & state);

	static const Eigen::Vector3d getFluidFluxFromState(const Eigen::Vector3d & q);

	static const Eigen::VectorXd getRusanovFlux(const Eigen::VectorXd & qL, const Eigen::VectorXd & qR, const bool showOutput = false);
};

