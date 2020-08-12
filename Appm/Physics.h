#pragma once

#include <Eigen/Dense>
#include <iostream>

#include <exception>

class Physics
{
public:
	Physics();
	~Physics();

	/** Ratio of heat capacities */
	static const double gamma;

	static const Eigen::VectorXd primitive2state(const double massRatio, const double n, const double p, const Eigen::Vector3d u);
	static void state2primitive(const double massRatio, const Eigen::VectorXd & state, double & n, double & p, Eigen::Vector3d & u);
	
	static const double getMaxWavespeed(const Eigen::Vector3d & qL, const Eigen::Vector3d & qR);
	static double getMaxWavespeed(const Eigen::Vector3d & state);

	/** Get temperature from fluid state in a single cell. */
	static double getTemperature(const Eigen::VectorXd & state, const double massRatio);

	/** Get temperature from fluid states in a list of cells. */
	static Eigen::VectorXd getTemperature(const Eigen::MatrixXd & states, const double massRatio);


	static const Eigen::Vector3d getFluidFluxFromState(const Eigen::Vector3d & q);

	static const Eigen::Vector3d getRusanovFlux(const Eigen::Vector3d & qL, const Eigen::Vector3d & qR, const bool showOutput = false);


	static const double thermionicEmissionCurrentDensity(const double Ts, const double W);
};

