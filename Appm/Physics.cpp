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

double Physics::getMaxWavespeed(const Eigen::Vector3d & state)
{
	const double n = state(0);
	const double u = state(1) / state(0);
	const double etot = state(2) / state(0);
	const double s2 = gamma * (gamma - 1) * (etot - 0.5 * pow(u, 2));
	assert(std::isfinite(s2));
	assert(s2 > 0);
	const double s = sqrt(s2);
	assert(s > 0);
	const double smax = (Eigen::Vector3d(-1, 0, 1).array() * s + u).cwiseAbs().maxCoeff();
	return smax;
}

const Eigen::Vector3d Physics::getFluidFluxFromState(const Eigen::Vector3d & q)
{
	assert(q.norm() > 0);
	const double q1_q0 = q(1) / q(0);
	Eigen::Vector3d flux;
	flux(0) = q(1);
	flux(1) = 0.5 * (3 - gamma) * q1_q0 * q(1) + (gamma - 1) * q(2);
	flux(2) = gamma * q1_q0 * q(2) - 0.5 * (gamma - 1) * pow(q1_q0, 2) * q(1);
	return flux;
}

const double Physics::getMaxWavespeed(const Eigen::Vector3d & qL, const Eigen::Vector3d & qR) 
{
	const double sL = Physics::getMaxWavespeed(qL);
	const double sR = Physics::getMaxWavespeed(qR);
	assert(isfinite(sL));
	assert(isfinite(sR));
	const double s = std::max(sL, sR);
	assert(s > 0);
	assert(isfinite(s));
	return s;
}

const Eigen::Vector3d Physics::getRusanovFlux(const Eigen::Vector3d & qL, const Eigen::Vector3d & qR, const bool showOutput)
{
	const Eigen::Vector3d fL = Physics::getFluidFluxFromState(qL);
	const Eigen::Vector3d fR = Physics::getFluidFluxFromState(qR);
	const double s = Physics::getMaxWavespeed(qL, qR);
	assert(s > 0);
	const Eigen::Vector3d flux = 0.5 * (fL + fR) - 0.5 * s * (qR - qL);
	if (showOutput) {
		std::cout << "qL:\t" << qL.transpose() << std::endl;
		std::cout << "qR:\t" << qR.transpose() << std::endl;
		std::cout << "fL:\t" << fL.transpose() << std::endl;
		std::cout << "fR:\t" << fR.transpose() << std::endl;
		std::cout << "flux:\t" << flux.transpose() << std::endl;
	}
	return flux;
}

const double Physics::thermionicEmissionCurrentDensity(const double Ts, const double W)
{
	const double A = 0.5 * 1.20173e6;    // universal constant (units: A m^-2 K^-2)
	const double kB = 8.617333262145e-5; // Boltzmann constant (units: eV K^-1 )
	const double j_em = A * pow(Ts, 2) * exp(-W / (kB * Ts));
	return j_em;
}
