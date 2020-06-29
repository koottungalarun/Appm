#include "Physics.h"

/** Ratio of heat capacities */
const double Physics::gamma = 1.4;

Physics::Physics()
{
}


Physics::~Physics()
{
}

const Eigen::VectorXd Physics::primitive2state(const double massRatio, const double n, const double p, const Eigen::Vector3d u)
{
	assert(n > 0);
	assert(p > 0);
	assert(massRatio > 0);
	const double n_e = 1. / (gamma - 1) * 1. / massRatio * p; // ideal gas law
	assert(n_e > 0);
	const double n_etot = n_e + 0.5 * n * u.squaredNorm();
	assert(n_etot > 0);
	Eigen::VectorXd state = Eigen::VectorXd::Zero(5);
	state(0) = n;
	state.segment(1, 3) = n * u;
	state(4) = n_etot;
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

double Physics::getTemperature(const Eigen::VectorXd & state, const double massRatio)
{
	double n = state(0);
	double ne_tot = state(5);
	assert(n > 0);
	assert(ne_tot > 0);
	const Eigen::Vector3d uvec = state.segment(1, 3) / n;
	double ne = ne_tot - 0.5 * n * uvec.norm();
	double T = (Physics::gamma - 1) * massRatio * ne;
	assert(T > 0);
	return T;
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

void Physics::state2primitive(const double massRatio, const Eigen::VectorXd & state, double & n, double & p, Eigen::Vector3d & u)
{
	assert(state.size() == 5);
	assert(massRatio > 0);
	n = state(0);
	if (!(n > 0)) {
		std::cout << n << std::endl;
		std::stringstream ss;
		ss << "Non-positive number density: n = " << n << std::endl;
		ss << "State: " << state.transpose() << std::endl;
		std::string str = ss.str();
		throw std::domain_error(str);
	}
	assert(n > 0);
	u = state.segment(1, 3) / n;
	const double n_etot = state(4);
	const double nu_sq = state.segment(1, 3).squaredNorm();
	const double ne = n_etot - 0.5 * nu_sq / n;
	p = (gamma - 1) * massRatio * ne;
	if (!(p > 0)) {
		std::stringstream ss;
		ss << "Non-positive pressure: p = " << p << std::endl;
		ss << "ne: " << ne << std::endl;
		ss << "State: " << state.transpose() << std::endl;
		ss << "mass ratio: " << massRatio << std::endl;
		std::string str = ss.str();
		throw std::domain_error(str);
	}
	assert(p > 0);
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
