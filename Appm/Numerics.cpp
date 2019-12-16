#include "Numerics.h"



Numerics::Numerics()
{
}


Numerics::~Numerics()
{
}

Eigen::VectorXd Numerics::fluidFlux_rusanov(const Eigen::VectorXd & qL_3d, const Eigen::VectorXd & qR_3d, const Eigen::Vector3d & fn, const double dx, double & dt_loc)
{
	assert(qL_3d.size() == 5);
	assert(qR_3d.size() == 5);
	assert(dx > 0);
	Eigen::VectorXd flux_3d(5);

	Eigen::Vector3d qL;
	Eigen::Vector3d qR;

	qL(0) = qL_3d(0);
	qL(1) = fn.dot(qL_3d.segment(1, 3));
	qL(2) = qL_3d(4);

	qR(0) = qR_3d(0);
	qR(1) = fn.dot(qR_3d.segment(1, 3));
	qR(2) = qR_3d(4);

	const Eigen::Vector3d flux = fluidFlux_rusanov(qL, qR);
	flux_3d(0) = flux(0);
	flux_3d.segment(1, 3) = flux(1) * fn;
	flux_3d(4) = flux(2);

	const double sL = maxWaveSpeed(qL);
	const double sR = maxWaveSpeed(qR);
	const double s = std::max(sL, sR);

	dt_loc = dx / s;
	return flux_3d;
}

Eigen::VectorXd Numerics::fluidFlux_rusanov(const Eigen::Vector3d & qL, const Eigen::Vector3d & qR)
{
	const Eigen::Vector3d fL = fluidFlux(qL);
	const Eigen::Vector3d fR = fluidFlux(qR);
	const double sL = maxWaveSpeed(qL);
	const double sR = maxWaveSpeed(qR);
	double s = std::max(sL, sR);

	Eigen::Vector3d flux;
	flux = 0.5 * (fL + fR) - 0.5 * s * (qR - qL);
	return flux;
}


Eigen::Vector3d Numerics::fluidFlux(const Eigen::Vector3d & q)
{
	const double gamma = 1.4;
	Eigen::Vector3d flux;
	flux(0) = q(1);
	flux(1) = (gamma - 1) * q(2) + (3. - gamma) / 2. * pow(q(1), 2) / q(0);
	flux(2) = gamma * q(2) * q(1) / q(0) + (1. - gamma) / 2. * pow(q(1), 3) / pow(q(0), 2);
	return flux;
}

double Numerics::maxWaveSpeed(const Eigen::Vector3d & q)
{
	const double gamma = 1.4;
	const double rho = q(0);
	const double u = q(1) / q(0);
	const double p = (gamma - 1) * (q(2) - 0.5 * pow(q(1), 2) / q(0));
	const double c = sqrt(gamma * p / rho);
	return Eigen::Vector3d(u - c, u, u + c).cwiseAbs().maxCoeff();
}

Eigen::Vector3d Numerics::raviartThomasBasis(const int idx, const Eigen::Vector3d & x)
{
	switch (idx) {
	case 0:
		return Eigen::Vector3d(x(0), x(1), 0);
	case 1:
		return Eigen::Vector3d(x(0) - 1, x(1), 0);
	case 2:
		return Eigen::Vector3d(x(0), x(1) - 1, 0);
	case 3:
		return Eigen::Vector3d(0, 0, x(2) - 1);
	case 4:
		return Eigen::Vector3d(0, 0, x(2));
	default:	
		return Eigen::Vector3d::Zero();
	}
	return Eigen::Vector3d::Zero();
}


