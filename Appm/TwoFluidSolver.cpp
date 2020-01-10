#include "TwoFluidSolver.h"



TwoFluidSolver::TwoFluidSolver()
{
}

TwoFluidSolver::TwoFluidSolver(const DualMesh * dualMesh)
	: FluidSolver(dualMesh)
{
}


TwoFluidSolver::~TwoFluidSolver()
{
}

void TwoFluidSolver::writeStates(H5Writer & writer) const
{
	const int nDualCells = mesh->getNumberOfCells();

	Eigen::VectorXd pressure_A(nDualCells);
	Eigen::VectorXd density_A(nDualCells);
	Eigen::Matrix3Xd velocity_A(3, nDualCells);
	
	Eigen::VectorXd pressure_B(nDualCells);
	Eigen::VectorXd density_B(nDualCells);
	Eigen::Matrix3Xd velocity_B(3, nDualCells);

	for (int i = 0; i < nDualCells; i++) {
		Eigen::VectorXd state = fluidStates.col(i);
		const PrimitiveState primitive = getPrimitiveState(state);

		pressure_A(i) = primitive.p_a;
		density_A(i) = massRatio_a * primitive.n_a;
		velocity_A.col(i) = primitive.u_a;

		pressure_B(i) = primitive.p_b;
		density_B(i) = massRatio_b * primitive.n_b;
		velocity_B.col(i) = primitive.u_b;
	}
	writer.writeData(pressure_A, "/pressureA");
	writer.writeData(density_A, "/densityA");
	writer.writeData(velocity_A, "/velocityA");
	writer.writeData(pressure_B, "/pressureB");
	writer.writeData(density_B, "/densityB");
	writer.writeData(velocity_B, "/velocityB");
}

const std::string TwoFluidSolver::getXdmfOutput(const int iteration) const
{
	// TODO
	return std::string();
}

void TwoFluidSolver::init()
{
	std::cout << "TwoFluidSolver::init()" << std::endl;
	const int nCells = mesh->getNumberOfCells();
	assert(nCells > 0);
	fluidStates = Eigen::MatrixXd::Zero(stateVectorLength, nCells);
	fluidFluxes = Eigen::MatrixXd::Zero(stateVectorLength, nCells);

	const Eigen::VectorXd qL = getInitStateLeft();
	const Eigen::VectorXd qR = getInitStateRight();

	for (int i = 0; i < nCells; i++) {
		const Cell * cell = mesh->getCell(i);
		const Eigen::Vector3d cellCenter = cell->getCenter();
		fluidStates.col(i) = (cellCenter(2) < 0) ? qL : qR;
	}

	std::cout << "Test for fluid flux: " << std::endl;
	double dt_loc = 0;
	const Eigen::Vector3d fn = Eigen::Vector3d::UnitZ();
	Eigen::VectorXd flux = getRusanovFlux(qL, qR, fn, 1, dt_loc);
	std::cout << "Flux:   " << flux.transpose() << std::endl;
	std::cout << "dt_loc: " << dt_loc << std::endl;

}

Eigen::VectorXd TwoFluidSolver::getRusanovFlux(const Eigen::VectorXd & qL, const Eigen::VectorXd & qR, const Eigen::Vector3d & fn, const double dx, double & dt_loc)
{
	// TODO
	return Eigen::VectorXd();
}

Eigen::VectorXd TwoFluidSolver::getFlux(const Eigen::VectorXd & q)
{
	// TODO
	Eigen::VectorXd flux(6);
	return flux;
}

const double TwoFluidSolver::maxWaveSpeed(const Eigen::Vector3d & q) const
{
	const double u = q(1) / q(0);
	const double s2 = gamma * (gamma - 1) * ( q(2) / q(0) - 0.5 *pow(q(1) / q(0), 2));
	assert(s2 > 0);
	const double s = sqrt(s2);
	return (u * Eigen::Vector3d::Ones() + s * Eigen::Vector3d(-1, 0, 1)).cwiseAbs().maxCoeff();
}

const double TwoFluidSolver::maxWaveSpeed(const Eigen::VectorXd & q, const Eigen::Vector3d & fn)
{
	assert(q.size() == stateVectorLength);
	const Eigen::Vector3d q_a = getState1d(q.segment(0, 5), fn);
	const double s_a = maxWaveSpeed(q_a);

	const Eigen::Vector3d q_b = getState1d(q.segment(5, 5), fn);
	const double s_b = maxWaveSpeed(q_b);
	return std::max(s_a, s_b);
}


TwoFluidSolver::PrimitiveState TwoFluidSolver::getPrimitiveState(const Eigen::VectorXd & q) const
{
	assert(q.size() == stateVectorLength);
	PrimitiveState primitive;
	primitive.n_a = q(0) / massRatio_a;
	primitive.u_a = q.segment(1, 3) / (q(0) * massRatio_a);
	primitive.p_a = (gamma - 1) * (q(4) - 0.5 * q.segment(1, 3).squaredNorm() / q(0));
	
	primitive.n_b = q(5) / massRatio_b;
	primitive.u_b = q.segment(6, 3) / (q(5) * massRatio_b);
	primitive.p_b = (gamma - 1) * (q(9) - 0.5 * q.segment(6, 3).squaredNorm() / q(9));
	return primitive;
}

Eigen::VectorXd TwoFluidSolver::getConservationState(const PrimitiveState & primitive) const
{
	Eigen::VectorXd q(10);
	q(0) = massRatio_a * primitive.n_a;
	q.segment(1, 3) = massRatio_a * primitive.u_a;
	q(4) = primitive.p_a / (gamma - 1) + 0.5 * q.segment(1, 3).squaredNorm() / q(0);

	q(5) = massRatio_b * primitive.n_b;
	q.segment(6, 3) = massRatio_b * primitive.u_b;
	q(9) = primitive.p_b / (gamma - 1) + 0.5 * q.segment(6, 3).squaredNorm() / q(5);
	return q;
}

const Eigen::VectorXd TwoFluidSolver::getInitStateLeft() const
{
	PrimitiveState primitive;
	primitive.p_a = 1.0;
	primitive.n_a = 1.0;
	primitive.u_a = Eigen::Vector3d::Zero();

	primitive.p_b = 1.0;
	primitive.n_b = 1.0;
	primitive.u_b = Eigen::Vector3d::Zero();

	return getConservationState(primitive);
}

const Eigen::VectorXd TwoFluidSolver::getInitStateRight() const
{
	PrimitiveState primitive;
	primitive.p_a = 0.1;
	primitive.n_a = 0.125;
	primitive.u_a = Eigen::Vector3d::Zero();

	primitive.p_b = 0.1;
	primitive.n_b = 0.125;
	primitive.u_b = Eigen::Vector3d::Zero();

	return getConservationState(primitive);
}

const Eigen::VectorXd TwoFluidSolver::getState1d(const Eigen::VectorXd & q_3d, const Eigen::Vector3d & n)
{
	Eigen::VectorXd q_1d(5);
	q_1d(0) = q_3d(0);
	q_1d(1) = q_3d.segment(1, 3).dot(n);
	q_1d(2) = q_3d(4);
	return q_1d;
}

std::ostream & operator<<(std::ostream & os, const TwoFluidSolver::PrimitiveState & obj)
{
	os << "p_a:\t" << obj.p_a << std::endl;
	os << "n_a:\t" << obj.n_a << std::endl;
	os << "u_a:\t" << obj.u_a.transpose() << std::endl;
	os << "p_b:\t" << obj.p_b << std::endl;
	os << "n_b:\t" << obj.n_b << std::endl;
	os << "u_b:\t" << obj.u_b << std::endl;
	return os;
}
