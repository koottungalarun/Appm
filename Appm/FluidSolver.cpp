#include "FluidSolver.h"



FluidSolver::FluidSolver()
{
}


FluidSolver::FluidSolver(const DualMesh * mesh)
{
	this->mesh = mesh;
	init();
}

FluidSolver::~FluidSolver()
{
}

const double FluidSolver::updateFluidState()
{
	std::cout << "You should call the inherited function instead of this one" << std::endl;
	return 0.0;
}

void FluidSolver::writeStates(H5Writer & writer) const
{
	const int nDualCells = mesh->getNumberOfCells();

	Eigen::VectorXd q_density = fluidStates.row(0);
	Eigen::MatrixXd q_momentum = fluidStates.middleRows(1, 3);
	Eigen::VectorXd q_energy = fluidStates.row(4);
	writer.writeData(q_density, "/qDensity");
	writer.writeData(q_momentum, "/qMomentum");
	writer.writeData(q_energy, "/qEnergy");

	// Fluid states in primitive variables
	Eigen::VectorXd pressure(nDualCells);
	Eigen::VectorXd density(nDualCells);
	Eigen::Matrix3Xd velocity(3, nDualCells);
	for (int i = 0; i < nDualCells; i++) {
		const FluidPrimitiveState primitiveState = FluidState(fluidStates.col(i)).getPrimitiveState();
		pressure(i) = primitiveState.p;
		density(i) = primitiveState.rho;
		velocity.col(i) = primitiveState.u;
	}
	writer.writeData(pressure, "/pressure");
	writer.writeData(density, "/density");
	writer.writeData(velocity, "/velocity");
}

void FluidSolver::init()
{
	const int nCells = mesh->getNumberOfCells();
	assert(nCells > 0);

	fluidStates = Eigen::MatrixXd::Zero(5, nCells);
	fluidFluxes = Eigen::MatrixXd::Zero(5, nCells);

	const double rho_L = 1.0;
	const double p_L = 1.0;
	const double rho_R = 0.125;
	const double p_R = 0.1;
	const Eigen::Vector3d uZero = Eigen::Vector3d::Zero();

	FluidPrimitiveState primitiveL;
	primitiveL.rho = rho_L;
	primitiveL.p = p_L;
	primitiveL.u = uZero;

	FluidPrimitiveState primitiveR;
	primitiveR.rho = rho_R;
	primitiveR.p = p_R;
	primitiveR.u = uZero;

	const FluidState qL(primitiveL);
	const FluidState qR(primitiveR);

	const Eigen::VectorXd qL_vec = qL.getVecState();
	const Eigen::VectorXd qR_vec = qR.getVecState();

	std::cout << "qL: " << qL_vec.transpose() << std::endl;
	std::cout << "qR: " << qR_vec.transpose() << std::endl;
	std::cout << std::endl;

	for (int i = 0; i < nCells; i++) {
		const Cell * cell = mesh->getCell(i);
		const double idxC = cell->getIndex();
		const Eigen::Vector3d cellCenter = cell->getCenter();

		fluidStates.col(i) = (cellCenter(2) < 0.5) ? qL_vec : qR_vec;
	}

	std::cout << "Test for fluid flux: " << std::endl;
	double dt_loc = 0;
	Eigen::VectorXd fluidFlux = Numerics::fluidFlux_rusanov(qL_vec, qR_vec, Eigen::Vector3d(0, 0, 1), 1, dt_loc);
	std::cout << "Flux: " << fluidFlux.transpose() << std::endl;
	std::cout << "dt_loc: " << dt_loc << std::endl;
	
}
