#pragma once

#include "PrimalMesh.h"
#include "DualMesh.h"
#include "FluidState.h"
#include "Numerics.h"

#include <Eigen/SparseLU>


class AppmSolver
{
public:
	AppmSolver();
	~AppmSolver();

	void run();


private:
	PrimalMesh primalMesh;
	DualMesh dualMesh;

	Eigen::VectorXd B_h, E_h, J_h;
	Eigen::MatrixXd fluidStates;
	Eigen::MatrixXd fluidFluxes;

	Eigen::VectorXd x_m, x_mm1;

	Eigen::SparseMatrix<double> A;
	Eigen::SparseMatrix<double> C;
	Eigen::SparseMatrix<double> M;
	Eigen::SparseMatrix<double> M_d, M_f;
	Eigen::SparseMatrix<double> Q;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> maxwellSolver;

	std::vector<double> timeStamps;

	void init_meshes();
	void init_fluid();
	void init_maxwell(const double dt);


	const double update_fluid();
	void update_maxwell(const double dt, const double time);
	
	void writeXdmf();
	void writeXdmfDualVolume();

	void writeOutput(const int iteration, const double time);

	XdmfGrid getOutputPrimalEdgeGrid(const int iteration, const double time, const std::string & dataFilename);
	XdmfGrid getOutputPrimalSurfaceGrid(const int iteration, const double time, const std::string & dataFilename);
	
	XdmfGrid getOutputDualEdgeGrid(const int iteration, const double time, const std::string & dataFilename);
	XdmfGrid getOutputDualSurfaceGrid(const int iteration, const double time, const std::string & dataFilename);
	XdmfGrid getOutputDualVolumeGrid(const int iteration, const double time, const std::string & dataFilename);

	Eigen::SparseMatrix<int> setupOperatorQ();
	Eigen::SparseMatrix<double> setupOperatorMeps();
	Eigen::SparseMatrix<double> setupOperatorMnu();
};

