#pragma once

#include "PrimalMesh.h"
#include "DualMesh.h"
#include "FluidState.h"

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

	Eigen::VectorXd bvec, dvec, E_h, hvec, jvec;
	Eigen::MatrixXd fluidStates;
	Eigen::MatrixXd fluidFluxes;

	std::vector<double> timeStamps;

	void init_meshes();
	void init_fluid();

	const double update_fluid();

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

