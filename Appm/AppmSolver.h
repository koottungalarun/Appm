#pragma once

#include "PrimalMesh.h"
#include "DualMesh.h"
#include "Numerics.h"
#include "FluidSolver.h"
#include "SingleFluidSolver.h"
#include "TwoFluidSolver.h"
#include "MultiFluidSolver.h"

#include "MaxwellSolver.h"
#include "MaxwellSolverCrankNicholson.h"
#include "MaxwellSolverImplicitEuler.h"

#include "Physics.h"

#include <Eigen/SparseLU>


#define _RT_ONECELL 
#undef _RT_ONECELL


class AppmSolver
{

public:
	AppmSolver();
	AppmSolver(const PrimalMesh::PrimalMeshParams & primalMeshParams);
	~AppmSolver();

	void run();

protected:
	PrimalMesh primalMesh;
	DualMesh dualMesh;

	Eigen::Matrix3Xd B_vertex;

	bool isMaxwellCurrentSource = false;

private:

	struct ParticleParameters {
		std::string name = "neutral";
		double mass = 1.0;
		int electricCharge = 0;
	};
	std::vector<ParticleParameters> particleParams;

	/** Get number of fluids. */
	const int getNFluids() const;

	bool isStateWrittenToOutput = false;



	bool isMaxwellEnabled = false;
	bool isFluidEnabled = true;
	double timestepSize = 1.0;
	int maxIterations = 0;
	double maxTime = 0;
	double lambdaSquare = 1.0;

	bool isWriteEfield = false;
	bool isWriteBfield = false;
	bool isWriteHfield = false;
	bool isWriteJfield = false;

	Eigen::MatrixXd fluidStates;
	Eigen::MatrixXd fluidStates_new;
	Eigen::MatrixXd fluidSources;
	Eigen::MatrixXd fluidFluxes;
	Eigen::MatrixXd faceFluxes;

	int faceIdxRef = -1;

	// Isentropic expansion coefficient, aka ratio of heat capacities
	//const double gamma = 1.4; 

	void init_multiFluid(const std::string & filename);

	void init_SodShockTube(const double zRef);
	void init_Uniformly(const double n, const double p, const double u);
	void init_Explosion();

	
	const int getFluidStateLength() const;
	const double getNextFluidTimestepSize() const;
	const double getWaveSpeed(const Eigen::VectorXd & state, const Eigen::Vector3d & fn) const;
	const double getWaveSpeed(const Eigen::Vector3d & state) const;

	const Eigen::VectorXd getFluidState(const int cellIdx, const int fluidIdx) const;
	const Eigen::Vector3d getFluidState(const int cellIdx, const int fluidIdx, const Eigen::Vector3d & faceNormal) const;
	
	const int getOrientation(const Cell * cell, const Face * face) const;
	//const bool isFaceCellsReversed(const int faceIdx) const;

	const Eigen::Vector3d getFluidFluxFromState(const Eigen::Vector3d & q) const;
	//const Eigen::Vector3d getRusanovFluxExplicit(const Eigen::Vector3d & qL, const Eigen::Vector3d & qR) const;
	const Eigen::Vector3d getRusanovFluxExplicit(const int faceIdx, const int fluidIdx) const;
	const double getMomentumUpdate(const int k, const Eigen::Vector3d & nvec, const int fluidIdx) const;
	
	const Eigen::Vector3d getFluidStateProjected(const Eigen::VectorXd & state, const Eigen::Vector3d & fn) const;
	void interpolateMagneticFluxToPrimalVertices();

	// void test_raviartThomas();
	const Eigen::Matrix3Xd getPrismReferenceCoords(const int nSamples);

	std::vector<double> timeStamps;

	void init_meshes(const PrimalMesh::PrimalMeshParams & primalParams);

	void writeXdmf();
	void writeXdmfDualVolume();
	void writeXdmfDualFaceFluxes();

	void writeOutput(const int iteration, const double time);

	void writeFluidStates(H5Writer & writer);
	void writeMaxwellStates(H5Writer & writer);

	XdmfGrid getOutputPrimalEdgeGrid(const int iteration, const double time, const std::string & dataFilename);
	XdmfGrid getOutputPrimalSurfaceGrid(const int iteration, const double time, const std::string & dataFilename);
	
	XdmfGrid getOutputDualEdgeGrid(const int iteration, const double time, const std::string & dataFilename);
	XdmfGrid getOutputDualSurfaceGrid(const int iteration, const double time, const std::string & dataFilename);
	XdmfGrid getOutputDualVolumeGrid(const int iteration, const double time, const std::string & dataFilename);

	void init_RaviartThomasInterpolation();

	std::vector<Eigen::Matrix3d> rt_piolaMatrix;
	std::vector<Eigen::Vector3d> rt_piolaVector;

	void readParameters(const std::string & filename);

	const std::string xdmf_GridPrimalEdges(const int iteration) const;
	const std::string xdmf_GridPrimalFaces(const int iteration) const;
	const std::string xdmf_GridDualEdges(const int iteration) const;
	const std::string xdmf_GridDualFaces(const int iteration) const;
	const std::string xdmf_GridDualCells(const int iteration) const;

	const std::string fluidXdmfOutput(const std::string & filename) const;
};

