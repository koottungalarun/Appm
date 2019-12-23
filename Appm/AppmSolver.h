#pragma once

#include "PrimalMesh.h"
#include "DualMesh.h"
#include "FluidState.h"
#include "Numerics.h"

#include <Eigen/SparseLU>


#define _RT_ONECELL 
#undef _RT_ONECELL


class AppmSolver
{

public:
	AppmSolver();
	//AppmSolver(const AppmSolver & other);
	virtual ~AppmSolver();

	void run();

protected:
	struct MeshInfo {
		int nVertices = 0;         // number of vertices
		int nVerticesBoundary = 0; // number of vertices on domain boundary
		int nVerticesTerminal = 0; // number of degrees of freedom with Dirichlet conditions

		int nEdges = 0;      // number of edges
		int nEdgesInner = 0; // number of edges in interior of domain

		int nFaces = 0;      // number of faces
		int nFacesInner = 0; // number of faces in interior of domain

		int nCells = 0; // number of cells
	} primalMeshInfo, dualMeshInfo;

	PrimalMesh primalMesh;
	DualMesh dualMesh;

	Eigen::VectorXd B_h, E_h, H_h, J_h;
	Eigen::Matrix3Xd B_vertex;

	Eigen::MatrixXd fluidStates;
	Eigen::MatrixXd fluidFluxes;

	bool isMaxwellCurrentSource = false;

	virtual void init_maxwell(const double dt) = 0;
	virtual void init_maxwell() = 0;
	virtual void update_maxwell(const double dt, const double time) = 0;

	Eigen::SparseMatrix<int> setupOperatorQ();
	Eigen::SparseMatrix<double> setupOperatorMeps();
	Eigen::SparseMatrix<double> setupOperatorMnu();

	Eigen::VectorXd electricPotentialTerminals(const double time);

	Eigen::SparseMatrix<double> speye(const int rows, const int cols);

	Eigen::SparseMatrix<double> hodgeOperatorPrimalEdgeToDualFace();
	Eigen::SparseMatrix<double> hodgeOperatorDualEdgeToPrimalFace();
	Eigen::SparseMatrix<double> hodgeOperatorElectricalConductivity();

	/** Inclusion operator of boundary vertices into all vertices */
	Eigen::SparseMatrix<double> inclusionOperatorBoundaryVerticesToAllVertices();


	AppmSolver::MeshInfo setMeshInfo(const Mesh & mesh);

private:
	void interpolateMagneticFluxToPrimalVertices();

	// void test_raviartThomas();
	const Eigen::Matrix3Xd getPrismReferenceCoords(const int nSamples);

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

	void setAzimuthalMagneticFluxField();
	void setUniformMagneticFluxField(const Eigen::Vector3d & fieldVector);
	void init_RaviartThomasInterpolation();

	std::vector<Eigen::Matrix3d> rt_piolaMatrix;
	std::vector<Eigen::Vector3d> rt_piolaVector;

	void setTorusCurrent();

};

