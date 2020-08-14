#pragma once

#include "EigenAuxiliaries.h"
#include "Mesh.h"
#include "PrimalMesh.h"

class DualMesh :
	public Mesh
{
public:
	DualMesh();
	DualMesh(const std::string & meshPrefix);
	~DualMesh();

	void init_dualMesh(const PrimalMesh & primal, const double terminalRadius);
	const int getAssociatedPrimalEdgeIndex(const int dualFaceIndex) const;

	const int getNumberFluidCells() const;

protected:

	XdmfGrid getXdmfSurfaceGrid() const;
	XdmfGrid getXdmfVolumeGrid() const;


private:
	Eigen::SparseVector<int> primalFaceToDualVertex;
	Eigen::SparseVector<int> primalEdgeToDualVertex;
	Eigen::SparseVector<int> primalVertexToDualVertex;


	Eigen::VectorXi dualFaceToPrimalEdgeList;
	Eigen::VectorXi associateDualFacesWithPrimalEdges(const PrimalMesh & primal);

	void init_cellFluidType(const PrimalMesh & primal);
	void init_faceFluidType(const double terminalRadius);

	void init_dualVertices(const PrimalMesh & primal);
	void init_dualEdges(const PrimalMesh & primal);
	void init_dualFaces(const PrimalMesh & primal);
	void init_dualCells(const PrimalMesh & primal);
};

