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

protected:

	XdmfGrid getXdmfSurfaceGrid() const;
	XdmfGrid getXdmfVolumeGrid() const;


private:
	Eigen::VectorXi dualFaceToPrimalEdgeList;
	Eigen::VectorXi associateDualFacesWithPrimalEdges(const PrimalMesh & primal);

	void init_cellFluidType();
	void init_faceFluidType(const double terminalRadius);
};

