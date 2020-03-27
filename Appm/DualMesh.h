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

protected:

	XdmfGrid getXdmfSurfaceGrid() const;
	XdmfGrid getXdmfVolumeGrid() const;


private:

	void init_cellFluidType();
	void init_faceFluidType(const double terminalRadius);
};

