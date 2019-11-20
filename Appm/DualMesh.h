#pragma once

#include "Mesh.h"
#include "PrimalMesh.h"

class DualMesh :
	public Mesh
{
public:
	DualMesh();
	DualMesh(const std::string & meshPrefix);
	DualMesh(const PrimalMesh & primal);
	~DualMesh();
};

