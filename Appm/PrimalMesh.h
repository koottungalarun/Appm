#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include "Mesh.h"
class PrimalMesh :
	public Mesh
{
public:
	PrimalMesh();
	PrimalMesh(const std::string & meshPrefix);
	~PrimalMesh();

private:
	void init();

};

