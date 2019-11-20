#include "DualMesh.h"



DualMesh::DualMesh()
	: Mesh("dual")
{
}

DualMesh::DualMesh(const std::string & meshPrefix) 
	: Mesh(meshPrefix)
{
}

DualMesh::DualMesh(const PrimalMesh & primal)
	: Mesh("dual")
{

}


DualMesh::~DualMesh()
{
}
