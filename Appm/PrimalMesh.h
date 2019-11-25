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

	void init();

private:
	void init_hexagon();
	void init_triangle();
	void refineMesh(const int nRefinements);
	void outerMeshExtrude(const int nLayers);
	void outerMeshExtrude();
	void extrudeMesh(const int nLayers);

	Eigen::Matrix3Xi refine_triangles();
	Eigen::Matrix3Xi refine_triangles_specialCorners();

	void test_quadFace();
};

