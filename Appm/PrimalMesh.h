#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include "Mesh.h"
class PrimalMesh :
	public Mesh
{
public:
	struct PrimalMeshParams {
		int nRefinements = 0;
		int nAxialLayers = 1;
		int nOuterLayers = 0;
		double electrodeRadius = 0.35;
	};

	PrimalMesh();
	PrimalMesh(const std::string & meshPrefix);
	PrimalMesh(const PrimalMeshParams & p);
	~PrimalMesh();

	void init();

private:
	// Mesh parameters
	PrimalMeshParams params;

	void init_hexagon();
	void init_triangle();
	void refineMesh(const int nRefinements);
	void outerMeshExtrude(const int nLayers);
	void outerMeshExtrude();
	void extrudeMesh(const int nLayers, const double zmax);

	Eigen::Matrix3Xi refine_triangles();
	Eigen::Matrix3Xi refine_triangles_specialCorners();

	void test_quadFace();


	void sortVertices(const double electrodeRadius);
	void sortEdges();
	void sortFaces();

	void validateParameters();
};

