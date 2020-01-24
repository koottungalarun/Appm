#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include "Mesh.h"
class PrimalMesh :
	public Mesh
{
public:
	class PrimalMeshParams {
	public:
		PrimalMeshParams();
		PrimalMeshParams(const std::string & filename);

		const int getRefinements() const;
		const int getAxialLayers() const;
		const int getOuterLayers() const;
		const double getElectrodeRadius() const;
		const double getZmax() const;

		friend std::ostream & operator<<(std::ostream & os, const PrimalMeshParams & obj);

	private:
		int nAxialLayers = 10;
		int nRefinements = 2;
		int nOuterLayers = 2;
		double electrodeRadius = 0.35;
		double zmax = 1;

		void readParameters(const std::string & filename);
	};

	PrimalMesh();
	PrimalMesh(const std::string & meshPrefix);
	PrimalMesh(const PrimalMeshParams & p);
	~PrimalMesh();

	void init();

private:
	// Mesh parameters
	PrimalMeshParams params;

	void init_hexagon(const double zValue);
	void init_triangle();
	void refineMesh(const int nRefinements);
	void outerMeshExtrude(const int nLayers);
	void outerMeshExtrude_triangles();
	void outerMeshExtrude_prisms();
	void extrudeMesh(const int nLayers, const double zmax);

	Eigen::Matrix3Xi refine_triangles();
	Eigen::Matrix3Xi refine_triangles_specialCorners();

	void test_quadFace();


	void sortVertices(const double electrodeRadius);
	void sortEdges();
	void sortFaces();
	void sortCells();

	void validateParameters();
	void check_zCoord(const double z0);
};

