#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include "H5Writer.h"
#include "H5Reader.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "EigenAuxiliaries.h"


#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "Cell.h"
#include "TriPrism.h"
#include "XmlElement.h"
#include "XdmfAttribute.h"
#include "XdmfDataItem.h"
#include "XdmfDomain.h"
#include "XdmfGeometry.h"
#include "XdmfGrid.h"
#include "XdmfRoot.h"
#include "XdmfTime.h"
#include "XdmfTopology.h"


class Mesh
{
public:
	struct MeshInfo {
		int nVertices = 0;         // number of vertices
		int nVerticesBoundary = 0; // number of vertices on domain boundary
		int nVerticesTerminal = 0; // number of degrees of freedom with Dirichlet conditions

		int nEdges = 0;      // number of edges
		int nEdgesInner = 0; // number of edges in interior of domain

		int nFaces = 0;      // number of faces
		int nFacesInner = 0; // number of faces in interior of domain

		int nCells = 0; // number of cells
	};



	Mesh();
	Mesh(const std::string & meshPrefix);
	~Mesh();

	void writeToFile();
	void writeXdmf();

	Vertex * addVertex(const Eigen::Vector3d & position);
	Edge * addEdge(Vertex * A, Vertex * B);
	Face * addFace(const std::vector<Edge*> & faceEdges);
	Face * addFace(const std::vector<Vertex*> & faceVertices);

	TriPrism * addTriPrism(const std::vector<Face*> & sideFaces, Face * bottomFace, Face * topFace);
	Cell * addCell(const std::vector<Face*> & cellFaces);
	Cell * getCell(const std::vector<Face*> & cellFaces);

	Vertex * getVertex(const int index) const;
	Edge * getEdge(const int index) const;
	Face * getFace(const int index) const;
	Face * getFace(const std::vector<Edge*> & faceEdges);
	Cell * getCell(const int index) const;

	const std::string getPrefix() const;

	void createIncidenceMaps();

	const int getNumberOfVertices() const;
	const int getNumberOfEdges() const;
	const int getNumberOfFaces() const;
	const int getNumberOfCells() const;

	const std::vector<int> getXdmfTopology_cell2vertexIndices() const;
	const std::vector<int> getXdmfTopology_face2vertexIndices() const;

	const std::vector<Cell*> getCells() const;
	const std::vector<Face*> getFaces() const;

	void check() const;

	const Eigen::SparseMatrix<int> & get_f2eMap() const;
	const Eigen::SparseMatrix<int> & get_e2vMap() const;

	const Eigen::VectorXi getVertexTypes() const;
	const Eigen::VectorXi getEdgeTypes() const;
	const Eigen::VectorXi getFaceTypes() const;

	MeshInfo getMeshInfo() const;


protected:
	std::vector<Vertex*> vertexList;
	std::vector<Edge*> edgeList;
	std::vector<Face*> faceList;
	std::vector<Cell*> cellList;
	Eigen::Matrix3Xd vertexCoordinates;

	std::vector<Edge*> makeContinuousLoop(std::vector<Edge*> edges);

	Eigen::SparseMatrix<int> edge2vertexMap;
	Eigen::SparseMatrix<int> face2edgeMap;
	Eigen::SparseMatrix<int> cell2faceMap;

	virtual XdmfGrid getXdmfSurfaceGrid() const;
	virtual XdmfGrid getXdmfVolumeGrid() const;


private:
	std::string meshPrefix = "mesh";



	bool isConnected(const Vertex * A, const Vertex * B) const;
	bool isConnected(const std::vector<Edge*> & faceEdges) const;
	bool isConnected(const std::vector<Face*> & cellFaces) const;
	Edge * getEdge(Vertex * A, Vertex * B);

	void create_vertexCoordinates();
	void create_edge2vertex_map();
	void create_face2edge_map();
	void create_cell2face_map();

	XdmfGrid getXdmfVertexGrid() const;
	XdmfGrid getXdmfEdgeGrid() const;
	
	void writeXdmfVolumeMesh() const;

	//void writeXdmf_surface();
	//void writeXdmf_volume();

};

