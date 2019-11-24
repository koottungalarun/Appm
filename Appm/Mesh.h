#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

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

class Mesh
{
public:
	Mesh();
	Mesh(const std::string & meshPrefix);
	~Mesh();

	void writeToFile();
	void writeXdmf();

	Vertex * addVertex(const Eigen::Vector3d & position);
	Edge * addEdge(Vertex * A, Vertex * B);
	Face * addFace(const std::vector<Edge*> & faceEdges);
	//Face * addFace(const std::vector<Vertex*> & faceVertices);

	TriPrism * addTriPrism(const std::vector<Face*> & sideFaces, Face * bottomFace, Face * topFace);
	Cell * addCell(const std::vector<Face*> & cellFaces);
	Cell * getCell(const std::vector<Face*> & cellFaces);

	Vertex * getVertex(const int index);
	Edge * getEdge(const int index);
	Face * getFace(const int index);
	Face * getFace(const std::vector<Edge*> & faceEdges);
	Cell * getCell(const int index) const;

	void createIncidenceMaps();

	const int getNumberOfVertices() const;

	const std::vector<Cell*> getCells() const;
	const std::vector<Face*> getFaces() const;

	void check() const;

protected:
	std::vector<Vertex*> vertexList;
	std::vector<Edge*> edgeList;
	std::vector<Face*> faceList;
	std::vector<Cell*> cellList;
	Eigen::Matrix3Xd vertexCoordinates;


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

	Eigen::SparseMatrix<int> edge2vertexMap;
	Eigen::SparseMatrix<int> face2edgeMap;
	Eigen::SparseMatrix<int> cell2faceMap;

	void writeXdmf_surface();
	void writeXdmf_volume();

};

