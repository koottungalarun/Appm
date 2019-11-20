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

	Vertex * getVertex(const int index);
	Face * getFace(const std::vector<Edge*> & faceEdges);

	void createIncidenceMaps();


private:
	std::string meshPrefix = "mesh";

	std::vector<Vertex*> vertexList;
	std::vector<Edge*> edgeList;
	std::vector<Face*> faceList;

	Eigen::Matrix3Xd vertexCoordinates;

	bool isConnected(const Vertex * A, const Vertex * B) const;
	bool isConnected(const std::vector<Edge*> & faceEdges);
	Edge * getEdge(Vertex * A, Vertex * B) const;

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

