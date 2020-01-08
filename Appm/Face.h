#pragma once
#include "GeometryItem.h"
#include "Edge.h"
class Edge;
#include "Cell.h"
class Cell;

class Face :
	public GeometryItem
{
public:
	Face();
	Face(const std::vector<Edge*> & faceEdges);
	Face(const std::vector<Vertex*> & faceVertices);
	~Face();

	std::vector<Edge*> getEdgeList() const;
	std::vector<Vertex*> getVertexList() const;
	std::vector<Cell*> getCellList() const;

	bool hasFaceEdges(const std::vector<Edge*> faceEdges) const;
	

	friend std::ostream & operator<<(std::ostream & os, const Face & obj);

	const int getOrientation(const Edge * edge);

	bool isBoundary() const;
	bool hasBoundaryEdges() const;

	const Eigen::Vector3d getCenter() const;
	const double getArea() const;

	bool hasCommonCell(const Face * other) const;
	Cell * getCommonCell(const Face * other) const;

	void setAdjacient(Cell * cell);

	const Eigen::Vector3d getNormal() const;
	void setNormal(const Eigen::Vector3d & fn);

private:
	std::vector<Edge*> edgeList;
	std::vector<Vertex*> vertexList;
	std::vector<Cell*> cellList;

	Eigen::Vector3d center;
	Eigen::Vector3d faceNormal;
	double area = 0;

	void init();
	bool isListOfVerticesUnique() const;
	bool isListOfEdgesUnique() const;

	const Eigen::Vector3d getCircumCenter() const;
	const double computeArea() const;
};

