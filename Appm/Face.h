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

	enum class Type {
		DEFAULT, INTERIOR, OPENING, TERMINAL, WALL
	};

	friend std::ostream & operator<<(std::ostream & os, const Face::Type & obj);
	//enum class EmagType {
	//	DEFAULT
	//};

	Face();
	Face(const std::vector<Edge*> & faceEdges);
	Face(const std::vector<Vertex*> & faceVertices);
	~Face();

	std::vector<Edge*> getEdgeList() const;
	std::vector<Vertex*> getVertexList() const;
	std::vector<Cell*> getCellList() const;

	bool hasFaceEdges(const std::vector<Edge*> faceEdges) const;
	bool hasFluidCells() const;
	bool isFluidBoundary() const;
	

	friend std::ostream & operator<<(std::ostream & os, const Face & obj);

	const int getOrientation(const Edge * edge);

	const bool isAdjacient(const Cell * cell) const;

	bool isBoundary() const;
	bool hasBoundaryEdges() const;

	const Eigen::Vector3d getCenter() const;
	const double getArea() const;

	bool hasCommonCell(const Face * other) const;
	Cell * getCommonCell(const Face * other) const;

	void setAdjacient(Cell * cell);

	const Eigen::Vector3d getNormal() const;
	//void setNormal(const Eigen::Vector3d & fn);
	void flipNormal();

	void setType(const Type & fluidType);
	const Type getType() const;


private:
	std::vector<Edge*> edgeList;
	std::vector<Vertex*> vertexList;
	std::vector<Cell*> cellList;

	Eigen::Vector3d center;
	Eigen::Vector3d faceNormal;
	double area = 0;

	Type fluidType = Type::DEFAULT;

	void init();
	bool isListOfVerticesUnique() const;
	bool isListOfEdgesUnique() const;

	const Eigen::Vector3d getCircumCenter() const;
	const double computeArea() const;

	const Eigen::Vector3d computeFaceNormal(const int version) const;
	const Eigen::Vector3d computeFaceNormal_version1() const;
	const Eigen::Vector3d computeFaceNormal_version2() const;


};

