#pragma once

#include "GeometryItem.h"
#include "Edge.h"
class Edge;

class Vertex
	: public GeometryItem
{
public:
	Vertex();
	Vertex(const int index);
	Vertex(const Eigen::Vector3d & position);
	Vertex(const Eigen::Vector3d & position, const int index);
	~Vertex();

	const Eigen::Vector3d getPosition() const;

	void setAdjacient(Edge * edge);
	bool isAdjacientTo(const Vertex * other) const;
	Edge * getAdjacientEdge(const Vertex * other);

	std::vector<Edge*> getEdges();

	friend std::ostream & operator<<(std::ostream & os, const Vertex & obj);

private:
	Eigen::Vector3d position;
	std::vector<Edge*> adjacientEdges;
};

