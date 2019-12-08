#pragma once
#include "GeometryItem.h"
#include "Vertex.h"
class Vertex;

#include "Face.h"
class Face;

class Edge :
	public GeometryItem
{
public:

	enum class Type {
		Interior, InteriorToBoundary, Boundary
	};

	Edge();
	Edge(const int index);
	Edge(Vertex * A, Vertex * B);
	~Edge();

	Vertex * getVertexA();
	const Vertex * getVertexA() const;
	Vertex * getVertexB();
	const Vertex * getVertexB() const;
	const Eigen::Vector3d getDirection() const;
	const double getLength() const;

	bool isAdjacient(const Vertex * A, const Vertex * B) const;
	bool isAdjacient(const Vertex * A, const Vertex * B);

	void setAdjacient(Face * face);

	bool hasConnectedFace(const std::vector<Edge*> & faceEdges) const;
	Face * getConnectedFace(const std::vector<Edge*> & faceEdges) const;

	friend std::ostream & operator<<(std::ostream & os, const Edge & obj);

	Vertex * getOppositeVertex(const Vertex * v) const;
	bool hasVertex(const Vertex * v) const;
	bool isBoundary() const;

	const Eigen::Vector3d getHalfwayPosition() const;

	Vertex * getCoincidentVertex(const Edge * other) const;
	bool hasCoincidentVertex(const Edge * other) const;

	int getIncidence(const Vertex * v) const;
	const std::vector<Face*> getFaceList() const;

	void setType(const Edge::Type & type);
	const Edge::Type getType() const;


private:
	Vertex * A = nullptr;
	Vertex * B = nullptr;
	Eigen::Vector3d edgeCenter;
	std::vector<Face*> faceList;
	Type type;
};

