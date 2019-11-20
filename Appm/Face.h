#pragma once
#include "GeometryItem.h"
#include "Edge.h"
class Edge;

class Face :
	public GeometryItem
{
public:
	Face();
	Face(std::vector<Edge*> faceEdges);
	~Face();

	std::vector<Edge*> getEdgeList() const;
	std::vector<Vertex*> getVertexList() const;

	bool hasFaceEdges(const std::vector<Edge*> faceEdges) const;

	friend std::ostream & operator<<(std::ostream & os, const Face & obj);

	const int getOrientation(const Edge * edge);

private:
	std::vector<Edge*> edgeList;
};

