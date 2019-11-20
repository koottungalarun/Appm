#include "Face.h"



Face::Face()
{
}

Face::Face(std::vector<Edge*> faceEdges)
{
	this->edgeList = faceEdges;
	for (auto edge : edgeList) {
		edge->setAdjacient(this);
	}
}


Face::~Face()
{
}

std::vector<Edge*> Face::getEdgeList() const
{
	return this->edgeList;
}

std::vector<Vertex*> Face::getVertexList() const
{
	std::vector<Vertex*> faceVertices;

	Edge * e = edgeList[0];
	Vertex * v = e->getVertexB();

	const int nEdges = edgeList.size();
	for (int i = 0; i < nEdges; i++) {
		faceVertices.push_back(v);
		Vertex * otherVertex = e->getOppositeVertex(v);
		assert(otherVertex != nullptr);
		Edge * otherEdge = nullptr;
		for (auto edge : edgeList) {
			if (edge != e && edge->hasVertex(otherVertex)) {
				otherEdge = edge;
			}
		}
		assert(otherEdge != nullptr);
		e = otherEdge;
		v = otherVertex;
	}

	assert(faceVertices.size() == edgeList.size());
	return faceVertices;
}

bool Face::hasFaceEdges(const std::vector<Edge*> faceEdges) const
{
	std::vector<int>::iterator it;
	std::vector<int> edgeIdx;
	for (auto edge : faceEdges) {
		edgeIdx.push_back(edge->getIndex());
	}
	std::sort(edgeIdx.begin(), edgeIdx.end());
	// remove duplicate edge indices 
	it = std::unique(edgeIdx.begin(), edgeIdx.end());
	edgeIdx.resize( std::distance(edgeIdx.begin(), it) );

	if (faceEdges.size() != edgeList.size()) {
		return false;
	}

	std::vector<int> thisFaceEdgeIdx;
	for (auto edge : edgeList) {
		thisFaceEdgeIdx.push_back(edge->getIndex());
	}
	std::sort(thisFaceEdgeIdx.begin(), thisFaceEdgeIdx.end());
	
	std::vector<int> result(edgeIdx.size());
	it = std::set_intersection(edgeIdx.begin(), edgeIdx.end(), thisFaceEdgeIdx.begin(), thisFaceEdgeIdx.end(), result.begin());
	const int distance = std::distance(result.begin(), it);
	if (distance < edgeList.size()) {
		return false;
	}

	return true;
}

const int Face::getOrientation(const Edge * edge)
{
	return 0;
}

std::ostream & operator<<(std::ostream & os, const Face & obj)
{
	os << "Face " << obj.getIndex() << ": {";
	for (auto edge : obj.edgeList) {
		os << edge->getIndex() << ",";
	}
	os << "}";
	return os;
}
