#include "Edge.h"



Edge::Edge()
{
}

Edge::Edge(const int index)
	: GeometryItem(index)
{
}

Edge::Edge(Vertex * A, Vertex * B)
{
	assert(A != nullptr);
	assert(B != nullptr);
	this->A = A;
	this->B = B;
	A->setAdjacient(this);
	B->setAdjacient(this);
}


Edge::~Edge()
{
}

Vertex * Edge::getVertexA()
{
	return A;
}

Vertex * Edge::getVertexB()
{
	return B;
}

bool Edge::isAdjacient(const Vertex * A, const Vertex * B) const
{
	assert(A != nullptr);
	assert(B != nullptr);
	assert(A != B);
	return (A == this->A && B == this->B) || (A == this->B && B == this->A);
}

bool Edge::isAdjacient(const Vertex * A, const Vertex * B) 
{
	assert(A != nullptr);
	assert(B != nullptr);
	assert(A != B);
	return (A == this->A && B == this->B) || (A == this->B && B == this->A);
}

void Edge::setAdjacient(Face * face)
{
	assert(face != nullptr);
	this->faceList.push_back(face);
}

bool Edge::hasConnectedFace(const std::vector<Edge*> & faceEdges) const
{
	for (auto face : faceList) {
		if (face->hasFaceEdges(faceEdges)) {
			return true;
		}
	}
	return false;
}

Face * Edge::getConnectedFace(const std::vector<Edge*>& faceEdges) const
{
	assert(hasConnectedFace(faceEdges));
	for (auto face : faceList) {
		if (face->hasFaceEdges(faceEdges)) {
			return face;
		}
	}
	assert(false);
	return nullptr;
}

Vertex * Edge::getOppositeVertex(const Vertex * v)
{
	assert(v != nullptr);
	assert(hasVertex(v));
	return (v == A) ? B : A;
}

bool Edge::hasVertex(const Vertex * v)
{
	assert(v != nullptr);
	return v == A || v == B;
}

std::ostream & operator<<(std::ostream & os, const Edge & obj)
{
	os << "Edge " << obj.getIndex() << ": ";
	os << "{" << obj.A->getIndex() << ", ";
	os << obj.B->getIndex() << "}";
	return os;
}
