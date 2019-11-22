#include "PrimalMesh.h"



PrimalMesh::PrimalMesh() :
	PrimalMesh("primal")
{
	std::cout << "Call to PrimalMesh()" << std::endl;
}

PrimalMesh::PrimalMesh(const std::string & meshPrefix)
	: Mesh(meshPrefix)
{
	std::cout << "Call to PrimalMesh(string)" << std::endl;
}


PrimalMesh::~PrimalMesh()
{
	std::cout << "Call to ~PrimalMesh()" << std::endl;
}

void PrimalMesh::init()
{
	init_hexagon();
	const int nRefinements = 0;
	const int nOuterMeshLayers = 0;
	if (nOuterMeshLayers > 0) {
		assert(nRefinements > 1);
	}
	refineMesh(nRefinements);
	outerMeshExtrude(nOuterMeshLayers);

	const int nLayers = 1;
	extrudeMesh(Eigen::Vector3d(0, 0, 1), nLayers);
}


void PrimalMesh::init_hexagon()
{
	std::cout << "Initialize with hexagon" << std::endl;
	Vertex * origin = addVertex(Eigen::Vector3d(0, 0, 0));
	const int corners = 6;
	for (int k = 0; k < corners; k++) {
		const double phi = 2 * M_PI * k / corners;
		const Eigen::Vector3d pos(cos(phi), sin(phi), 0);
		Vertex * v = addVertex(pos);
		addEdge(origin, v);
	}

	for (int k = 0; k < corners; k++) {
		Vertex * A = getVertex(k + 1);
		Vertex * B = getVertex(((k + 1) % corners) + 1);
		auto edgeA = addEdge(origin, A);
		auto edgeB = addEdge(origin, B);
		auto edgeC = addEdge(A, B);
		std::vector<Edge*> faceEdges = { edgeA, edgeB, edgeC };
		Face * face = addFace(faceEdges);
	}
}

void PrimalMesh::refineMesh(const int nRefinements)
{
	// Strategy: 
	// construct refined f2v-map from current mesh 
	// and reconstruct from that map

	for (int level = 0; level < nRefinements; level++) {
		std::cout << "Mesh refinement level: " << level << std::endl;
		Eigen::Matrix3Xi f2v;

		// Get face-to-vertex list of refined faces
		if (level != 1) {
			f2v = refine_triangles();
		}
		else {
			f2v = refine_triangles_specialCorners();
		}
		//std::cout << "f2v: " << std::endl;
		//std::cout << f2v.transpose() << std::endl;
		assert((f2v.array() >= 0).all());

		std::ofstream file;
		std::stringstream ss;
		ss << "level" << level << "-coords.dat";
		file = std::ofstream(ss.str());
		file << vertexCoordinates.transpose() << std::endl;
		ss = std::stringstream();
		ss << "level" << level << "-f2v.dat";
		file = std::ofstream(ss.str());
		file << f2v.transpose() << std::endl;

		//if (level == 1) { break; }

		// clear mesh elements
		for (auto v : vertexList) {
			delete v;
		}
		vertexList.resize(0);
		for (auto e : edgeList) {
			delete e;
		}
		edgeList.resize(0);
		for (auto f : faceList) {
			delete f;
		}
		faceList.resize(0);

		// reconstruct mesh
		int nVertices = vertexCoordinates.cols();
		assert(vertexCoordinates.rows() == 3);
		for (int i = 0; i < nVertices; i++) {
			Eigen::Vector3d pos = vertexCoordinates.col(i);
			addVertex(pos);
		}
		assert(f2v.rows() == 3);
		int nFaces = f2v.cols();
		for (int i = 0; i < nFaces; i++) {
			Eigen::Vector3i triangleIdx = f2v.col(i);
			std::vector<Edge*> faceEdges;
			for (int j = 0; j < 3; j++) {
				assert((triangleIdx.array() >= 0).all());
				const int idxA = triangleIdx(j);
				const int idxB = triangleIdx((j + 1) % 3);

				Vertex * A = getVertex(idxA);
				Vertex * B = getVertex(idxB);
				Edge * edge = addEdge(A, B);
				faceEdges.push_back(edge);
			}
			Face * face = addFace(faceEdges);
		}

		// Check
		for (auto face : faceList) {
			assert(face->getEdgeList().size() == 3);
			assert(face->getVertexList().size() == 3);
		}
	}

}

void PrimalMesh::outerMeshExtrude(const int nLayers)
{
	for (int layer = 0; layer < nLayers; layer++) {
		outerMeshExtrude();
	}
}

void PrimalMesh::outerMeshExtrude()
{
	const int nVertices = vertexList.size();
	assert(nVertices == vertexCoordinates.cols());

	const int nEdges = edgeList.size();

	Eigen::SparseVector<int> boundaryVertices(nVertices);
	Eigen::SparseVector<int> boundaryEdges(nEdges);
	Eigen::SparseVector<int> boundaryEdgeVertexIdx(nEdges);

	// get boundary vertices and boundary edges
	for (int i = 0; i < nEdges; i++) {
		Edge * edge = edgeList[i];
		if (edge->isBoundary()) {
			boundaryEdges.coeffRef(i) = 1;
			const int idxA = edge->getVertexA()->getIndex();
			const int idxB = edge->getVertexB()->getIndex();
			boundaryVertices.coeffRef(idxA) = 1;
			boundaryVertices.coeffRef(idxB) = 1;
		}
	}

	// create outer vertices for each boundary edge
	const int nnzE = boundaryEdges.nonZeros();
	Eigen::Matrix3Xd boundaryEdgeVertexCoord(3, nnzE);
	const int offset_E = vertexCoordinates.cols();
	int idx = 0;
	for (int i = 0; i < nnzE; i++) {
		const int idxE = boundaryEdges.innerIndexPtr()[i];
		const Edge * edge = edgeList[idxE];
		assert(edge != nullptr);
		//const std::vector<Face*> edgeFaces = edge->getFaceList();
		//assert(edgeFaces.size() == 1);
		//const Face * face = edgeFaces[0];
		//assert(face != nullptr);
		//const std::vector<Vertex*> faceVertices = face->getVertexList();
		//assert(faceVertices.size() == 3);
		//const Vertex * A = edge->getVertexA();
		//const Vertex * B = edge->getVertexB();
		//Vertex * C = nullptr;
		//for (int j = 0; j < 3; j++) {
		//	if (faceVertices[j] != A && faceVertices[j] != B) {
		//		C = faceVertices[j];
		//		break;
		//	}
		//}
		//assert(C != nullptr);
		//const Eigen::Vector3d posC = C->getPosition();
		Eigen::Vector3d pos = edge->getHalfwayPosition();
		pos *= 1.25;
		boundaryEdgeVertexCoord.col(idx++) = pos;
		boundaryEdgeVertexIdx.coeffRef(edge->getIndex()) = i + offset_E;

		const Vertex * V = addVertex(pos);
		assert(V->getIndex() == i + offset_E);
	}
	assert(idx == nnzE);
	vertexCoordinates.conservativeResize(3, offset_E + nnzE);
	vertexCoordinates.rightCols(nnzE) = boundaryEdgeVertexCoord;

	// create outer vertices for each boundary vertex
	const int nnzV = boundaryVertices.nonZeros();
	Eigen::SparseVector<int> boundaryVertexVertexIdx(nVertices);
	Eigen::Matrix3Xd boundaryVertexVertexCoord(3, nnzV);
	const int offset_V = vertexCoordinates.cols();
	idx = 0;
	for (int i = 0; i < nnzV; i++) {
		const int idxV = boundaryVertices.innerIndexPtr()[i];
		boundaryVertexVertexIdx.coeffRef(idxV) = i + offset_V;
		const Vertex * vertex = getVertex(idxV);
		Eigen::Vector3d pos = vertex->getPosition();
		pos *= 1.5;
		boundaryVertexVertexCoord.col(idx++) = pos;

		const Vertex * V = addVertex(pos);
		assert(V->getIndex() == i + offset_V);
	}
	assert(idx == nnzV);
	vertexCoordinates.conservativeResize(3, offset_V + nnzV);
	vertexCoordinates.rightCols(nnzV) = boundaryVertexVertexCoord;

	// create faces at boundary edge
	Eigen::Matrix3Xi f2v_edges(3, 2*nnzE);
	f2v_edges.array() = -1;
	idx = 0;
	for (int i = 0; i < nnzE; i++) {
		const int idxE = boundaryEdges.innerIndexPtr()[i];
		Edge * edge = edgeList[idxE];
		const int idxC = boundaryEdgeVertexIdx.coeff(edge->getIndex());
		Vertex * C1 = getVertex(idxC);
		assert(C1 != nullptr);
		Vertex * A = edge->getVertexA();
		Vertex * B = edge->getVertexB();
		f2v_edges.col(idx++) = Eigen::Vector3i(A->getIndex(), B->getIndex(), C1->getIndex());
	
		const int idxA1 = boundaryVertexVertexIdx.coeff(A->getIndex());
		const int idxB1 = boundaryVertexVertexIdx.coeff(B->getIndex());
		Vertex * A1 = getVertex(idxA1);
		Vertex * B1 = getVertex(idxB1);
		f2v_edges.col(idx++) = Eigen::Vector3i(A1->getIndex(), B1->getIndex(), C1->getIndex());
	}
	assert(f2v_edges.cols() == idx);

	// create faces at boundary vertices
	Eigen::Matrix3Xi f2v_vertices(3, 2*nnzV);
	f2v_vertices.array() = -1;
	idx = 0;
	for (int i = 0; i < nnzV; i++) {
		const int idxV = boundaryVertices.innerIndexPtr()[i];
		Vertex * bV = getVertex(idxV);

		// get adjacient boundary edges
		std::vector<Edge*> vertexEdges = bV->getEdges();
		std::vector<Edge*> boundaryEdgeList;
		for (auto edge : vertexEdges) {
			if ((edge->getIndex() < nEdges) && edge->isBoundary()) {
				boundaryEdgeList.push_back(edge);
			}
		}
		assert(boundaryEdgeList.size() == 2);

		Vertex * C1 = getVertex(boundaryEdgeVertexIdx.coeff(boundaryEdgeList[0]->getIndex()));
		Vertex * D1 = getVertex(boundaryEdgeVertexIdx.coeff(boundaryEdgeList[1]->getIndex()));
		Vertex * B1 = getVertex(boundaryVertexVertexIdx.coeff(bV->getIndex()));

		f2v_vertices.col(idx++) = Eigen::Vector3i(bV->getIndex(), C1->getIndex(), D1->getIndex());
		f2v_vertices.col(idx++) = Eigen::Vector3i(B1->getIndex(), C1->getIndex(), D1->getIndex());
	}
	assert(f2v_vertices.cols() == idx);

	// add faces
	for (int i = 0; i < f2v_edges.cols(); i++) {
		Vertex * A = getVertex(f2v_edges(0, i));
		Vertex * B = getVertex(f2v_edges(1, i));
		Vertex * C = getVertex(f2v_edges(2, i));
		addFace({ addEdge(A, B), addEdge(B, C), addEdge(C, A) });
	}
	for (int i = 0; i < f2v_vertices.cols(); i++) {
		Vertex * A = getVertex(f2v_vertices(0, i));
		Vertex * B = getVertex(f2v_vertices(1, i));
		Vertex * C = getVertex(f2v_vertices(2, i));
		addFace({ addEdge(A, B), addEdge(B, C), addEdge(C, A) });
	}
}

void PrimalMesh::extrudeMesh(const Eigen::Vector3d & extrudeDirection, const int nLayers)
{
	assert(extrudeDirection.norm() > 0);
	assert(extrudeDirection(2) > 0);

	const int nVertices_2d = vertexList.size();
	const int nEdges_2d = edgeList.size();
	const int nFaces_2d = faceList.size();

	// Create vertices
	for (int layer = 1; layer <= nLayers; layer++) {
		for (int i = 0; i < nVertices_2d; i++) {
			Eigen::Vector3d pos = vertexList[i]->getPosition() + layer * ((1./nLayers) * extrudeDirection);
			addVertex(pos);
		}
	}
	
	// Create faces 
	for (int layer = 1; layer <= nLayers; layer++) {
		// Create faces parallel to initial mesh
		for (int j = 0; j < nFaces_2d; j++) {
			const Face * face = faceList[j];
			std::vector<Edge*> faceEdges = face->getEdgeList();
			std::vector<Edge*> newFaceEdges;
			for (int k = 0; k < faceEdges.size(); k++) {
				Edge * refEdge = faceEdges[k];

				Vertex * A = refEdge->getVertexA();
				Vertex * B = refEdge->getVertexB();
				Vertex * Anew = getVertex(A->getIndex() + layer * nVertices_2d);
				Vertex * Bnew = getVertex(B->getIndex() + layer * nVertices_2d);
				Edge * e = addEdge(Anew, Bnew);
				newFaceEdges.push_back(e);
			}
			addFace(newFaceEdges);
		}
	}
	for (int layer = 1; layer <= nLayers; layer++) {
		// Create faces parallel to extrude direction
		for (int i = 0; i < nEdges_2d; i++) {
			Edge * refEdge = edgeList[i];
			Vertex * Aref = refEdge->getVertexA();
			Vertex * Bref = refEdge->getVertexB();

			Vertex * A0 = getVertex(Aref->getIndex() + (layer - 1)* nVertices_2d);
			Vertex * A1 = getVertex(Aref->getIndex() + layer * nVertices_2d);
			Vertex * B0 = getVertex(Bref->getIndex() + (layer - 1) * nVertices_2d);
			Vertex * B1 = getVertex(Bref->getIndex() + layer * nVertices_2d);
			addFace({ addEdge(A0, B0), addEdge(B0, B1), addEdge(B1, A1), addEdge(A1, A0) });
		}
	}

}

Eigen::Matrix3Xi PrimalMesh::refine_triangles()
{
	const int nVertices = vertexList.size();
	assert(nVertices > 0);

	// copy vertices
	vertexCoordinates = Eigen::Matrix3Xd(3, nVertices);
	for (int i = 0; i < nVertices; i++) {
		vertexCoordinates.col(i) = vertexList[i]->getPosition();
	}

	// define edge midpoint vertices
	const int nEdges = edgeList.size();
	Eigen::Matrix3Xd edgeMidpoints(3, nEdges);
	for (int i = 0; i < nEdges; i++) {
		const Edge * edge = edgeList[i];
		Eigen::Vector3d pos = edge->getHalfwayPosition();
		if (edge->isBoundary()) {
			pos /= pos.norm();
		}
		edgeMidpoints.col(i) = pos;
	}
	vertexCoordinates.conservativeResize(3, nVertices + nEdges);
	vertexCoordinates.rightCols(nEdges) = edgeMidpoints;

	// Refine faces into four subfaces, defined by edge midpoints
	const int nFaces = faceList.size();
	Eigen::Matrix3Xi f2v(3, 4*nFaces);
	f2v.array() = -1;

	int idx = 0;
	for (int i = 0; i < nFaces; i++) {
		Face * face = faceList[i];
		std::vector<Edge*> faceEdges = face->getEdgeList();
		assert(faceEdges.size() == 3);
		Edge * edge0 = faceEdges[0];
		Edge * edge1 = faceEdges[1];
		Edge * edge2 = faceEdges[2];

		assert(edge1->hasCoincidentVertex(edge2));
		Vertex * A = edge1->getCoincidentVertex(edge2);
		assert(edge2->hasCoincidentVertex(edge0));
		Vertex * B = edge2->getCoincidentVertex(edge0);
		assert(edge0->hasCoincidentVertex(edge1));
		Vertex * C = edge0->getCoincidentVertex(edge1);

		int v0 = A->getIndex();
		int v1 = B->getIndex();
		int v2 = C->getIndex();
		int e0 = faceEdges[0]->getIndex() + nVertices;
		int e1 = faceEdges[1]->getIndex() + nVertices;
		int e2 = faceEdges[2]->getIndex() + nVertices;

		f2v.col(idx++) = Eigen::Vector3i(v0, e1, e2);
		f2v.col(idx++) = Eigen::Vector3i(e0, v1, e2);
		f2v.col(idx++) = Eigen::Vector3i(e0, e1, v2);
		f2v.col(idx++) = Eigen::Vector3i(e0, e1, e2);
	}
	assert(idx == f2v.cols());
	assert((f2v.array() >= 0).all());
	return f2v;
}

Eigen::Matrix3Xi PrimalMesh::refine_triangles_specialCorners()
{
	const int nVertices = vertexList.size();
	assert(nVertices > 0);
	const int nEdges = edgeList.size();
	const int nFaces = faceList.size();
	Eigen::Matrix3Xi f2v(3, 0);


	// copy vertices
	vertexCoordinates = Eigen::Matrix3Xd(3, nVertices);
	for (int i = 0; i < nVertices; i++) {
		vertexCoordinates.col(i) = vertexList[i]->getPosition();
	}

	// Indicator for special edges (i.e., those that are split in half instead of thirds)
	Eigen::VectorXi specialEdges = Eigen::VectorXi::Zero(nEdges);
	for (int i = 0; i < nEdges; i++) {
		Edge * edge = edgeList[i];
		const int nA = edge->getVertexA()->getEdges().size();
		const int nB = edge->getVertexB()->getEdges().size();
		specialEdges(i) = (nA == 3 || nB == 3);
	}

	// Create vertices at one-third and two-third of each edge;
	// except for special edges, where we create vertex at midpoint.
	// We keep track of the mapping from edge to new vertices.
	Eigen::Matrix3Xd edgeInnerCoords(3, 2 * nEdges - specialEdges.array().sum());
	assert(edgeInnerCoords.cols() > nEdges);

	Eigen::Matrix2Xi edge2innerVertices(2, nEdges);
	edge2innerVertices.array() = -1;
	int idx = 0;
	for (int i = 0; i < nEdges; i++) {
		Edge * edge = edgeList[i];
		const Eigen::Vector3d A = edge->getVertexA()->getPosition();
		const Eigen::Vector3d B = edge->getVertexB()->getPosition();
		const Eigen::Vector3d edgeVector = edge->getDirection();

		if (specialEdges(i)) {
			Eigen::Vector3d pos = A + 0.5 * edgeVector;
			if (edge->isBoundary()) {
				pos /= pos.norm();
			}
			edgeInnerCoords.col(idx) = pos;
			edge2innerVertices(0, i) = nVertices + idx++;
		}
		else {
			const Eigen::Vector3d pos0 = A + 1. / 3. * edgeVector;
			edgeInnerCoords.col(idx) = pos0;
			edge2innerVertices(0, i) = nVertices + idx++;
			const Eigen::Vector3d pos1 = A + 2. / 3. * edgeVector;
			edgeInnerCoords.col(idx) = pos1;
			edge2innerVertices(1, i) = nVertices + idx++;
		}
	}
	vertexCoordinates.conservativeResize(3, nVertices + edgeInnerCoords.cols());
	vertexCoordinates.rightCols(edgeInnerCoords.cols()) = edgeInnerCoords;

	// new vertices at face center 
	int faceOffset = vertexCoordinates.cols();
	Eigen::Matrix3Xd fc(3, nFaces);
	for (int i = 0; i < nFaces; i++) {
		Eigen::Vector3d center;
		center.setZero();
		const Face * face = faceList[i];

		// face center by arithmetic mean of vertices
		std::vector<Vertex*> faceVertices = face->getVertexList();
		assert(faceVertices.size() == 3);
		for (int j = 0; j < 3; j++) {
			center += faceVertices[j]->getPosition();
		}
		center /= 3;
		fc.col(i) = center;
	}
	vertexCoordinates.conservativeResize(3, faceOffset + nFaces);
	vertexCoordinates.rightCols(nFaces) = fc;

	// For each face ...
	f2v = Eigen::Matrix3Xi(3, 9 * nFaces);
	f2v.array() = -1;
	int fidx = 0;
	for (int i = 0; i < nFaces; i++) {
		const Face * face = faceList[i];
		const std::vector<Edge*> faceEdges = face->getEdgeList();
		if (face->isBoundary()) {
			// ... boundary face ...
			// find local index of not-special edge
			int idx = 0;
			for (idx = 0; idx < 3; idx++) {
				if (specialEdges(faceEdges[idx]->getIndex()) == 0) {
					break;
				}
			}
			// edge0 is not a special edge, but edge1 and edge2 are a special edge
			const Edge * edge0 = faceEdges[idx];
			const Edge * edge1 = faceEdges[(idx + 1) % 3];
			const Edge * edge2 = faceEdges[(idx + 2) % 3];
			assert(specialEdges(edge0->getIndex()) == 0);
			assert(specialEdges(edge1->getIndex()) == 1);
			assert(specialEdges(edge2->getIndex()) == 1);

			int v0 = edge1->getCoincidentVertex(edge2)->getIndex();
			int v1 = edge2->getCoincidentVertex(edge0)->getIndex();
			int v2 = edge0->getCoincidentVertex(edge1)->getIndex();

			// points at 1/3 and 2/3 of edge0 
			const int incidence = edge0->getIncidence(getVertex(v1));
			const int e0 = edge0->getIndex();
			int v01 = (incidence == 1) ? edge2innerVertices(0, e0) : edge2innerVertices(1, e0);
			int v02 = (incidence == 1) ? edge2innerVertices(1, e0) : edge2innerVertices(0, e0);
			assert(v01 >= 0);
			assert(v02 >= 0);

			// Midpoint vertices on special edges
			int v11 = edge2innerVertices(0, edge1->getIndex());
			int v22 = edge2innerVertices(0, edge2->getIndex());

			int vc = face->getIndex() + faceOffset;

			// Define refined faces
			f2v.col(fidx++) = Eigen::Vector3i(vc, v01, v02);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v02, v11);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v11, v0);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v0, v22);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v22, v01);
			f2v.col(fidx++) = Eigen::Vector3i(v1, v01, v22);
			f2v.col(fidx++) = Eigen::Vector3i(v2, v11, v02);
		}
		else {
			// ... not boundary face ...
			const Edge * edge0 = faceEdges[0];
			const Edge * edge1 = faceEdges[1];
			const Edge * edge2 = faceEdges[2];

			const int v0 = edge1->getCoincidentVertex(edge2)->getIndex();
			const int v1 = edge2->getCoincidentVertex(edge0)->getIndex();
			const int v2 = edge0->getCoincidentVertex(edge1)->getIndex();
			
			const int incidence0 = edge0->getIncidence(getVertex(v1));
			const int incidence1 = edge1->getIncidence(getVertex(v0));
			const int incidence2 = edge2->getIncidence(getVertex(v0));

			const int e0 = edge0->getIndex();
			const int e1 = edge1->getIndex();
			const int e2 = edge2->getIndex();
			int v01 = (incidence0 == 1) ? edge2innerVertices(0, e0) : edge2innerVertices(1, e0);
			int v02 = (incidence0 == 1) ? edge2innerVertices(1, e0) : edge2innerVertices(0, e0);
			int v10 = (incidence1 == 1) ? edge2innerVertices(0, e1) : edge2innerVertices(1, e1);
			int v12 = (incidence1 == 1) ? edge2innerVertices(1, e1) : edge2innerVertices(0, e1);
			int v20 = (incidence2 == 1) ? edge2innerVertices(0, e2) : edge2innerVertices(1, e2);
			int v21 = (incidence2 == 1) ? edge2innerVertices(1, e2) : edge2innerVertices(0, e2);

			const int vc = face->getIndex() + faceOffset;
			// Define refined faces
			f2v.col(fidx++) = Eigen::Vector3i(vc, v10, v20);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v20, v21);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v21, v01);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v01, v02);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v02, v12);
			f2v.col(fidx++) = Eigen::Vector3i(vc, v12, v10);
			f2v.col(fidx++) = Eigen::Vector3i(v0, v10, v20);
			f2v.col(fidx++) = Eigen::Vector3i(v1, v01, v21);
			f2v.col(fidx++) = Eigen::Vector3i(v2, v02, v12);
		}
	}
	f2v.conservativeResize(3, fidx);
	assert((f2v.array() >= 0).all());
	return f2v; 
}

void PrimalMesh::test_quadFace()
{
	Vertex * A = addVertex(Eigen::Vector3d(0, 0, 0));
	Vertex * B = addVertex(Eigen::Vector3d(1, 0, 0));
	Vertex * C = addVertex(Eigen::Vector3d(1, 1, 0));
	Vertex * D = addVertex(Eigen::Vector3d(0, 1, 0));
	Edge * edge0 = addEdge(A, B);
	Edge * edge1 = addEdge(B, C);
	Edge * edge2 = addEdge(C, D);
	Edge * edge3 = addEdge(D, A);
	addFace({ edge0, edge1, edge2, edge3 });
}
