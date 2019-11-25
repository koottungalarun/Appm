#include "DualMesh.h"



DualMesh::DualMesh()
	: Mesh("dual")
{
}

DualMesh::DualMesh(const std::string & meshPrefix) 
	: Mesh(meshPrefix)
{
}

DualMesh::~DualMesh()
{
}

void DualMesh::init_dualMesh(const PrimalMesh & primal)
{
	// add dual vertices ... 
	
	// ... at primal cell centers
	const std::vector<Cell*> primalCells = primal.getCells();
	for (auto cell : primalCells) {
		addVertex(cell->getCenter());
	}
	// ... at primal boundary face centers
	const int nPrimalFaces = primal.getNumberOfFaces();
	Eigen::SparseVector<int> primalFaceToDualVertex(nPrimalFaces);
	for (int i = 0; i < nPrimalFaces; i++) {
		const Face * face = primal.getFace(i);
		if (face->isBoundary()) {
			const Vertex * V = addVertex(face->getCenter());
			primalFaceToDualVertex.coeffRef(i) = V->getIndex();
		}
	}
	// ... at primal boundary edge centers
	const int nPrimalEdges = primal.getNumberOfEdges();
	Eigen::SparseVector<int> primalEdgeToDualVertex(nPrimalEdges);
	for (int i = 0; i < nPrimalEdges; i++) {
		const Edge * edge = primal.getEdge(i);
		if (edge->isBoundary()) {
			// get adjacient boundary faces; 
			// if their face normals are not parallel, add a vertex
			std::vector<Face*> boundaryFaces;
			for (auto face : edge->getFaceList()) {
				if (face->isBoundary()) {
					boundaryFaces.push_back(face);
				}
			}
			assert(boundaryFaces.size() == 2);
			Eigen::Vector3d fn0 = boundaryFaces[0]->getNormal();
			Eigen::Vector3d fn1 = boundaryFaces[1]->getNormal();
			if (fn0.cross(fn1).norm() > 4*std::numeric_limits<double>::epsilon()) {
				const Vertex * V = addVertex(edge->getHalfwayPosition());
				primalEdgeToDualVertex.coeffRef(i) = V->getIndex();
			}
		}
	}
	// ... at primal boundary vertices
	const int nPrimalVertices = primal.getNumberOfVertices();
	Eigen::SparseVector<int> primalVertexToDualVertex(nPrimalVertices);
	for (int i = 0; i < nPrimalVertices; i++) {
		const Vertex * vertex = primal.getVertex(i);
		if (vertex->isBoundary()) {
			const Vertex * V = addVertex(vertex->getPosition());
			primalVertexToDualVertex.coeffRef(i) = V->getIndex();
		}
	}

	// add dual edges across primal faces ...
	const std::vector<Face*> primalFaces = primal.getFaces();
	for (auto face : primalFaces) {
		const std::vector<Cell*> adjCells = face->getCellList();

		int idxA = -1;
		int idxB = -1;
		if (face->isBoundary()) {
			assert(adjCells.size() == 1);
			const Cell * cell = adjCells[0];
			const int orientation = cell->getOrientation(face);
			assert(abs(orientation) == 1);
			if (orientation > 0) {
				idxA = cell->getIndex();
				idxB = primalFaceToDualVertex.coeff(face->getIndex());
			}
			else {
				idxA = primalFaceToDualVertex.coeff(face->getIndex());
				idxB = cell->getIndex();
			}
		}
		else {
			assert(adjCells.size() == 2);
			for (auto cell : adjCells) {
				const int orientation = cell->getOrientation(face);
				if (orientation > 0) {
					idxA = cell->getIndex();
				}
				else {
					idxB = cell->getIndex();
				}				
			}
		}
		assert(idxA >= 0);
		assert(idxB >= 0);
		Vertex * A = getVertex(idxA);
		Vertex * B = getVertex(idxB);
		addEdge(A, B);
	}

	// Define dual faces around primal edges
	for (int i = 0; i < nPrimalEdges; i++) {
		const Edge * primalEdge = primal.getEdge(i);
		if (primalEdge->isBoundary()) {
			// get adjacient primal boundary faces connect their face centers
			const std::vector<Face*> primalEdgeFaces = primalEdge->getFaceList();
			std::vector<Face*> primalBoundaryFaces;
			for (auto face : primalEdgeFaces) {
				if (face->isBoundary()) {
					primalBoundaryFaces.push_back(face);
				}
			}
			assert(primalBoundaryFaces.size() == 2);
			const int idxA = primalFaceToDualVertex.coeff(primalBoundaryFaces[0]->getIndex());
			const int idxB = primalFaceToDualVertex.coeff(primalBoundaryFaces[1]->getIndex());
			assert(idxA >= 0);
			assert(idxB >= 0);
			Vertex * A = getVertex(idxA);
			Vertex * B = getVertex(idxB);
			std::vector<Edge*> dualEdges;

			// Check if the faces have collinear face normals; 
			// if yes, connect face centers.
			// if not, connect face centers with primal edge center
			const Eigen::Vector3d fn0 = primalBoundaryFaces[0]->getNormal();
			const Eigen::Vector3d fn1 = primalBoundaryFaces[1]->getNormal();
			if (fn0.cross(fn1).norm() < 4 * std::numeric_limits<double>::epsilon()) {
				Edge * dualEdge = addEdge(A, B);
				dualEdges.push_back(dualEdge);
			}
			else {
				const int idxC = primalEdgeToDualVertex.coeff(primalEdge->getIndex());
				assert(idxC > 0);
				Vertex * C = getVertex(idxC);
				Edge * dualEdgeAC = addEdge(A, C);
				Edge * dualEdgeBC = addEdge(C, B);
				dualEdges.push_back(dualEdgeAC);
				dualEdges.push_back(dualEdgeBC);
			}

			// add dual edges given by primal face normals
			for (auto primalFace : primalEdgeFaces) {
				const int pFidx = primalFace->getIndex();
				Edge * dualEdge = getEdge(pFidx);
				dualEdges.push_back(dualEdge);
			}

			//std::cout << "dual edges: " << std::endl;
			//for (auto edge : dualEdges) {
			//	std::cout << *edge << std::endl;
			//}
			// sort list of edges such that it forms a continuous loop
			std::vector<Edge*> dualEdgesLoop = makeContinuousLoop(dualEdges);
			addFace(dualEdgesLoop);
		}
		else {
			const std::vector<Face*> primalEdgeFaces = primalEdge->getFaceList();
			std::vector<Edge*> dualEdges;
			for (auto primalFace : primalEdgeFaces) {
				const int pFidx = primalFace->getIndex();
				dualEdges.push_back(getEdge(pFidx));
			}
			//std::cout << "dual face edges: {";
			//for (auto edge : dualEdges) {
			//	std::cout << edge->getIndex() << ",";
			//}
			//std::cout << "}" << std::endl;
			std::vector<Edge*> dualEdgesLoop = makeContinuousLoop(dualEdges);
			addFace(dualEdgesLoop);
		}
	}

	// add dual cells
	for (int i = 0; i < nPrimalVertices; i++) {
		const Vertex * primalVertex = primal.getVertex(i);
		std::cout << "Create dual cell at primal vertex (idx = " << primalVertex->getIndex() << ")" << std::endl;
		Vertex * dualVertex = getVertex(primalVertexToDualVertex.coeffRef(primalVertex->getIndex()));
		std::cout << "(= dual vertex " << dualVertex->getIndex() << ")" << std::endl;

		const std::vector<Edge*> primalVertexEdges = primalVertex->getEdges();
		bool hasPrimalBoundaryEdges = false;
		for (auto edge : primalVertexEdges) {
			hasPrimalBoundaryEdges |= edge->isBoundary();
		}
		
		if (hasPrimalBoundaryEdges) {
			std::vector<Face*> dualFaces;
			std::vector<Edge*> dualBoundaryEdges;

			// For each primal edge ...
			for (auto primalEdge : primalVertexEdges) {
				const int pEidx = primalEdge->getIndex();
				// ...  get dual faces ...
				Face * face = getFace(pEidx);
				dualFaces.push_back(face);

				// ... if primal boundary edge ...
				if (primalEdge->isBoundary()) {
					// ... get adjacient primal boundary faces ...
					std::vector<Face*> primalFaces = primalEdge->getFaceList();
					std::vector<Face*> primalBoundaryFaces;
					for (auto primalFace : primalFaces) {
						if (primalFace->isBoundary()) {
							primalBoundaryFaces.push_back(primalFace);
						}
					}
					assert(primalBoundaryFaces.size() == 0 || primalBoundaryFaces.size() == 2);
					// ... if it has primal boundary faces ... 
					if (primalBoundaryFaces.size() == 2) {
						// ... get dual edges ...
						// if primal buondary faces have parallel face normals, get dual edge between the face centers
						// otherwise, get two edges (each of them connecting edge center and a face center)
						const Eigen::Vector3d fn0 = primalBoundaryFaces[0]->getNormal();
						const Eigen::Vector3d fn1 = primalBoundaryFaces[1]->getNormal();
						const bool isFacesParallel = fn0.cross(fn1).norm() < 4 * std::numeric_limits<double>::epsilon();
						Vertex * A = getVertex(primalFaceToDualVertex.coeff(primalBoundaryFaces[0]->getIndex()));
						Vertex * B = getVertex(primalFaceToDualVertex.coeff(primalBoundaryFaces[1]->getIndex()));
						if (isFacesParallel) {
							Edge * dualBoundaryEdge = addEdge(A, B);
							dualBoundaryEdges.push_back(dualBoundaryEdge);
						}
						else {
							Vertex * C = getVertex(primalEdgeToDualVertex.coeff(primalEdge->getIndex()));
							Edge * dualBoundaryEdge0 = addEdge(A, C);
							Edge * dualBoundaryEdge1 = addEdge(C, B);
							dualBoundaryEdges.push_back(dualBoundaryEdge0);
							dualBoundaryEdges.push_back(dualBoundaryEdge1);
						}
					} // end if(primalBoundaryFaces.size() == 2)
				} // end if(primalEdge->isBoundary())
			} // end for (auto primalEdge : primalVertexEdges)

			//if (i >= 1) {
			//	break;
			//}

				
			// We are not sure if the dual boundary edges are all in a single plane.
			// (for instance: a corner in primal mesh)
			// Therefore, determine a reference normal vector and select all dual edges that yield to same normal vector.
			// This procedure may be repeated several times
			while (dualBoundaryEdges.size() > 0) {
				std::cout << "Dual boundary edges: " << std::endl;
				for (auto edge : dualBoundaryEdges) {
					std::cout << *edge << std::endl;
				}
				std::cout << std::endl;

				std::cout << "Add auxiliary dual face" << std::endl;
				// get reference normal direction
				const Eigen::Vector3d primalVertexPos = primalVertex->getPosition();
				Eigen::Vector3d posA = dualBoundaryEdges[0]->getVertexA()->getPosition();
				Eigen::Vector3d edgeDir = dualBoundaryEdges[0]->getDirection();
				Eigen::Vector3d a = (posA - primalVertexPos).normalized();
				Eigen::Vector3d b = edgeDir.normalized();
				const Eigen::Vector3d nVec_ref = a.cross(b).normalized();

				// select edges that make same normal vector
				std::vector<Edge*> selectedEdges;
				for (auto dualEdge : dualBoundaryEdges) {
					Eigen::Vector3d posA = dualEdge->getVertexA()->getPosition();
					Eigen::Vector3d edgeDir = dualEdge->getDirection();
					Eigen::Vector3d a = (posA - primalVertexPos).normalized();
					Eigen::Vector3d b = edgeDir.normalized();
					const Eigen::Vector3d nVec = a.cross(b).normalized();
					if (nVec.cross(nVec_ref).norm() < 4 * std::numeric_limits<double>::epsilon()) {
						selectedEdges.push_back(dualEdge);
					}
				}
				assert(selectedEdges.size() >= 2);

				std::cout << "Selected edges: " << std::endl;
				for (auto edge : selectedEdges) {
					std::cout << *edge << std::endl;
				}

				// remove selected edges from list
				const int nDualBoundaryEdges = dualBoundaryEdges.size();
				for (auto dualEdge : selectedEdges) {
					std::vector<Edge*>::iterator it;
					for (it = dualBoundaryEdges.begin(); it != dualBoundaryEdges.end(); it++) {
						if ((*it) == dualEdge) {
							dualBoundaryEdges.erase(it);
							break;
						}
					}
				}
				assert(dualBoundaryEdges.size() == nDualBoundaryEdges - selectedEdges.size());

				// Check if this loop is closed; if not, add auxiliary edges between dual vertex and open end vertices
				Eigen::SparseVector<int> visitorCounts(getNumberOfVertices());
				for (auto edge : selectedEdges) {
					const int idxA = edge->getVertexA()->getIndex();
					const int idxB = edge->getVertexB()->getIndex();
					visitorCounts.coeffRef(idxA) += 1;
					visitorCounts.coeffRef(idxB) += 1;
				}
				std::vector<Vertex*> openEndVertices;
				for (int k = 0; k < visitorCounts.nonZeros(); k++) {
					const int idx = visitorCounts.innerIndexPtr()[k];
					if (visitorCounts.coeff(idx) == 1) {
						openEndVertices.push_back(getVertex(idx));
					}
				}
				assert(openEndVertices.size() == 0 || openEndVertices.size() == 2);
				if (openEndVertices.size() == 2) {
					Edge * dualEdge0 = addEdge(openEndVertices[0], dualVertex);
					Edge * dualEdge1 = addEdge(openEndVertices[1], dualVertex);

					selectedEdges.push_back(dualEdge0);
					selectedEdges.push_back(dualEdge1);
					std::cout << "Open end vertices: " << openEndVertices[0]->getIndex() << ", " << openEndVertices[1]->getIndex() << std::endl;
					std::cout << "Add edges to close: " << std::endl;
					std::cout << *dualEdge0 << std::endl;
					std::cout << *dualEdge1 << std::endl;
					std::cout << "Selected edges (updated): " << std::endl;
					for (auto edge : selectedEdges) {
						std::cout << *edge << std::endl;
					}
				}

				std::vector<Edge*> selectedEdgesLoop = makeContinuousLoop(selectedEdges);
				Face * face = addFace(selectedEdgesLoop);
				dualFaces.push_back(face);
			} // end while (dualBoundaryEdges.size() > 0)
			std::cout << "Create dual cell from faces: " << "(n = " << dualFaces.size() << ") " << "{";
			for (auto face : dualFaces) {
				std::cout << face->getIndex() << ",";
			}
			std::cout << "}" << std::endl;
			addCell(dualFaces);
			//if (i >= 1) { break; }
		}
		else {
			std::vector<Face*> dualFaces;
			for (auto primalEdge : primalVertexEdges) {
				const int pEidx = primalEdge->getIndex();
				Face * face = getFace(pEidx);
				dualFaces.push_back(face);
			}
			std::cout << "Create dual cell from faces: {";
			for (auto face : dualFaces) {
				std::cout << face->getIndex() << ",";
			}
			std::cout << "}" << std::endl;
			addCell(dualFaces);
		}



	}

}
