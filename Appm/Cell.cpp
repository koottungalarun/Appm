#include "Cell.h"



Cell::Cell()
{
}

Cell::Cell(const std::vector<Face*>& faces)
{
	//std::cout << "Create cell from faces: {";
	//for (auto face : faces) {
	//	std::cout << face->getIndex() << ",";
	//}
	//std::cout << "}" << std::endl;
	//std::cout << "Face vertices: " << std::endl;
	//for (auto face : faces) {
	//	std::cout << "face idx " << face->getIndex() << ": ";
	//	const std::vector<Vertex*> f2v = face->getVertexList();
	//	for (auto v : f2v) {
	//		std::cout << v->getIndex() << ",";
	//	}
	//	std::cout << std::endl;
	//}

	this->faceList = faces;
	for (auto face : faceList) {
		face->setAdjacient(this);
	}

	// determine cell center
	this->center = computeCellCenter();
	assert(validateCellGeometry());

	// Set cell volume
	volume = 0;
	for (auto face : faceList) {
		const double fA = face->getArea();
		const Eigen::Vector3d fn = face->getNormal();
		double h = std::abs(fn.dot(center - face->getCenter()));
		double dV = 1. / 3. * fA * h;
		assert(dV > 0);
		volume += dV;
	}
	assert(volume > 0);

	//Eigen::Matrix3d M;
	//M.col(0) = zFaces[0]->getCenter();
	//M.col(1) = zFaces[1]->getCenter();
	//M.col(2) = this->center;
	//std::cout << "face/cell center: " << std::endl << M << std::endl;
}


Cell::~Cell()
{
}

const double Cell::getVolume() const
{
	return volume;
}

const std::vector<Face*> Cell::getFaceList() const
{
	return faceList;
}

bool Cell::hasFace(const Face * face) const
{
	assert(face != nullptr);
	for (auto f : faceList) {
		if (f == face) {
			return true;
		}
	}
	return false;
}

const int Cell::getOrientation(const Face * face) const
{
	assert(face != nullptr);
	bool isMember = false;
	for (auto f : faceList) {
		if (f == face) {
			isMember |= true;
			break;
		}
	}
	assert(isMember);
	const Eigen::Vector3d a = face->getCenter() - this->center;
	const Eigen::Vector3d fn = face->getNormal();
	const int orientation = (a.dot(fn) > 0) ? 1 : -1;
	return orientation;
}

const Eigen::Vector3d & Cell::getCenter() const
{
	return center;
}

const Eigen::Matrix3Xd Cell::getVertexCoordinates() const
{
	std::vector<Vertex*> vertexList;
	for (auto face : faceList) {
		for (auto vertex : face->getVertexList()) {
			vertexList.push_back(vertex);
		}
	}
	std::sort(vertexList.begin(), vertexList.end());
	std::vector<Vertex*>::iterator it;
	it = std::unique(vertexList.begin(), vertexList.end());
	int nVertices = std::distance(vertexList.begin(), it);
	assert(nVertices >= 4);

	Eigen::Matrix3Xd coords(3, nVertices);
	for (int i = 0; i < nVertices; i++) {
		const Vertex * vertex = vertexList[i];
		coords.col(i) = vertex->getPosition();
	}
	return coords;
}

void Cell::setType(const Type & type)
{
	this->fluidType = type;
}

const Cell::Type Cell::getType() const
{
	return this->fluidType;
}

bool Cell::validateCellGeometry() const
{
	bool isValid = true;
	auto cc = this->getCenter();
	auto cellFaces = this->getFaceList();
	for (auto face : cellFaces) {
		auto fc = face->getCenter();
		auto posVec = fc - cc;
		const double x = posVec(0);
		const double y = posVec(1);
		const double z = posVec(2);
		bool isParallel = (x == 0) && (y == 0) && (z != 0);
		bool isPerpendicular = ((x != 0) || (y != 0)) && (z == 0);
		isValid &= isParallel ^ isPerpendicular;// boolean xor operator (^)
	}
	if (!isValid) {
		std::cout << "Validate cel geometry" << std::endl;
		std::cout << "Cell idx: " << this->getIndex() << std::endl;
		std::cout << "cell center: " << cc.transpose() << std::endl;

		for (auto face : cellFaces) {
			std::cout << "face idx, center: ";
			std::cout << face->getIndex() << ",\t";
			std::cout << face->getCenter().transpose() << std::endl;
		}

		Eigen::VectorXd zValues(cellFaces.size());
		std::cout << "Position vector:" << std::endl;
		
		for (int idx = 0; idx < cellFaces.size(); idx++) {
			auto face = cellFaces[idx];
			auto fc = face->getCenter();
			auto posVec = fc - cc;
			zValues(idx) = posVec(2);
			std::cout << std::scientific << posVec.transpose() << std::endl;
		}
		std::cout << "zValues:" << std::endl;
		for (int i = 0; i < zValues.size(); i++) {
			std::cout << std::scientific << zValues(i) << std::endl;
		}
	}
	return isValid;
}

Eigen::Vector3d Cell::computeCellCenter() const
{
	//std::vector<Face*> zFaces;
	//for (auto face : faceList) {
	//	const Eigen::Vector3d fn = face->getNormal();
	//	if (fn.cross(Eigen::Vector3d(0, 0, 1)).norm() < 100 * std::numeric_limits<double>::epsilon()) {
	//		zFaces.push_back(face);
	//	}
	//}
	//const int n_zFaces = zFaces.size();
	//if (n_zFaces < 2) {
	//	for (auto face : faceList) {
	//		const Eigen::Vector3d fn = face->getNormal();
	//		const double fn_x_zUnit_norm = fn.cross(Eigen::Vector3d(0, 0, 1)).norm();
	//		std::cout << "Face " << face->getIndex() << ": " << "(fn x zUnit).norm() = " << fn_x_zUnit_norm << std::endl;
	//		if (fn.cross(Eigen::Vector3d(0, 0, 1)).norm() < 100 * std::numeric_limits<double>::epsilon()) {
	//			zFaces.push_back(face);
	//		}
	//	}
	//}
	//assert(zFaces.size() == 2);
	//Eigen::Vector3d cc = 0.5 * (zFaces[0]->getCenter() + zFaces[1]->getCenter());

	// Show list of face centers
	bool showListOfFaceCenters = false;
	if (showListOfFaceCenters) {
		std::cout << "List of face centers:" << std::endl;
		for (auto face : faceList) {
			std::cout << face->getCenter().transpose() << std::endl;
		}
		std::cout << std::endl;
	}

	bool isValid = true;
	Eigen::Vector3d cc;
	cc.setZero();

	// Create lists of faces that are parallel and perpendicular to z-axis
	std::vector<Face*> parallelFaces, perpendicularFaces;
	for (auto face : faceList) {
		const Eigen::Vector3d fn = face->getNormal();
		if (fn(2) != 0) {
			const double tol = 2 * std::numeric_limits<double>::epsilon();
			assert(abs(abs(fn(2)) - 1) < tol); // z-component is +/- 1.0
			perpendicularFaces.push_back(face);
		}
		else {
			parallelFaces.push_back(face);
		}
	}
	assert(perpendicularFaces.size() == 2);

	// Check that perpendicular faces have same x-y coordinates in face centers
	Eigen::Vector3d posVec1, posVec2;
	posVec1 = perpendicularFaces[0]->getCenter();
	posVec2 = perpendicularFaces[1]->getCenter();
	Eigen::Vector3d posVecDelta = posVec2 - posVec1;
	Eigen::Vector2d posVec2d = posVecDelta.segment(0, 2);
	if (posVec2d.norm() != 0) { // this should not be true
		std::cout << "posVec2d: " << posVec2d.transpose() << std::endl;
		isValid &= false;
	}
	cc = 0.5 * (posVec1 + posVec2);

	// Check that parallel faces have same z-value in face centers
	Eigen::VectorXd zValuesParallel(parallelFaces.size());
	for (int i = 0; i < parallelFaces.size(); i++) {
		auto face = parallelFaces[i];
		zValuesParallel(i) = face->getCenter()(2);
	}
	isValid &= zValuesParallel.array().isApproxToConstant(zValuesParallel(0));
	isValid &= (zValuesParallel.array() == zValuesParallel(0)).all();
	//cc = posVec1;
	cc(2) = zValuesParallel(0);

	// Show data
	if (!isValid) {
		std::cout << "is not valid" << std::endl;
		std::cout << "face center coords:" << std::endl;
		for (auto face : faceList) {
			std::cout << std::scientific << face->getCenter().transpose() << std::endl;
		}
	}
	assert(isValid);

	// Check if local position vector from cell center to face center
	// is either parallel or perpendicular to z-axis
	isValid = true;
	for (auto face : faceList) {
		const Eigen::Vector3d fc = face->getCenter();
		const Eigen::Vector3d posVec = fc - cc;

		if (posVec(2) != 0) {
			isValid = posVec.segment(0, 2).norm() == 0;
		}
		else {
			isValid = posVec.segment(0, 2).norm() != 0;
		}
	}

	if (!isValid) {
		std::cout << "Cell center:" << std::endl;
		std::cout << cc.transpose() << std::endl;
		std::cout << "List of face centers:" << std::endl;
		for (auto face : faceList) {
			std::cout << face->getCenter().transpose() << std::endl;
		}
	}
	assert(isValid);



	//std::cout << "z-Faces: ";
	//for (auto f : zFaces) {
	//	std::cout << f->getIndex() << ",";
	//}
	//std::cout << std::endl;
	
	//this->center = Eigen::Vector3d::Zero();

	// Check position of cell-center relative to face centers





	return cc;
}

