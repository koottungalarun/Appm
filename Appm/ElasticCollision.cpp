#include "ElasticCollision.h"

ElasticCollision::ElasticCollision()
{
}

ElasticCollision::ElasticCollision(const std::string & filename, const int idxA, const int idxB)
	: ElasticCollision(filename, idxA, idxB, 1., 1.)
{
}

ElasticCollision::ElasticCollision(const std::string & filename, const int idxA, const int idxB, const double Tscale, const double crossSecScale)
{
	assert(filename.size() > 3);
	assert(idxA >= 0);
	assert(idxB >= 0);
	assert(Tscale > 0);
	assert(crossSecScale > 0);

	this->fluidxA = idxA;
	this->fluidxB = idxB;

	std::cout << "Collisions: " << filename << std::endl;
	std::cout << "fluid idx A, B: " << idxA << ", " << idxB << std::endl;

	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cout << "File is not open: " << filename << std::endl;
		assert(false);
		exit(-1);
	}

	std::string line;
	std::vector<double> Tdata, Qdata;
	while (std::getline(file, line)) {
		if (line.substr(0, 1) == "#") { continue; } // skip comment lines
		double T, Q;
		char c;
		std::istringstream iss(line);
		iss >> T >> c >> Q;
		Tdata.push_back(T);
		Qdata.push_back(Q);
	}

	// Scaling of input data to dimensionless quantities
	for (int i = 0; i < Tdata.size(); i++) {
		Tdata[i] /= Tscale;
		Qdata[i] /= crossSecScale;
	}

	// create interpolation table
	table = new InterpolationTable(Tdata, Qdata);
}

ElasticCollision::~ElasticCollision()
{
	if (table != nullptr) {
		delete table;
		table = nullptr;
	}
}

const int ElasticCollision::getFluidIdxA() const
{
	assert(fluidxA >= 0);
	return this->fluidxA;
}

void ElasticCollision::setFluidIdxA(const int idxA)
{
	assert(idxA >= 0);
	this->fluidxA = idxA;
}

const int ElasticCollision::getFluidIdxB() const
{
	assert(fluidxB >= 0);
	return this->fluidxB;
}

void ElasticCollision::setFluidIdxB(const int idxB)
{
	assert(idxB >= 0);
	this->fluidxB = idxB;
}

const Eigen::VectorXd ElasticCollision::getAvgMomCrossSection(const Eigen::VectorXd & T) const 
{
	const int N = T.size();
	Eigen::VectorXd result = table->interpolate(T);
	assert(result.size() == N);
	//assert((result.array() > 0).all());
	return result;
}

const Eigen::MatrixXd ElasticCollision::getData() const
{
	Eigen::VectorXd x = table->getXdata();
	Eigen::VectorXd y = table->getYdata();
	Eigen::MatrixXd data(x.size(), 2);
	data.col(0) = x;
	data.col(1) = y;
	return data;
}
