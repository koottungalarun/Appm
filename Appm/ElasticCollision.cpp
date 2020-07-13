#include "ElasticCollision.h"

ElasticCollision::ElasticCollision()
{
}

ElasticCollision::ElasticCollision(const std::string & filename, const int idxA, const int idxB)
{
	assert(filename.size() > 3);
	assert(idxA >= 0);
	assert(idxB >= 0);
	//std::cout << "Read elastic collisions data: " << filename << "; idxA = " << idxA << "; idxB = " << idxB << std::endl;
	readFile(filename);
	this->fluidxA = idxA;
	this->fluidxB = idxB;
}

ElasticCollision::~ElasticCollision()
{
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

const Eigen::VectorXd ElasticCollision::getAvgMomCrossSection(const Eigen::VectorXd & T)
{
	const int N = T.size();
	Eigen::VectorXd result(N);
	result.setZero();
	return result;
}

void ElasticCollision::readFile(const std::string & filename)
{
	std::ifstream file(filename);
	if (!file.is_open()) { 
		std::cout << "File is not open: " << filename << std::endl;
		assert(false);
		exit(-1);
	}

	std::string line;
	std::vector<double> xdata, ydata;
	while (std::getline(file, line)) {
		if (line.substr(0, 1) == "#") { continue; } // skip comment lines
		double x, y;
		char c;
		std::istringstream iss(line);
		iss >> x >> c >> y;
		xdata.push_back(x);
		ydata.push_back(y);
	}

	//for (int i = 0; i < xdata.size(); i++) {
	//	std::cout << xdata[i] << "\t" << ydata[i] << std::endl;
	//}

	// create interpolation table
	table = InterpolationTable(xdata, ydata);
}
