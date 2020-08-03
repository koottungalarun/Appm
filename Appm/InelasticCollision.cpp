#include "InelasticCollision.h"

InelasticCollision::InelasticCollision()
{
}

InelasticCollision::InelasticCollision(const std::string & filename)
{
	assert(filename.size() > 0);
	
	std::vector<double> Tvec, kiVec;

	// Open data file
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cout << "File is not open: " << filename << std::endl;
		assert(false);
	}
	std::string line;

	int idx = 0;
	// Read data file line by line
	while (std::getline(file, line)) {
		if (line.substr(0, 1) == "#") { continue; } // skip comment lines

		// Read one line
		double T, ki, kr;
		char c1, c2;
		std::istringstream iss(line);
		iss >> T >> c1 >> ki >> c2 >> kr;

		// Store data from file in data vectors
		Tvec.push_back(T);
		kiVec.push_back(ki);
	}
	table = new InterpolationTable(Tvec, kiVec);
}

InelasticCollision::~InelasticCollision()
{
	if (table != nullptr) {
		delete table;
		table = nullptr;
	}
}

/**
* Get ionization rate at temperature T.
*
* @param T electron temperature
* @return ionization rate k_i (T).
*/
Eigen::VectorXd InelasticCollision::getIonizationRate(const Eigen::VectorXd & T)
{
	const int N = T.size();
	Eigen::VectorXd result = table->interpolate(T);
	assert(result.size() == N);
	return result;
}

Eigen::MatrixXd InelasticCollision::getData() const
{
	Eigen::VectorXd x = table->getXdata();
	Eigen::VectorXd y = table->getYdata();
	Eigen::MatrixXd data(x.size(), 2);
	data.col(0) = x;
	data.col(1) = y;
	return data;
}
