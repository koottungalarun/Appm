#include "InelasticCollision.h"

InelasticCollision::InelasticCollision()
{
}

InelasticCollision::InelasticCollision(const std::string & filename) 
	: InelasticCollision(filename, 1.0, 1.0)
{
}

InelasticCollision::InelasticCollision(const std::string & filename, const double kScale, const double Tscale)
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
		iss >> T >> c1 >> ki;

		// Store data from file in data vectors
		Tvec.push_back(T / Tscale);
		kiVec.push_back(ki / kScale);
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
Eigen::VectorXd InelasticCollision::getIonizationRate(const Eigen::VectorXd & T) const
{
	const int N = T.size();
	Eigen::VectorXd result = table->interpolate(T);
	assert(result.size() == N);
	return result;
}

Eigen::VectorXd InelasticCollision::getRecombinationRate_Saha(const Eigen::VectorXd & ki, const Eigen::VectorXd & Te) const
{
	const double a = getRecombSahaCoeff();

	Eigen::VectorXd kr(ki.size());
	kr = a * Te.array().pow(-3. / 2.) * ki.array();
	return kr;
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

int InelasticCollision::getElectronFluidx() const
{
	return fluidxElectrons;
}

void InelasticCollision::setElectronFluidx(const int idx)
{
	assert(idx >= 0);
	this->fluidxElectrons = idx;
}

int InelasticCollision::getAtomFluidx() const
{
	return fluidxAtoms;
}

void InelasticCollision::setAtomFluidx(const int idx)
{
	assert(idx >= 0);
	this->fluidxAtoms = idx;
}

int InelasticCollision::getIonFluidx() const
{
	return fluidxIons;
}

void InelasticCollision::setIonFluidx(const int idx)
{
	assert(idx >= 0);
	this->fluidxIons = idx;
}

void InelasticCollision::setScalingParameters(const ScalingParameters & params, const double electronMassRatio)
{
	const double mbar = params.getMassScale();
	const double nbar = params.getNumberDensityScale();
	const double Tbar = params.getTemperatureScale();

	PhysicsConstants & pc = PhysicsConstants::instance();

	const double h2 = pow(pc.planckConstant(), 2);
	const double kB = pc.kB();
	const double twoPi_inv = 0.5 * M_1_PI;
	const double temp = twoPi_inv * h2 / (mbar * kB * Tbar);
	const double deBroglieScaled_pow3 = pow(temp, 3./2.);
	const double g0 = 1;
	const double g1 = 6;
	const double a = 2*g1/g0 * 1./(deBroglieScaled_pow3 * nbar) * pow(electronMassRatio,-3./2.);
	setRecombSahaCoeff(a);
}

const double InelasticCollision::getRecombSahaCoeff() const
{
	assert(this->recomb_Saha_coeff > 0);
	return this->recomb_Saha_coeff;
}

void InelasticCollision::setRecombSahaCoeff(const double a)
{
	assert(a > 0);
	assert(isfinite(a));
	this->recomb_Saha_coeff = a;
}
