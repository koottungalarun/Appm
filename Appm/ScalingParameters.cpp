#include "ScalingParameters.h"

const int ScalingParameters::undefinedScales()
{
	return (x0 <= 0) + (T0 <= 0) + (n0 <= 0) + (lambdaSq <= 0);
}

ScalingParameters::ScalingParameters()
{
	setLengthScale(1.0);
	setTemperatureScale(10e3);
	setScaledDebyeLength(1e-4);
	compute();
}

/**
* Read scaling parameters from input file. 
*/
ScalingParameters::ScalingParameters(const std::string & filename)
{
	assert(filename.length() > 0);

	// Open input file
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cout << "File could not be not opened: " << filename << std::endl;
	}
	else {
		std::cout << "Read file: " << filename << std::endl;
	}

	// Expected delimiter '=' at which the line is split 
	const char delim = '=';

	// Read contents
	std::string line;
	bool isValid = true;
	while (std::getline(file, line) && isValid) {
		//std::cout << line << std::endl;
		if (line.length() == 0) { continue; } 		// Skip empty line
		if (line.front() == '#') { continue; }		// Skip comment line 

		int pos = line.find(delim);
		std::string tag = line.substr(0, pos);

		if (tag == "x") {
			double x = 0;
			std::istringstream(line.substr(pos + 1)) >> x;
			setLengthScale(x);
		}
		if (tag == "T") {
			double T = 0;
			std::istringstream(line.substr(pos + 1)) >> T;
			setTemperatureScale(T);
		}
		if (tag == "p") {
			double p = 0;
			std::istringstream(line.substr(pos + 1)) >> p;
			setPressureScale(p);
		}
		if (tag == "lambdaSq") {
			double lambdasq = 0;
			std::istringstream(line.substr(pos + 1)) >> lambdasq;
			setScaledDebyeLength(lambdasq);
		}
	}
	// TODO Show error message if input has unexpected format

	applyIdealGasLaw();
	compute();
	std::cout << *this << std::endl;
}

ScalingParameters::~ScalingParameters()
{
}

void ScalingParameters::applyIdealGasLaw()
{
	assert(n0 <= 0 && p0 > 0 && T0 > 0);
	const double kB = PhysicsConstants::instance().kB();
	n0 = p0 / (kB * T0);
	assert(n0 > 0 && p0 > 0 && T0 > 0);
}

/**
* Compute scaling parameter according to equation
* lambdaSq = const * T0 / (n0 * x0^2), where const = eps0 * kB / q^2.
*
* Assume that only one parameter is not given, identified by a non-positive value.
*/
void ScalingParameters::compute()
{

	assert(undefinedScales() == 1);
	PhysicsConstants & pc = PhysicsConstants::instance();
	const double eps0 = pc.eps0();
	const double kB = pc.kB();
	const double q = pc.q();
	const double a = eps0 * kB / pow(q, 2);

	if (x0 <= 0) {
		assert(lambdaSq > 0);
		assert(T0 > 0);
		assert(n0 > 0);
		x0 = a * sqrt(T0 / n0) / lambdaSq;
	}
	if (T0 <= 0) {
		assert(lambdaSq > 0);
		assert(x0 > 0);
		assert(n0 > 0);
		T0 = a * lambdaSq * n0 * pow(x0, 2);
	}
	if (n0 <= 0) {
		assert(lambdaSq > 0);
		assert(x0 > 0);
		assert(T0 > 0);
		n0 = a * T0 / (lambdaSq * pow(x0, 2));
	}
	if (lambdaSq <= 0) {
		assert(T0 > 0);
		assert(x0 > 0);
		assert(n0 > 0);
		lambdaSq = a * T0 / (n0 * pow(x0, 2));
	}
	assert(undefinedScales() == 0);
}

void ScalingParameters::setTemperatureScale(const double T)
{
	assert(T > 0);
	this->T0 = T;
}

void ScalingParameters::setNumberDensityScale(const double n)
{
	assert(n > 0);
	this->n0 = n;
}

void ScalingParameters::setLengthScale(const double x)
{
	assert(x > 0);
	this->x0 = x;
}

void ScalingParameters::setScaledDebyeLength(const double lambdasq)
{
	assert(lambdasq > 0);
	this->lambdaSq = lambdasq;
}

void ScalingParameters::setPressureScale(const double p)
{
	assert(p > 0);
	this->p0 = p;
}

double ScalingParameters::getTemperatureScale() const
{
	return T0;
}

double ScalingParameters::getLengthScale() const
{
	return x0;
}

double ScalingParameters::getNumberDensityScale() const
{
	return n0;
}

double ScalingParameters::getScaledDebyeLengthSquared() const
{
	return lambdaSq;
}

double ScalingParameters::getCrossSectionsScale() const
{
	return  1. / (getNumberDensityScale() * getLengthScale());
}

std::ostream & operator<<(std::ostream & os, const ScalingParameters & obj)
{
	os << "Scaling parameters: " << std::endl;
	os << "x: " << obj.x0 << std::endl;
	os << "n: " << obj.n0 << std::endl;
	os << "T: " << obj.T0 << std::endl;
	os << "lambdaSq: " << obj.lambdaSq;
	return os;
}