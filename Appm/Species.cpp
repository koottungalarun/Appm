#include "Species.h"



Species::Species()
{
}

Species::Species(const std::string & symbol, const std::string & name, const double massRatio, const int charge)
{
	this->symbol = symbol;
	this->name = name;
	this->massRatio = massRatio;
	this->charge = charge;
}

Species::Species(const std::string & str)
{
	const char delim = ',';
	std::string temp = str;
	std::vector<std::string> chunks;
	int counter = 0; // protect against infinite loop
	while (temp.size() > 0 && counter < 10) {
		const int pos = temp.find(delim);
		std::string chunk = temp.substr(0, pos);
		chunks.push_back(chunk);
		if (pos == std::string::npos) { break; } // delimiter not found
		temp = temp.substr(pos + 1);
		counter++;
	}
	assert(chunks.size() >= 4);

	const std::string symbol = trim(chunks[0]);
	const std::string name = trim(chunks[1]);
	const double massRatio = std::stod(chunks[2]); // string-to-double
	const int charge = std::stoi(chunks[3]);       // string-to-int
	
	assert(symbol.size() > 0);
	this->symbol = symbol;
	assert(name.size() > 0);
	this->name = name;
	assert(massRatio > 0);
	this->massRatio = massRatio;
	this->charge = charge;
}


Species::~Species()
{
}

const int Species::getCharge() const
{
	return this->charge;
}

const double Species::getMassRatio() const
{
	return this->massRatio;
}

const std::string Species::getName() const
{
	return this->name;
}

std::ostream & operator<<(std::ostream & os, const Species & obj)
{
	os << "name:       " << obj.name <<  std::endl;
	os << "symbol:     " << obj.symbol << std::endl;
	os << "mass ratio: " << obj.massRatio << std::endl;
	os << "charge:     " << obj.charge;
	return os;
}

