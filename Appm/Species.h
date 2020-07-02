#pragma once

#include <cassert>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>

class Species
{
public:
	Species();
	Species(const std::string & symbol, const std::string & name, const double massRatio, const int charge);
	Species(const std::string & str);
	~Species();

	friend std::ostream & operator<<(std::ostream & os, const Species & obj);

	const int getCharge() const;
	const double getMassRatio() const;

private:
	std::string symbol;
	std::string name;
	double massRatio = 0;
	int charge = 0;

	/**
	* Trim white space of a string. 
	* @see https://stackoverflow.com/a/17976541
	*/
	inline std::string trim(const std::string &s)
	{
		auto  wsfront = std::find_if_not(s.begin(), s.end(), [](int c) {return std::isspace(c); });
		return std::string(wsfront, std::find_if_not(s.rbegin(), std::string::const_reverse_iterator(wsfront), [](int c) {return std::isspace(c); }).base());
	}
};



