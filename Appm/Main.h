#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "AppmSolver.h"
#include "Species.h"


class Main
{
public:
	Main();
	~Main();

	void run();

	void setInputFilename(const std::string & filename);
	void readInputFile();

private:
	std::string showVersionInfo();

	std::string inputFilename;
	std::string meshFilename;

	std::vector<Species> speciesList;
	//std::vector<ElasticCollision> elasticCollisions;
	std::vector<std::string> elasticCollisionList;

	AppmSolver::SolverParameters solverParameters;

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

