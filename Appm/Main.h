#pragma once

#include <iostream>
#include "AppmSolver.h"

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
};

