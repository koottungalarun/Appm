#include "H5Writer.h"



H5Writer::H5Writer()
{
}

H5Writer::H5Writer(const std::string & filename)
{
	std::cout << "Create HDF5 file: " << filename << std::endl;
	const int length = 3;
	assert(filename.size() > 3);
	const int pos = filename.size() - 3;
	std::string substr = filename.substr(pos, length);
	assert(filename.compare(pos, length, ".h5") == 0);

	this->filename = filename;
	this->file = H5File(filename.c_str(), H5F_ACC_TRUNC);
}

H5Writer::H5Writer(const std::string & filename, const bool showOutput)
	: H5Writer(filename)
{
	this->showOutputInfo = showOutput;
}

H5Writer::~H5Writer()
{
}

void H5Writer::setShowOutput(const bool isShowOutput) 
{
	this->showOutputInfo = isShowOutput;
}
