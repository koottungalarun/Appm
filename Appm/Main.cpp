#include "Main.h"

int main() {
	std::cout << "***********************" << std::endl;
	std::cout << "*    APPM             *" << std::endl;
	std::cout << "***********************" << std::endl;
	Main main;
	main.run();
	std::cout << "TERMINATED" << std::endl;
	return EXIT_SUCCESS;
}


Main::Main()
{	
	std::cout << showVersionInfo() << std::endl;
	std::cout << std::endl;;
}


Main::~Main()
{
}

void Main::run()
{
	PrimalMesh::PrimalMeshParams primalMeshParams;
	primalMeshParams = PrimalMesh::PrimalMeshParams("primalMeshParams.txt");

	AppmSolver appm(primalMeshParams);
	appm.run();	
}

/**
* @return a text showing compile time and versions of used libraries.
*/
std::string Main::showVersionInfo()
{
	std::stringstream ss;
	ss << "Compiled on " << __TIMESTAMP__ << std::endl;
	ss << "with libraries: " << std::endl;
	
	// Version of Eigen library
	ss << Eigen::showVersion() << std::endl;

	// Version of MKL library
	const int len = 200;
	char mklVersion[len];
	mkl_get_version_string(mklVersion, len);
	ss << mklVersion << std::endl;

	// Version of HDF5 library
	ss << H5_VERS_INFO << std::endl;
	return ss.str();
}
