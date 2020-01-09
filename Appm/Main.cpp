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
}


Main::~Main()
{
}

void Main::run()
{
	AppmSolver appm;
	appm.run();	
}
