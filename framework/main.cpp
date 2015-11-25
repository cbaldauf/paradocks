#include <string>
#include <iostream>
#include <stdexcept>
#include "Paradocks.hpp"

using namespace std;
using namespace paradocks;

int main(int argc, char *argv[])
{
	cout<<"ParaDockS 0.1.2"<<endl;
	
	// get config file content
	ParadocksConfig cfg;
	// on rank 0 parse the configfile
	try
	{
		cfg=ParadocksConfig(((argc<=1) ? "paradocks.xml" : argv[1]));
	}
	catch(std::runtime_error e)
	{
		cout<<e.what()<<endl;
	}
	
	// set up application
	Paradocks paradocks;
	try
	{
		paradocks=Paradocks(&cfg);
	}
	catch(std::runtime_error e)
	{
		cout<<e.what()<<endl;
	}
	
	// execute
	try
	{
		paradocks.run();
	}
	catch(std::runtime_error e)
	{
		cout<<e.what()<<endl;
	}
	
	return 0;
}
