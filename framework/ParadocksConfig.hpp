#ifndef PARADOCKSCONFIG_H
#define PARADOCKSCONFIG_H

#include <string>
#include <vector>
#include <stdexcept>

namespace paradocks
{
using namespace std;

class ParadocksConfigError : public runtime_error
{
public:
	ParadocksConfigError(string err) : runtime_error("ParadocksConfigError: "+err) { }
};

/**
\brief class to hold configuration data.
This class parses the configuration file and holds the configuration data.
With the help of this object, data is transmitted between participating nodes.
\ingroup framework
**/
class ParadocksConfig
{
public:
	ParadocksConfig() {}
	ParadocksConfig(string configfile);
	virtual ~ParadocksConfig() {}
	
	string& X() { return _x; }
	string& Y() { return _y; }
	string& Z() { return _z; }
	string& Rad() { return _rad; }
	string& Seed() { return _seed; }
	string& Runs() { return _runs; }
	string& Protein() { return _protein; }
	string& FitnessType() { return _fitnesstype; }
	vector<string>& FitnessPar() { return _fitnesspar; }
	string& OptimizerType() { return _optimizertype; }
	vector<string>& OptimizerPar() { return _optimizerpar; }
	vector<string>& Ligand() { return _ligand; }
	vector<string>& LigandIdx() { return _idx; }
	string& OutIteration() { return _iterout; }
	string& OutEnd() { return _endout; }
	string& OutPrefix() { return _prefixout; }
	
protected:
	string _optimizertype, _fitnesstype;
	vector<string> _optimizerpar, _fitnesspar;
	string _protein, _x, _y, _z, _rad;
	vector<string> _ligand, _idx;
	string _seed;
	string _runs;
	string _iterout;
	string _endout;
	string _prefixout;
};

}

#endif /*PARADOCKSONFIG_H*/
