#ifndef PARADOCKS_H
#define PARADOCKS_H

#include <boost/shared_ptr.hpp>

#include "ParadocksConfig.hpp"
#include "fitness/fitness.hpp"
#include "optimizer/optimizer.hpp"

namespace paradocks
{
using namespace std;

class ParadocksError : public runtime_error
{
public:
	ParadocksError(string err) : runtime_error("ParadocksError: "+err) {}
};

class Paradocks
{
public:
	Paradocks() : _cfg(NULL) {}
	Paradocks(ParadocksConfig* cfg);
	virtual ~Paradocks() {}
	
	void run();
	
protected:
	ParadocksConfig* _cfg;
	vector<unsigned int> _ligandIdx;
	
	boost::shared_ptr<Fitness> _fit;
	boost::shared_ptr<Optimizer> _opt;
	Vec3_t _as;
	Float_t _rad;
	unsigned int _runs;
	RNG _rng;
	unsigned int _outiter;
	unsigned int _outend;
	string _outprefix;	
};
}
#endif /*PARADOCKS_H*/
