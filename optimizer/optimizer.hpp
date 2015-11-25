#ifndef OPTIMZER_H
#define OPTIMZER_H

#include <boost/random.hpp>
#include "fitness/fitness.hpp"


namespace paradocks
{
using namespace std;

typedef boost::mt19937 RNG;

/**
\defgroup optimizer ParaDockS Optimizer function Interface 
**/

/**
\brief Optimizer error class.
This class should be thrown by all implemented Optimizer classes in case of an error.
\ingroup optimizer
\ingroup framework
**/
class OptimizerError : public runtime_error
{
public:
	OptimizerError(string err) : runtime_error("OptimizerError: "+err) { }
};

/**
\brief Optimizer class.
Virtual base class interface definition for all Optimizer class implementations.
\ingroup optimizers
\ingroup framework
**/
class Optimizer
{
public:
	/// Constructor
	Optimizer(vector<string>& parameter, paradocks::Fitness* fit, RNG& rng, Vec3_t& center, Float_t& radius)
		: _parameter(parameter), _fit(fit), _rng(rng), _center(center), _radius(radius) {}
	/// Destructor
	virtual ~Optimizer() {}
	
	/**
	\brief Start the optimization run.
	\param output The results should go into this stream.
	\param outiter Save the structure from every iteration: 0 - no structure, 1 - best structure, 2 - all structures
	\param outend Save the structure at the end: 1 - best structure, 2 - all structures
	**/
	/// start a optimization
	virtual void run(ostream& output, unsigned int outiter, unsigned int outend) = 0;
	
	
protected:
	// parameter for the optimizer
	vector<string>& _parameter;
	
	// fitness object
	Fitness* _fit;
	// random number generator
	RNG& _rng;
	
	// active site
	Vec3_t& _center;
	Float_t& _radius;
};

}

#endif /* OPTIMZER_H */
