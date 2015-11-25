#ifndef RANDOM_HILL_CLIMBING_H
#define RANDOM_HILL_CLIMBING_H

#include "optimizer/optimizer.hpp"
#include "mgraph/mgraphtypes.hpp"

namespace paradocks
{
/**
\brief random hill climbing optimization class.
\verbatim
<optimizer type="random">  <!-- name of the optimizer-->
<par val="1000"/>       <!-- number of iterations -->
<par val="1000"/>       <!-- number of iterations in local optimization -->
</optimizer>
\endverbatim
\ingroup optimizer
\ingroup framework
**/
class RANDOM_HILL_CLIMBING : public Optimizer
{
public:
	struct Parameter
	{
		Vec3_t _pos; // position
		Quat_t _ori; // orientation
		vector<Float_t> _bond; // rotatable bonds
		Float_t _score; // score
	};
	
	RANDOM_HILL_CLIMBING(vector<string>& parameter, paradocks::Fitness* fit, RNG& rng, Vec3_t& center, Float_t& radius);
	
	void run(ostream& output, unsigned int outiter, unsigned int outend);

protected:
	// derived from base class
	// parameter for the optimizer
	// vector<string>& _parameter;
	// fitness object
	// Fitness* _fit;
	// random number generator
	// RNG _rng;
	
	// active site
	// Vec3_t& _center;
	// Float_t& _radius;
	
	unsigned int _iter;
	unsigned int _local_iter;
	
	Parameter _best;
	Parameter move(Parameter);
	Parameter init_new_Parameter();
};

}

#endif /*RANDOM_H*/
