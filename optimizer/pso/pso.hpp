#ifndef PSO_H
#define PSO_H

#include "optimizer/optimizer.hpp"
#include "mgraph/mgraphtypes.hpp"

namespace paradocks
{
/**
\brief Particle Swarm Optimization class.
This class implements particle swarm optimization from:\n
J. Kennedy and R. Eberhart \n
<em>Particle swarm optimization</em> \n
Proc. of the IEEE Int. Conf. on Neural Networks, Piscataway, NJ, pp. 1942-1948, 1995. \n \n
an example input section looks like this:
\verbatim
<optimizer type="pso">  <!-- name of the optimizer-->
<par val="100000"/>     <!-- iterations -->
<par val="30"/>         <!-- particle count -->
<par val="1"/>          <!-- inertia weight start -->
<par val="0.7"/>        <!-- inertia weight end -->
<par val="1.4"/>        <!-- cognitive weight -->
<par val="1.4"/>        <!-- social weight -->
<par val="5"/>          <!-- maximum velocity -->
<par val="0.79"/>       <!-- maximum angle velocity for quaternions -->
<par val="0.79"/>       <!-- maximum angle velocity for angles -->
</optimizer>
\endverbatim
\ingroup optimizer
\ingroup framework
**/
class PSO : public Optimizer
{
public:
	struct Particle
	{
		Vec3_t _pos, _pv, _pb;
		Quat_t _ori, _ov, _ob;
		vector<Float_t> _bond, _bv, _bb;
		Float_t _score, _bscore;
	};
	
	PSO(vector<string>& parameter, paradocks::Fitness* fit, RNG& rng, Vec3_t& center, Float_t& radius);
	
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
	unsigned int _particle_count;
	
	Float_t _inertia_start;
	Float_t _inertia_end;
	
	Float_t _cognitive_weight;
	Float_t _social_weight;
	
	Particle _best;
	
	template<typename T>
	void constrict_to_PI(T& r);
	Particle init_new_Particle();
	void write_Particle(Particle& p);
	Particle best() {return _best;}
	
	struct Parameter
	{
		Vec3_t _pos;
		Quat_t _ori;
		vector<Float_t> _bond;
	};
	Parameter initParameter(Float_t factor);
	Parameter add(Parameter x, Parameter y);
	Parameter sub(Parameter x, Parameter y);
};

}

#endif /*PSO_H*/
