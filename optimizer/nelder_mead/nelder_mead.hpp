#ifndef NELDER_MEAD_H
#define NELDER_MEAD_H

#include "optimizer/optimizer.hpp"
#include "mgraph/mgraphtypes.hpp"
//#include "MinFinder.h"
//#include "Point.h"

namespace paradocks
{
	
class NELDER_MEAD : public Optimizer
{
    public:
		struct Parameter
		{
			Vec3_t _pos; // position
			Quat_t _ori; // orientation
			vector<Float_t> _bond; // rotatable bonds
			Float_t _score; // score
		};
		
		NELDER_MEAD(vector<string>& parameter, paradocks::Fitness* fit, RNG& rng, Vec3_t& center, Float_t& radius);
		
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
	Parameter init_new_Parameter(Vec3_t);
	
/*        virtual Point *findMinima(Point *p = NULL);
#ifdef TEST
        bool test();
#endif

private:
        Point *calcCentroid(vector<Point *> &);
        void updateCentroid(Point &, Point &, Point &);
        Point *reflect(Point &, Point &);
        Point *expand(Point &, Point &);
        Point *contract(Point &, Point &);
        void reduce(vector<Point *> &);
        double stdError(vector<Point *> &);*/
};



}

#endif
