#include "random.hpp"


namespace paradocks
{

RANDOM::RANDOM(vector<string>& parameter, paradocks::Fitness* fit, RNG& rng, Vec3_t& center, Float_t& radius)
	: Optimizer(parameter, fit, rng, center, radius)
{
	if (parameter.size()!=1)
	{
		string emsg("Optimizer random expects 1 parameter");
		throw OptimizerError(emsg);
	}
	try
	{
		_iter=boost::lexical_cast<unsigned int>(parameter[0]);
	}
	catch(boost::bad_lexical_cast e)
	{
		string emsg("Error in input conversion! optimizer iteration line\n"+string(e.what()));
		throw OptimizerError(emsg);
	}

}

RANDOM::Parameter RANDOM::init_new_Parameter()
{
	boost::variate_generator<RNG&, boost::uniform_real<Float_t> > r(_rng, boost::uniform_real<Float_t>(-1,1));
	Parameter result;
	// position
	do
	{
		result._pos=_center+Vec3_t(r() ,r(), r())*_radius;
	} while ((_center-result._pos).length()>_radius);
	
	// orientation
	// random quaternion
	// Ken Shoemake
	// UNIFORM RANDOM ROTATIONS
	// In D. Kirk, editor, Graphics Gems III, pages 124-132. Academic, New York, 1992.
	boost::variate_generator<RNG&, boost::uniform_real<Float_t> > 
	rsimple(_rng, boost::uniform_real<Float_t>(0,1));
	Float_t x0(rsimple());
	Float_t t1(2*osg::PI*rsimple()), t2(2*osg::PI*rsimple());
	Float_t s1(sin(t1)), s2(sin(t2)), c1(cos(t1)), c2(cos(t2));
	Float_t r1(sqrt(1-x0)), r2(sqrt(x0));
	result._ori=Quat_t(s1*r1, c1*r1, s2*r2, c2*r2);
	
	// bond rotation
	result._bond=vector<Float_t>(_fit->RotCount());
	for (unsigned int i=0; i!=_fit->RotCount(); ++i)
	{
		result._bond[i]=r()*osg::PI;
	}
	return result;
}

void RANDOM::run(ostream& outfile, unsigned int outiter, unsigned int outend)
{
	_best=init_new_Parameter();
	vector<Float_t>* score;
	_fit->eval(_best._pos, _best._ori, _best._bond, score);
	_best._score=accumulate(score->begin(), score->end(), (Float_t)0);
	for (unsigned int i=0; i!=_iter; ++i)
	{
		Parameter np=init_new_Parameter();
		_fit->eval(np._pos, np._ori, np._bond, score);
		np._score=accumulate(score->begin(), score->end(), (Float_t)0);
		if (np._score<_best._score)
		{
			_best=np;
			std::cout <<i<<" new best "<< _best._score<<endl;
		}
	}
	
	// output best
	if (outend==1 || outend==2)
	{
		outfile<<"#structure fitness =\""<<_best._score<<"\""<<endl;
		_fit->writeOptimum(_best._pos, _best._ori, _best._bond, outfile);
	}
}

}
