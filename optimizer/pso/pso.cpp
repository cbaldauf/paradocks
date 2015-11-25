#include "pso.hpp"


namespace paradocks
{

/// scale a rotation with scalar factor
/*inline Quat_t mul(const Quat_t& q, const Float_t f)
{
        Quat_t result;
        result.slerp(f,  Quat_t(), q);
        return result;
}*/

/*
/// catenate two rotations
inline Quat_t add(const Quat_t& q1, const Quat_t& q2)
{
        return q1*q2;
}

/// substract q2 from q1; the result is the rotation from q1 to q2
inline Quat_t sub(const Quat_t& q1, const Quat_t& q2)
{
        return q1*q2.conj();
}

PSO::Parameter PSO::initParameter(Float_t factor)
{
	boost::variate_generator<RNG&, boost::uniform_real<Float_t> > r(_rng, boost::uniform_real<Float_t>(-1,1));
	Parameter result;
	
	Float_t maxpos=0.1*factor;
	Float_t maxang=osg::PI/45*factor;
	do
	{
		result._pos=Vec3_t(r(), r(), r())*maxpos;
	} while (result._pos.length()>maxpos);
	
	boost::variate_generator<RNG&, boost::uniform_real<Float_t> > rsimple(_rng, boost::uniform_real<Float_t>(0,1));
	Float_t x0(rsimple());
	Float_t t1(2*osg::PI*rsimple()), t2(2*osg::PI*rsimple());
	Float_t s1(sin(t1)), s2(sin(t2)), c1(cos(t1)), c2(cos(t2));
	Float_t r1(sqrt(1-x0)), r2(sqrt(x0));
	Quat_t randomori=Quat_t(s1*r1, c1*r1, s2*r2, c2*r2);
	result._ori.slerp(maxang, Quat_t(), randomori);
	
	// bond rotatio
	result._bond=vector<Float_t>(_fit->RotCount());
	for (unsigned int i=0; i!=_fit->RotCount(); ++i)
	{
		result._bond[i]=r()*maxang;
	}
	return result;
}
*/

/*
PSO::Parameter PSO::add(PSO::Parameter x, PSO::Parameter y)
{
	Parameter result;
	result._pos=x._pos+y._pos;
	
	//quaternion multiplication is concatanation(addition) of the represented orientation
	result._ori=x._ori*y._ori;
	
	result._bond=x._bond;
	for (unsigned int i=0; i!=result._bond.size(); ++i)
	{
		result._bond[i]+=y._bond[i];
	}
	
	return result;
}

PSO::Parameter PSO::sub(PSO::Parameter x, PSO::Parameter y)
{
	Parameter result;
	result._pos=x._pos-y._pos;
	
	//quaternion multiplication is concatanation(addition) of the represented orientation
	//quaternion inversion results in a rotation in the opposite direction
	Quat_t tmp=y._ori;
	result._ori=x._ori*tmp.inverse();
	
	result._bond=x._bond;
	for (unsigned int i=0; i!=result._bond.size(); ++i)
	{
		result._bond[i]-=y._bond[i];
	}
	
	return result;
}
*/

template<typename T>
void PSO::constrict_to_PI(T& r)
{
	while (r<(-osg::PI))
	{
		r+=(2*osg::PI);
	}
	while (r>osg::PI)
	{
		r-=(2*osg::PI);
	}
}

struct less_particle : public binary_function<PSO::Particle,PSO::Particle,bool> {
bool operator()(PSO::Particle x, PSO::Particle y) { return x._score < y._score; }
};

/*
struct mpi_min_particle : public binary_function<PSO::Particle, PSO::Particle, PSO::Particle> {
PSO::Particle operator()(PSO::Particle x, PSO::Particle y) 
{ 
	if (x._score < y._score) return x;
	else return y;
}
};
*/

PSO::PSO(vector<string>& parameter, paradocks::Fitness* fit, RNG& rng, Vec3_t& center, Float_t& radius)
	: Optimizer(parameter, fit, rng, center, radius)
{
	if (parameter.size()!=6)
	{
		string emsg("Optimizer cpso expects 10 parameter");
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

	try
	{
		_particle_count=boost::lexical_cast<unsigned int>(parameter[1]);
	}
	catch(boost::bad_lexical_cast e)
	{
		string emsg("Error in input conversion! optimizer particle count line\n"+string(e.what()));
		throw OptimizerError(emsg);
	}
	
	try
	{
		_inertia_start=boost::lexical_cast<Float_t>(parameter[2]);
	}
	catch(boost::bad_lexical_cast e)
	{
		string emsg("Error in input conversion! optimizer inertia start line\n"+string(e.what()));
		throw OptimizerError(emsg);
	}
	
	try
	{
		_inertia_end=boost::lexical_cast<Float_t>(parameter[3]);
	}
	catch(boost::bad_lexical_cast e)
	{
		string emsg("Error in input conversion! optimizer inertia end line\n"+string(e.what()));
		throw OptimizerError(emsg);
	}
	
	try
	{
		_cognitive_weight=boost::lexical_cast<Float_t>(parameter[4]);
	}
	catch(boost::bad_lexical_cast e)
	{
		string emsg("Error in input conversion! optimizer cognitive weight line\n"+string(e.what()));
		throw OptimizerError(emsg);
	}
		
	try
	{
		_social_weight=boost::lexical_cast<Float_t>(parameter[5]);
	}
	catch(boost::bad_lexical_cast e)
	{
		string emsg("Error in input conversion! optimizer social weight line\n"+string(e.what()));
		throw OptimizerError(emsg);
	}
}

PSO::Particle PSO::init_new_Particle()
{
	Float_t max_v=_radius/4;
	Float_t max_a=osg::PI/4;
	
	boost::variate_generator<RNG&, boost::uniform_real<Float_t> > r(_rng, boost::uniform_real<Float_t>(-1,1));
	Particle result;
	// position
	do
	{
		result._pos=_center+Vec3_t(r() ,r(), r())*_radius;
	} while ((_center-result._pos).length()>_radius);
	result._pb=result._pos;
	// position velocity
	do
	{
		result._pv=Vec3_t(r(), r(), r())*max_v;
	} while (result._pv.length()>max_v);
	
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
	result._ob=result._ori;
	// orientation velocity
	x0=rsimple();
	t1=2*osg::PI*rsimple();
	t2=2*osg::PI*rsimple();
	s1=sin(t1);
	s2=sin(t2);
	c1=cos(t1);
	c2=cos(t2);
	r1=sqrt(1-x0);
	r2=sqrt(x0);
	Quat_t randomori=Quat_t(s1*r1, c1*r1, s2*r2, c2*r2);
	// max possible orientation velocity is pi(180)
	// but we allow only initanglevel orientation velocity
	result._ov.slerp(max_a/osg::PI, Quat_t(), randomori);
	
	// bond rotation
	result._bond=vector<Float_t>(_fit->RotCount());
	result._bv=vector<Float_t>(_fit->RotCount());
	for (unsigned int i=0; i!=_fit->RotCount(); ++i)
	{
		result._bond[i]=r()*osg::PI;
		// scale by maximum angle step width
		result._bv[i]=r()*max_a;
	}
	result._bb=result._bond;
	return result;
}

void PSO::run(ostream& outfile, unsigned int outiter, unsigned int outend)
{
	boost::variate_generator<RNG&, boost::uniform_real<Float_t> > r(_rng, boost::uniform_real<Float_t> (0,1));
	
	vector<Particle> allparticles;
	allparticles.reserve(_particle_count);
	
	for (unsigned int i=0; i!=_particle_count; ++i)
	{
		allparticles.push_back(init_new_Particle());
	}
	
	// first time evaluation
	for (vector<Particle>::iterator j=allparticles.begin(); j!=allparticles.end(); ++j)
	{
		vector<Float_t>* score;
		_fit->eval(j->_pos, j->_ori, j->_bond, score);
		j->_score=accumulate(score->begin(), score->end(), (Float_t)0);
		j->_bscore=j->_score;
	}
	Particle globalbest=*(min_element(allparticles.begin(), allparticles.end(), less_particle()));
	
	// main loop
	for (unsigned int i=0; i!=_iter; ++i)
	{
		// inertia weight
		Float_t w=_inertia_start-(_inertia_start-_inertia_end)*i/_iter;
		
		// particle loop
		for (vector<Particle>::iterator j=allparticles.begin(); j!=allparticles.end(); ++j)
		{
			// select neighborhood best
			Particle best=*(min_element(allparticles.begin(), allparticles.end(), less_particle()));
			
			Float_t cr=_cognitive_weight*r();
			Float_t sr=_social_weight*r();
			
			//// update velocity
			// position
			// apply inertia
			j->_pv*=w;
			// calculate velosity toward personal best
			Vec3_t cogvel((j->_pb-j->_pos)*cr);
			j->_pv+=cogvel;
			// calculate velosity toward global best
			Vec3_t socialvel((best._pos-j->_pos)*sr);
			j->_pv+=socialvel;
			
			// orientation
			// apply inertia
			Quat_t rvo;
			rvo.slerp(w, Quat_t(), j->_ov);
			// calculate velosity toward personal best
			Quat_t rvc;
			rvc.slerp(cr, j->_ori, j->_ob);
			rvc=rvc*j->_ori.conj();
			// calculate velosity toward global best
			Quat_t rvs;	
			rvs.slerp(sr, j->_ori, best._ori);
			rvs=rvs*j->_ori.conj();
			j->_ov=rvo*rvc*rvs;
			normalize(j->_ov);
			
			// bond rotation
			for (unsigned int k=0; k!=j->_bond.size(); ++k)
			{
				// apply inertia
				j->_bv[k]*=w;
				// calculate velosity toward personal best
				Float_t cogbv((j->_bb[k]-j->_bond[k])*cr);
				constrict_to_PI(cogbv);
				// calculate velosity toward global best
				Float_t socialbv((best._bond[k]-j->_bond[k])*sr);
				constrict_to_PI(socialbv);
				j->_bv[k]+=socialbv;
				j->_bv[k]+=cogbv;
			}
			
			// move
			j->_pos+=j->_pv;
			j->_ori*=j->_ov;
			normalize(j->_ori);
			for (unsigned int k=0; k!=j->_bond.size(); ++k)
			{
				j->_bond[k]+=j->_bv[k];
				while (j->_bond[k]<0)
				{
					j->_bond[k]+=(2*osg::PI);
				}
				while (j->_bond[k]>(2*osg::PI))
				{
					j->_bond[k]-=(2*osg::PI);
				}
			}
			
			// reset if outside
			if (j->_pos.x()<_center.x()-_radius || j->_pos.x()>_center.x()+_radius ||
				j->_pos.y()<_center.y()-_radius || j->_pos.y()>_center.y()+_radius ||
				j->_pos.z()<_center.z()-_radius || j->_pos.z()>_center.z()+_radius)
			{
				Particle p(init_new_Particle());
				//keep the particles best position
				j->_pos=p._pos;
				j->_pv=p._pv;
				j->_ori=p._ori;
				j->_ov=p._ov;
				j->_bond=p._bond;
				j->_bv=p._bv;
			}
		}
		
		// evaluation
		#pragma omp parallel for
		for (int j=0; j<allparticles.size(); ++j)
		{
			vector<Float_t>* score;
			_fit->eval(allparticles[j]._pos, allparticles[j]._ori, allparticles[j]._bond, score);
			allparticles[j]._score=accumulate(score->begin(), score->end(), (Float_t)0);
			// if score is better update personal best
			if (allparticles[j]._score<allparticles[j]._bscore)
			{
				allparticles[j]._bscore=allparticles[j]._score;
				allparticles[j]._pb=allparticles[j]._pos;
				allparticles[j]._ob=allparticles[j]._ori;
				allparticles[j]._bb=allparticles[j]._bond;
			}
		}
		#pragma omp barrier
		Particle globalmin=*(min_element(allparticles.begin(), allparticles.end(), less_particle()));
		
		if (outiter)
		{
			if (outiter==1)
			{
				outfile<<"#structure fitness =\""<<globalmin._score<<"\""<<endl;
				_fit->writeOptimum(globalmin._pos, globalmin._ori, globalmin._bond, outfile);
			}
			if (outiter==2)
			{
				for (vector<Particle>::iterator k=allparticles.begin(); k!=allparticles.end(); ++k)
				{
					outfile<<"#structure fitness =\""<<k->_score<<"\""<<endl;
					_fit->writeOptimum(k->_pos, k->_ori, k->_bond, outfile);
				}
			}
		}
		
		// select neighborhood best
		if (globalmin._score < globalbest._score)
		{
			globalbest=globalmin;
			std::cout <<i<<" new global best "<< globalbest._score<<endl;
		}
	}
	
	// refine by local minimization
	/*
	Parameter current_state;
	current_state._pos=globalbest._pos;
	current_state._ori=globalbest._ori;
	current_state._bond=globalbest._bond;
	vector<Float_t>* result;
	_fit->eval(current_state._pos, current_state._ori, current_state._bond, result);
	cout << "offset "<<globalbest._pos.x()<<" "<< globalbest._pos.y()<<" "<< globalbest._pos.z()<<endl;
	Float_t bestscore=accumulate(result->begin(), result->end(), (Float_t)0);
	cout << "SCORE "<<bestscore<<endl;
	Float_t scaling=4;
	
	for (unsigned int i=0; i!=10000; ++i)
	{
		Parameter dev=initParameter(scaling);
		Parameter new_pos=add(current_state, dev);
	
		vector<Float_t>* result;
		_fit->eval(new_pos._pos, new_pos._ori, new_pos._bond, result);
		Float_t new_score=accumulate(result->begin(), result->end(), (Float_t)0);

		if (new_score<bestscore)
		{
			bestscore=new_score;
			current_state=new_pos;
			cout << "SCORE "<<bestscore<<endl;
		}
		else
		{
			new_pos=sub(current_state, dev);
			vector<Float_t>* result;
			_fit->eval(new_pos._pos, new_pos._ori, new_pos._bond, result);
			Float_t new_score=accumulate(result->begin(), result->end(), (Float_t)0);
			if (new_score<bestscore)
			{
				bestscore=new_score;
				current_state=new_pos;
				cout << "SCORE "<<bestscore<<endl;
			}
			else
			{
				scaling*=0.95;
			}
		}
	}
	*/

	// output best
	if (outend==1)
	{
		outfile<<"#structure fitness =\""<<globalbest._score<<"\""<<endl;
		_fit->writeOptimum(globalbest._pos, globalbest._ori, globalbest._bond, outfile);
	}
	//output all
	if (outend==2)
	{
		for (vector<Particle>::iterator j=allparticles.begin(); j!=allparticles.end(); ++j)
		{
			outfile<<"#structure fitness =\""<<j->_bscore<<"\""<<endl;
			_fit->writeOptimum(j->_pb, j->_ob, j->_bb, outfile);
		}
	}
	
	_best=globalbest;
}

}
