/*
 * HonneyBee.h
 *
 *  Created on: 24.02.2010
 *      Author: daniel
 */

#ifndef HONEYBEE_H_
#define HONEYBEE_H_

#include "optimizer/optimizer.hpp"
#include "mgraph/mgraphtypes.hpp"
#include <stdio.h>
namespace paradocks
{
class HoneyBee : public Optimizer
{
public:

	struct Elem{
		Vec3_t ort;
		Vec3_t dreh_x;
		Quat_t dreh;
		vector<Float_t> bond;
		Float_t score;

		Elem(){

		}
		Elem(const Elem&e){
			ort=e.ort;
			dreh=e.dreh;
			bond=e.bond;
			score=e.score;
		}
		void print(FILE*file=stderr){
		    fprintf(file,"(%f, %f, %f),(%f, %f),(%f, %f, %f, %f)",
			ort._v[0],ort._v[1],ort._v[2],
			dreh_x._v[0],dreh_x._v[1],
			dreh._v[0],dreh._v[1],dreh._v[2],dreh._v[3]
			);
		    for(size_t k=0;k<bond.size();k++){
			fprintf(file,", %f",bond[k]);
		    }
		    fprintf(stderr,"\n");
		}
	};

	HoneyBee(vector<string>& parameter, paradocks::Fitness* fit, RNG& rng, Vec3_t& center, Float_t& radius);
	void run(ostream& output, unsigned int outiter, unsigned int outend);
	virtual ~HoneyBee();
	Elem best(){return _best;}
	void setFollowCount(int fc);
	int getFollowCount();
	unsigned int zwischenschritte;
protected:
	Elem _best;
	unsigned int _iter;
	unsigned int _particle_count;
	int ScoutFollowers;
	int NoRandomEval;
	int FollowCount;
	Float_t _inertia_start;
	Float_t _inertia_end;

	Float_t _cognitive_weight;
	Float_t _social_weight;


	Float_t hbPosRange,hbDreRange,hbBondRange;
	Float_t hbfPosRange,hbfDreRange,hbfBondRange;


	template<typename T>
		void constrict_to_PI(T& r);
	template<typename T>
		void constrict_to_2PI(T& r);
	boost::variate_generator<RNG&, boost::uniform_real<Float_t> >
			r_0_1;
	boost::variate_generator<RNG&, boost::uniform_real<Float_t> >
			r_1_0_1;
	Float_t degDiff(Float_t d1,Float_t d2);
	Float_t r01();
	Float_t r101();
	Float_t _run(ostream& outfile, unsigned int outiter, unsigned int outend);
	void _eval(Elem& e);
	void _route(const Elem& nest,Elem&e);
	void _calculateDrehQuat(Elem & e);
	void _rnd(Elem& e,Vec3_t,Float_t);
	void _rnd(Elem& e,const Elem &Center,Float_t radius,Float_t radDreh,Float_t radBond);
	unsigned int GetRandomArrayElement(Float_t*array,int size);
	Float_t _runs(ostream& outfile, unsigned int outiter, unsigned int outend,unsigned int runs);
	int evalCalls;

};
class RandomLocalSearch : public HoneyBee
{
public:

	RandomLocalSearch(vector<string>& parameter, paradocks::Fitness* fit, RNG& rng, Vec3_t& center, Float_t& radius);
	void run(std::ostream& output, unsigned int outiter, unsigned int outend);

protected:
};

}
#endif /* HONEYBEE_H_ */
