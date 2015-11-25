/*
 * HoneyBee.cpp
 *
 *  Created on: 24.02.2010
 *      Author: daniel
 */

#include "HoneyBee.hpp"
#include <time.h>
#include <stdio.h>
namespace paradocks{
const double _1PI2=0.5/osg::PI;
const double _2PI=osg::PI*2;
//HoneyBee Algo genauer als PSO Algo
template<typename T>
void HoneyBee::constrict_to_PI(T& r)
{
	if(r<(-osg::PI)){
		int f=(r-osg::PI)*_1PI2;
		r-=f*_2PI;
	}else if(r>(osg::PI)){
		int f=(r+osg::PI)*_1PI2;
		r-=f*_2PI;
	}
}
template<typename T>
void HoneyBee::constrict_to_2PI(T& r)
{
	if(r<0){
		int f=(r-_2PI)*_1PI2;
		r-=f*_2PI;
	}else if(r>_2PI){
		int f=(r+osg::PI)*_1PI2;
		r-=f*_2PI;
	}
}
Float_t HoneyBee::degDiff(Float_t d1,Float_t d2){
	Float_t dd(d2-d1);
	constrict_to_2PI(dd);
	return dd;
}

Float_t HoneyBee::r01(){
	return r_0_1();
}
Float_t HoneyBee::r101(){
	return r_1_0_1();
}

unsigned int HoneyBee::GetRandomArrayElement(Float_t*array,int size){
	Float_t follower=r01()*array[size-1];
	if(follower<array[0]) return 0;
	unsigned int i=1;
	unsigned int j=size-1;
	unsigned int p=i;
	while(i<=j){
		p=(i+j)/2;
		if( follower<array[p-1])
			j=p;
		else if(follower>array[p])
			i=p+1;
		else
			break;
	}
	return p;
}
HoneyBee::HoneyBee(vector<string>& parameter, paradocks::Fitness* fit, RNG& rng, Vec3_t& center, Float_t& radius)
: Optimizer(parameter, fit, rng, center, radius),
  ScoutFollowers(0),
  NoRandomEval(1),
  FollowCount(0),
  r_0_1(_rng, boost::uniform_real<Float_t>(0,1)),
  r_1_0_1(_rng, boost::uniform_real<Float_t>(-1,1)),
  evalCalls(0)
  {
	zwischenschritte=1;
	if (parameter.size()<2)
	{
		string emsg("Optimizer HBee expects min 2 parameter\n"
					"1: Iterationcount\n"
					"2: particle count\n");
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
	for(unsigned int k=2;k<parameter.size();k++){
		if(parameter[k]=="NoRandomEval") NoRandomEval=1;
		if(parameter[k]=="ScoutFollowers")ScoutFollowers=1;
		if(parameter[k]=="optimizePara") NoRandomEval=2;
	}
}
void HoneyBee::setFollowCount(int fc){
	FollowCount=fc;
}
int HoneyBee::getFollowCount(){
	return FollowCount;
}
void HoneyBee::_eval(Elem& e){
	vector<Float_t>* result;
	_fit->eval(e.ort,e.dreh,e.bond,result);
	#pragma omp atomic
	evalCalls++;
	e.score=accumulate(result->begin(),result->end(),Float_t(0) );
}
void HoneyBee::_route(const Elem& nest,Elem&e){
	if(zwischenschritte==0)return;
	unsigned int steps=zwischenschritte+1;

	Elem step(e);   //schrittweite
	Elem strt(nest);//Startposition
	step.ort         = (e.ort-nest.ort      )          / steps;
	step.dreh_x      = (e.dreh_x-nest.dreh_x)          / steps;
	for(unsigned int k = 0;k!=e.bond.size();k++)
		step.bond[k] = degDiff(nest.bond[k],e.bond[k]) / steps;

	for(unsigned int k=0;k!=zwischenschritte;k++){
		strt.ort+=step.ort;
		strt.dreh_x+=step.dreh_x;
		for(unsigned int l=0;l!=e.bond.size();l++)
			strt.bond[l]+=step.bond[l];
		_calculateDrehQuat(strt);
		_eval(strt);
		if (strt.score>e.score)e=strt;
	}
}
void HoneyBee::_calculateDrehQuat(Elem & e){
	Float_t t1(2*osg::PI*e.dreh_x._v[0]),
			t2(2*osg::PI*e.dreh_x._v[1]);
	Float_t s1(sin(t1)),
			s2(sin(t2)),
			c1(cos(t1)),
			c2(cos(t2));
	Float_t r1(sqrt(1-e.dreh_x._v[2])),
			r2(sqrt(e.dreh_x._v[2]));
	e.dreh=Quat_t(s1*r1, c1*r1, s2*r2, c2*r2);
	normalize(e.dreh);
}
void HoneyBee::_rnd(Elem &e,Vec3_t center,Float_t radius){
	Float_t radius2=radius*radius;
	// position

	do
	{
		e.ort=center+Vec3_t(r_1_0_1() ,r_1_0_1(), r_1_0_1())*radius;
	} while ((center-e.ort).length2()>radius2);
	// orientation
	// random quaternion
	// Ken Shoemake
	// UNIFORM RANDOM ROTATIONS
	// In D. Kirk, editor, Graphics Gems III, pages 124-132. Academic, New York, 1992.
//	boost::variate_generator<RNG&, boost::uniform_real<Float_t> >
//	rsimple(_rng, boost::uniform_real<Float_t>(0,1));
	e.dreh_x=Vec3_t(r_0_1() ,r_0_1(), r_0_1());
	Float_t t1(2*osg::PI*e.dreh_x._v[0]),
			t2(2*osg::PI*e.dreh_x._v[1]);
	Float_t s1(sin(t1)),
			s2(sin(t2)),
			c1(cos(t1)),
			c2(cos(t2));
	Float_t r1(sqrt(1-e.dreh_x._v[2])),
			r2(sqrt(e.dreh_x._v[2]));
	e.dreh=Quat_t(s1*r1, c1*r1, s2*r2, c2*r2);

	// bond rotation
	if(e.bond.size()!=_fit->RotCount())
		e.bond=vector<Float_t>(_fit->RotCount());

	for (unsigned int i=0; i!=e.bond.size(); ++i)
	{
		e.bond[i]=r_0_1()*2*osg::PI;
		constrict_to_2PI(e.bond[i]);
	}
}
void HoneyBee::_rnd(Elem& e,const Elem &Center,Float_t radOrt,Float_t radDreh,Float_t radBond){
	Float_t _radius2=_radius*_radius;
//	do
//	{
//		e.ort=Center.ort+Vec3_t(r_1_0_1() ,r_1_0_1(), r_1_0_1())*radius;
//		//Radius einhalten und
//	} while ((Center.ort-e.ort).length2()>radius2 || e.ort.length2()>_radius2);
	Vec3_t tOrt;
	int exitcounter=0;
	if((Center.ort-_center).length2()>_radius2){
		e.ort=(_center-Center.ort);
		e.ort/=e.ort.length();
		e.ort+=_center;
	}else{
		do{
			do{
				tOrt=Vec3_t(r_1_0_1() ,r_1_0_1(), r_1_0_1());
				exitcounter++;
			}while (tOrt.length2()>1);
			e.ort=Center.ort+tOrt*radOrt;
			if (exitcounter>100){
			    fprintf(stderr,"massive problems generating acceptable vector\n");
			    break;
			}
		}while((_center-e.ort).length2()>_radius2);
	}
	Vec3_t tDreh;
	do{
		tDreh=Vec3_t(r101(),r101(),r101());

	}while(tDreh.length2()>1);
	e.dreh_x=Center.dreh_x+tDreh*radDreh;
	Float_t t1(2*osg::PI*e.dreh_x._v[0]),
			t2(2*osg::PI*e.dreh_x._v[1]);
	constrict_to_2PI(t1);
	constrict_to_2PI(t2);

	Float_t s1(sin(t1)),
			s2(sin(t2)),
			c1(cos(t1)),
			c2(cos(t2));
	Float_t r1(sqrt(1-e.dreh_x._v[2])),
			r2(sqrt(e.dreh_x._v[2]));
	e.dreh=Quat_t(s1*r1, c1*r1, s2*r2, c2*r2);
	normalize(e.dreh);
	if(e.bond.size()!=_fit->RotCount())
		e.bond=vector<Float_t>(_fit->RotCount());
	for(unsigned int k=0;k<e.bond.size();k++){
		e.bond[k]=Center.bond[k]+radBond*r101();
	}

}
Float_t HoneyBee::_run(ostream& outfile, unsigned int outiter, unsigned int outend){
	unsigned int cc=0;
	//Arrays für die Scouts und Follower
	size_t scout_size(_particle_count/3);
	size_t follo_size(_particle_count-scout_size);

	Elem   *scout    = new Elem  [scout_size];
	size_t *scoutcho = new size_t[scout_size];
	Elem   *follower = new Elem  [follo_size];
	//Nest Standort festlegen
	Elem nest;
	_rnd(nest,_center,_radius);
	_eval(nest);
	//Best und AllBest festlegen
	Elem best(nest);
	Elem ALLbest(nest);
	Float_t *Summi=new Float_t[scout_size];
	Summi[0]=0;
	unsigned int counts=0;
//	int printi=0;
	Float_t faktor;
	while(evalCalls<_iter*_particle_count){
		counts++;
		faktor=1.0-counts/_iter;
//		//scouts
		for (unsigned int k=0;k<scout_size;k++){
			scoutcho[k]=0;
			if(k)
			    Summi[k]=Summi[k-1];
			//for(int l=0;l<5;l++){
				_rnd(scout[k],nest,hbPosRange*faktor,faktor*hbDreRange,faktor*hbBondRange);
				//no local opt today ...
				_eval(scout[k]);
				float_t ffff(nest.score*0.95-scout[k].score);
				if(ffff>0){
					//negative energiewerte -> Faktor kleiner 1
					Summi[k]+=ffff;
				//	break;
				}
			//}
////
			if(scout[k].score <best.score){
				best=scout[k];
			}
		}
//		//follower
//		//choose scout to follow
#pragma omp parallel for
		for(unsigned int k=0;k<follo_size;k++){
			unsigned int foll=
					GetRandomArrayElement(Summi,scout_size);

			_rnd(follower[k],scout[foll],faktor*hbfPosRange,faktor*hbfDreRange,faktor*hbfBondRange);
			_eval(follower[k]);
			_route(scout[foll],follower[k]);
#pragma omp atomic
			scoutcho[foll]++;
#pragma omp critical
			if(follower[k].score <best.score)
				best=follower[k];
#pragma omp end critical
		}
		cc++;

		if(best.score<nest.score){
			//local optimization
			unsigned int MaxLO=4;
			for(unsigned int k=0;k<MaxLO;k++){
				Elem tmp(best);
				Float_t faktor=0.0625*(1-(Float_t)k/(Float_t)MaxLO);
				_rnd(tmp,best,hbPosRange*faktor,faktor*hbDreRange,faktor*hbBondRange);
				_eval(tmp);
				if(tmp.score<best.score)
					best=tmp;
			}
			if(best.score<ALLbest.score){
				//printf("%5u %3u %+12f => %+12f \n",counts,cc,nest.score,best.score);
				printf("%i new global best %lf\n",evalCalls,best.score);
				cc=0;
				ALLbest=best;
			}
			nest=best;
		}else if (cc>20){
		//wenn mehr als 20 Laeufe keine besserung->neue Posi
			_rnd(nest,_center,_radius);
			_eval(nest);
			best=nest;
			cc=0;
		}else{
			size_t  scoutidx(0);
			size_t  scoutcount(0);
			Float_t scoutval(0);
			for(size_t idx=0;idx!=scout_size;idx++){
				if(scoutcho[idx]>scoutcount || (scoutcho[idx]==scoutcount && scout[idx].score<scoutval)){
				//mehr Follower oder gleiche + kleinerer wert;
					scoutcount=scoutcho[idx];
					scoutidx=idx;
					scoutval=scout[idx].score;
				}
			}
			if (FollowCount>0 && scoutcount>=FollowCount){
				nest=scout[scoutidx];
			}else
			{
			//derzeit beste moeglichkeit absuchen
				_rnd(nest,ALLbest,hbfPosRange,hbfDreRange,hbfBondRange);
				_eval(nest);
				best=nest;

			}
		}
	}
	//local optimization
	unsigned int MaxLO=4096;
	for(unsigned int k=0;k<MaxLO;k++){
		Elem tmp(ALLbest);
		Float_t faktor=0.125*(1-(Float_t)k/(Float_t)MaxLO);
		_rnd(tmp,ALLbest,hbPosRange*faktor,faktor*hbDreRange,faktor*hbBondRange);
		_eval(tmp);
		if(tmp.score<ALLbest.score){
			printf("%i new global best %lf\n",evalCalls,tmp.score);
			ALLbest=tmp;
		}
	}
	_best=ALLbest;
	// output best
	if (outend==1)
	{
		outfile<<"#structure fitness =\""<<_best.score<<"\""<<endl;
		_fit->writeOptimum(_best.ort, _best.dreh, _best.bond, outfile);
	}
	//output all
	if (outend==2)
	{
		for (unsigned int i =0;i<scout_size;i++){
			Elem * j(&scout[i]);
			outfile<<"#structure fitness =\""<<j->score<<"\""<<endl;
			_fit->writeOptimum(j->ort, j->dreh, j->bond, outfile);
			
		}
		for (unsigned int i =0;i<follo_size;i++){
			Elem * j(&follower[i]);
			outfile<<"#structure fitness =\""<<j->score<<"\""<<endl;
			_fit->writeOptimum(j->ort, j->dreh, j->bond, outfile);
			
		}
		
	}
	delete []scout;
	delete []scoutcho;
	delete []follower;
	delete []Summi;
	return ALLbest.score;
}
Float_t HoneyBee::_runs(ostream& outfile, unsigned int outiter, unsigned int outend,unsigned int runs)
{
    if(!runs)return 0;
    Float_t erg(0);
#pragma omp parallel for
    for( int k=0;k<runs;k++){
	Float_t tmp(_run(outfile,outiter,outend));
#pragma omp atomic
	erg+=tmp;
    }
    erg/=runs;
    return erg;
}
void HoneyBee::run(ostream& outfile, unsigned int outiter, unsigned int outend)
{
	fprintf(stdout,"SingleCalc @%6u:%3u\n",_iter,_particle_count);
	Float_t bk=0.003616  ; Float_t bl=0.001084  ; Float_t bm=0.027854  ;
	Float_t bn=0.025366  ; Float_t bo=0.039257  ; Float_t bp=0.012289  ;
	evalCalls=0;
		hbPosRange   = bk * _radius;
		hbDreRange   = bl * _2PI   ;
		hbBondRange  = bm * _2PI   ;
		hbfPosRange  = bn * _radius;
		hbfDreRange  = bo * _2PI   ;
		hbfBondRange = bp * _2PI   ;
		Float_t erg=_run(outfile,outiter,outend);
		printf("score: %f\n",erg);
	return;

}
HoneyBee::~HoneyBee() {
	// TODO Auto-generated destructor stub
}
RandomLocalSearch::RandomLocalSearch(vector<string>& parameter, paradocks::Fitness* fit, RNG& rng, Vec3_t& center, Float_t& radius):
HoneyBee(parameter,fit,rng,center,radius){

}

void RandomLocalSearch::run(ostream& outfile, unsigned int outiter, unsigned int outend){

	Float_t bk=0.003616  ; Float_t bl=0.001084  ; Float_t bm=0.027854  ;
	Float_t bn=0.025366  ; Float_t bo=0.039257  ; Float_t bp=0.012289  ;

		hbPosRange   = bk * _radius;
		hbDreRange   = bl * _2PI   ;
		hbBondRange  = bm * _2PI   ;
		hbfPosRange  = bn * _radius;
		hbfDreRange  = bo * _2PI   ;
		hbfBondRange = bp * _2PI   ;
unsigned int cc=0;
	//Arrays für die Scouts und Follower
	Elem   *follower = new Elem  [_particle_count];
	//Nest Standort festlegen
	Elem nest;
	_rnd(nest,_center,_radius);
	_eval(nest);
	//Best und AllBest festlegen
	Elem best(nest);
	Elem ALLbest(nest);
	unsigned int counts=0;
//	int printi=0;
	Float_t faktor;
	while(counts<_iter){
		counts++;
		faktor=1.0-counts/_iter;
#pragma omp parallel for
		for(unsigned int k=0;k<_particle_count;k++){

			_rnd(follower[k],nest,faktor*hbfPosRange,faktor*hbfDreRange,faktor*hbfBondRange);
			_eval(follower[k]);

#pragma omp critical
			if(follower[k].score <best.score)
				best=follower[k];
#pragma omp end critical

		}
		cc++;

		if(best.score<nest.score){
			//local optimization
			unsigned int MaxLO=4;
			for(unsigned int k=0;k<MaxLO;k++){
				Elem tmp(best);
				Float_t faktor=0.0625*(1-(Float_t)k/(Float_t)MaxLO);
				_rnd(tmp,best,hbPosRange*faktor,faktor*hbDreRange,faktor*hbBondRange);
				_eval(tmp);
				if(tmp.score<best.score)
					best=tmp;
			}
			if(best.score<ALLbest.score){
				//printf("%5u %3u %+12f => %+12f \n",counts,cc,nest.score,best.score);
				cc=0;
				ALLbest=best;
			}
			nest=best;
		}else if (cc>20){
		//wenn mehr als 20 Laeufe keine besserung->neue Posi
			_rnd(nest,_center,_radius);
			_eval(nest);
			best=nest;
			cc=0;
		}else{
				_rnd(nest,ALLbest,hbfPosRange,hbfDreRange,hbfBondRange);
				_eval(nest);
				best=nest;
		}
	}
	//local optimization
	unsigned int MaxLO=4096;
	for(unsigned int k=0;k<MaxLO;k++){
		Elem tmp(ALLbest);
		Float_t faktor=0.125*(1-(Float_t)k/(Float_t)MaxLO);
		_rnd(tmp,ALLbest,hbPosRange*faktor,faktor*hbDreRange,faktor*hbBondRange);
		_eval(tmp);
		if(tmp.score<ALLbest.score)
			ALLbest=tmp;
	}
	delete []follower;
	_best=ALLbest;
	return;
}
}
