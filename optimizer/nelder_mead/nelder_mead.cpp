#include "nelder_mead.hpp"


namespace paradocks
{

NELDER_MEAD::Parameter NELDER_MEAD::move(NELDER_MEAD::Parameter par)
{
	boost::variate_generator<RNG&, boost::uniform_real<Float_t> > r(_rng, boost::uniform_real<Float_t>(-1,1));
	Parameter result;
	
	Float_t factor=4;
	
	Float_t maxpos=0.1*factor;
	Float_t maxang=osg::PI/45*factor;
	do
	{
		result._pos=Vec3_t(r(), r(), r())*maxpos;
	} while (result._pos.length()>maxpos);
	result._pos=par._pos+result._pos;
	
	boost::variate_generator<RNG&, boost::uniform_real<Float_t> > rsimple(_rng, boost::uniform_real<Float_t>(0,1));
	Float_t x0(rsimple());
	Float_t t1(2*osg::PI*rsimple()), t2(2*osg::PI*rsimple());
	Float_t s1(sin(t1)), s2(sin(t2)), c1(cos(t1)), c2(cos(t2));
	Float_t r1(sqrt(1-x0)), r2(sqrt(x0));
	Quat_t randomori=Quat_t(s1*r1, c1*r1, s2*r2, c2*r2);
	result._ori.slerp(maxang, Quat_t(), randomori);
	Quat_t temp=par._ori*result._ori;
	result._ori=temp;
	
	// bond rotatio
	result._bond=vector<Float_t>(_fit->RotCount());
	for (unsigned int i=0; i!=_fit->RotCount(); ++i)
	{
		result._bond[i]=par._bond[i]+r()*maxang;
	}
	return result;
}


NELDER_MEAD::NELDER_MEAD(vector<string>& parameter, paradocks::Fitness* fit, RNG& rng, Vec3_t& center, Float_t& radius)
	: Optimizer(parameter, fit, rng, center, radius)
{
	if (parameter.size()!=2)
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
	try
	{
		_local_iter=boost::lexical_cast<unsigned int>(parameter[1]);
	}
	catch(boost::bad_lexical_cast e)
	{
		string emsg("Error in input conversion! local optimizer iteration line\n"+string(e.what()));
		throw OptimizerError(emsg);
	}

}

NELDER_MEAD::Parameter NELDER_MEAD::init_new_Parameter(Vec3_t offset)
{
	boost::variate_generator<RNG&, boost::uniform_real<Float_t> > r(_rng, boost::uniform_real<Float_t>(-1,1));
	Parameter result;
	// position
	do
	{
		result._pos=offset+Vec3_t(r() ,r(), r());
	} while ((offset-result._pos).length()>1);
	
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
	Quat_t randomori=Quat_t(s1*r1, c1*r1, s2*r2, c2*r2);
	result._ori.slerp(osg::PI/6, Quat_t(), randomori); //30 degree max 
	
	// bond rotation
	result._bond=vector<Float_t>(_fit->RotCount());
	for (unsigned int i=0; i!=_fit->RotCount(); ++i)
	{
		result._bond[i]=r()*osg::PI/6;
	}
	return result;
	vector<Float_t>* score;
	_fit->eval(result._pos, result._ori, result._bond, score);
	result._score=accumulate(score->begin(), score->end(), (Float_t)0);
}


struct less_parameter : public binary_function<NELDER_MEAD::Parameter,NELDER_MEAD::Parameter, bool> {
bool operator()(NELDER_MEAD::Parameter x, NELDER_MEAD::Parameter y) { return x._score < y._score; }
};

void NELDER_MEAD::run(ostream& outfile, unsigned int outiter, unsigned int outend)
{
	Vec3_t offset=_fit->Offset();
	int simplexsize=6+_fit->RotCount();
	list<Parameter> simplex;
	for (unsigned int i=0; i!=simplexsize; ++i)
	{
		simplex.push_back(init_new_Parameter(offset));
	}
	
	simplex.sort(less_parameter());
	
	//less_parameter grr;
	
	//sort(simplex.begin(), simplex.end());
	
	_best=init_new_Parameter(offset);
	vector<Float_t>* score;
	_fit->eval(_best._pos, _best._ori, _best._bond, score);
	_best._score=accumulate(score->begin(), score->end(), (Float_t)0);

	for (unsigned int i=0; i!=_iter; ++i)
	{
		Parameter np=init_new_Parameter(offset);
		_fit->eval(np._pos, np._ori, np._bond, score);
		np._score=accumulate(score->begin(), score->end(), (Float_t)0);
		// refine by local minimization

		for (unsigned int j=0; j!=_local_iter; ++j)
		{
			Parameter x=move(np);
			_fit->eval(x._pos, x._ori, x._bond, score);
			x._score=accumulate(score->begin(), score->end(), (Float_t)0);
			if (x._score<np._score)
			{
				np=x;
			}
		}

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








/*
#include "Main.h"

#include "MinFinder.h"
#include "Tools.h"
#include "Point.h"
#include "NelderMead.h"*/


/*NelderMead::NelderMead(Config *c) : MinFinder(c)
{
}


NelderMead::~NelderMead()
{
}

Point *NelderMead::findMinima(Point *startPoint)
{
    vector<Point *> simplex;
    int runs, reduceCount;
    Point *point, *star1, *star2, *centroid;
    int highID, secondHighID, lowID;
    double el, es, eh; //low, secondHigh, high energies
    double startError; 


    highID       = cfg->varCount;
    secondHighID = cfg->varCount-1;
    lowID        = 0;

    //create n+1 simplex points:

    for(int i = 0; i <= cfg->varCount; i++) {
        point = new Point(cfg->varCount);

        for(int j = 0; j < cfg->varCount; j++) {
            if(i == 0) {
                if(startPoint == NULL)
                    (*point)[j] = Tools::randomDouble(cfg->lowerBounds[j]+cfg->stepSizes[j],cfg->upperBounds[j]-cfg->stepSizes[j]);
                else
                    (*point)[j] = (*startPoint)[j];
            }
            else {
                int stepsToHigh = (cfg->upperBounds[j] - simplex[0]->at(j)) / cfg->stepSizes[j];
                int stepsToLow = (simplex[0]->at(j) - cfg->lowerBounds[j]) / cfg->stepSizes[j];

                if(stepsToHigh >= stepsToLow)
                    (*point)[j] = Tools::randomDouble(simplex[0]->at(j), cfg->upperBounds[j]);
                else
                    (*point)[j] = Tools::randomDouble(simplex[0]->at(j), cfg->lowerBounds[j]);
            }   
        }

//        cfg->correctPoint(*point);
        point->setEnergy(cfg->evalBorderCheck(*point));
        //point->setDelta(deltaEval(*point));
        simplex.push_back(point);
    }



    //order by energy:
//    sort(simplex.begin(), simplex.end(), Point::NMCmp);

    sort(simplex.begin(), simplex.end(), Config::PointPtrLess());
    for(int i = 1; i < (signed int) simplex.size(); i++) {
        assert(*simplex[0] <= *simplex[i]);
    }
    centroid = calcCentroid(simplex); 
    runs = 0;
    reduceCount = 0;
    startError = stdError(simplex);
    while(stdError(simplex) > cfg->NMStdError && runs < cfg->NMMaxIterations ) { 
        runs += 1;
        //store values:
        el = simplex[lowID]->getEnergy(); 
        es = simplex[secondHighID]->getEnergy();
        eh = simplex[highID]->getEnergy();

        //Reflection
        star1 = reflect(*simplex[highID], *centroid);

        if(star1->getEnergy() < es && star1->getEnergy() >= el) { 

            //position new point:
            updateCentroid(*centroid, *simplex[secondHighID], *star1);
            int i = cfg->varCount-1; 
            delete simplex[cfg->varCount];
            while(i > 0 && *simplex[i] > *star1) { 
                simplex[i+1] = simplex[i];
                i--;
            }
            simplex[i+1] = star1;
        }
        else if(star1->getEnergy() < el) {//Expand
            star2 = expand(*star1, *centroid);

            if(star2->getEnergy() < star1->getEnergy()) {  //expand success
                delete star1;
                star1 = star2;
            }
            else
                delete star2;
            updateCentroid(*centroid, *simplex[secondHighID], *star1);

            //position star1 in front:
            delete simplex[cfg->varCount];
            for(int i = cfg->varCount - 1; i >= 0; i--)
                simplex[i+1] = simplex[i];
            simplex[0] = star1;

        }
        else { //eStar1 >= e[i] for all i != highID 
            if(star1->getEnergy() < eh) {
                delete simplex[highID];
                simplex[highID] = star1;
                eh = star1->getEnergy();
            }
            else 
                delete star1;
            //Contraction
            star2 = contract(*simplex[highID], *centroid);

            if(star2->getEnergy() < eh) { //cotraction success

                if(star2->getEnergy() < es) //FIXME order
                    updateCentroid(*centroid, *simplex[secondHighID], *star2);

                int i = cfg->varCount - 1;
                delete simplex[cfg->varCount];
                while(i >= 0 && simplex[i]->getEnergy() > star2->getEnergy()) { 
                    simplex[i+1] = simplex[i];
                    i--;
                }
                simplex[i+1] = star2;
            }
            else { //contraction failed
                delete star2;
                reduceCount++;
                reduce(simplex);
                sort(simplex.begin(), simplex.end(), Point::energyCmp); 
                delete centroid;
                centroid = calcCentroid(simplex); 
            }
        }
        //delete centroid;
        //centroid = calcCentroid(simplex); 
    } //while stderror 
    if(runs >= cfg->NMMaxIterations) {
        cout << "Warning: [Retry] To many iterations on NM: " << cfg->NMMaxIterations << endl;
        return findMinima();
    }

    simplex[0]->setReliability((double)(runs-reduceCount)*100/(double)runs);
    simplex[0]->minInitStdError = startError;
    simplex[0]->minRuns = runs;

    //clean up
    delete centroid;
    for(int i = 1; i <= cfg->varCount; i++) {
        delete simplex[i];
    }

        //cout << *simplex[0] << " " << simplex[0]->getEnergy() << endl << flush;
        return simplex[0];

}



/* Calculates the centroid of the n-1 best points
 * note: cfg->varcount + 1 = n
 */
/*Point *NelderMead::calcCentroid(vector<Point *> &simplex)
{
    Point *pos;

    pos = new Point(cfg->varCount);
    for(int i = 0; i < cfg->varCount; i++) {
        Point *tmp = simplex[i];
        for(int j = 0; j < cfg->varCount; j++) {
            (*pos)[j] += tmp->at(j);
        }
    }
    for(int j = 0; j < cfg->varCount; j++)
        (*pos)[j] /= cfg->varCount; 

    //cfg->correctPoint(*pos);
    pos->setEnergy(cfg->evalBorderCheck(*pos));
//    pos->setDelta(deltaEval(*pos));
    return pos;
}

void NelderMead::updateCentroid(Point &centroid, Point &remove, Point &insert)
{

    for(int i = 0; i < cfg->varCount; i++) {
        centroid[i] += (insert[i] - remove[i])/cfg->varCount;
    }
    //cfg->correctPoint(centroid);
    centroid.setEnergy(cfg->evalBorderCheck(centroid));
}






Point *NelderMead::reflect(Point &p, Point &cen)
{
    double alpha = cfg->NMReflection;
    Point *pos = new Point(0);

    for(int i = 0; i < cfg->varCount; i++)
        pos->push_back(cen[i] + alpha*(cen[i]-p[i]));

    //cfg->correctPoint(*pos);

    pos->setEnergy(cfg->evalBorderCheck(*pos));
//    pos->setDelta(deltaEval(*pos));

    return pos;
}


Point *NelderMead::expand(Point &p, Point &cen)
{
    double gamma = cfg->NMExpansion;
    Point *pos;

    pos = new Point(cfg->varCount);

    for(int i = 0; i < cfg->varCount; i++)
        (*pos)[i] = cen[i] + gamma*(p[i] - cen[i]);
    //cfg->correctPoint(*pos);

    pos->setEnergy(cfg->evalBorderCheck(*pos));
//    pos->setDelta(deltaEval(*pos));

    return pos;
}


Point *NelderMead::contract(Point& p, Point& cen)
{
    double beta = cfg->NMContraction;
    Point *pos;

    pos = new Point(cfg->varCount);

    for(int i = 0; i < cfg->varCount; i++)
        (*pos)[i] = cen[i] + beta*(p[i]-cen[i]);

    //cfg->correctPoint(*pos);

    pos->setEnergy(cfg->evalBorderCheck(*pos));
//    pos->setDelta(deltaEval(*pos));

    return pos;
}

void NelderMead::reduce(vector<Point *> &simplex)
{

    double lambda = cfg->NMReduction;

    for(int i = 1; i < cfg->varCount + 1; i++) {
        for(int j = 0; j < cfg->varCount; j++) {
            (*simplex[i])[j] = simplex[0]->at(j) + lambda*(simplex[i]->at(j) - simplex[0]->at(j)); 
        }

        //cfg->correctPoint(*simplex[i]);

        simplex[i]->setEnergy(cfg->evalBorderCheck(*simplex[i]));
    }

}



double NelderMead::stdError(vector<Point *> &simplex)
{
    double err;
    Point *cen = calcCentroid(simplex); //FIXME needed recalced? 

    err = 0;
    for(int i = 0; i < cfg->varCount + 1; i++)
        err    += pow(simplex[i]->getEnergy() - cen->getEnergy(), 2);  
    err /= (cfg->varCount + 1);
    delete cen;
    return err;
}
*/
