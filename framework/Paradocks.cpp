#include "Paradocks.hpp"

#include <sys/time.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <fstream>

// all available fitness functions
#include "fitness/pscore/pscore.hpp"
#include "fitness/pmf04/pmf04.hpp"

// all available optimizer
#include "optimizer/pso/pso.hpp"
#include "optimizer/random/random.hpp"
#include "optimizer/random_hill_climbing/random_hill_climbing.hpp"
#include "optimizer/hbee/HoneyBee.hpp"

using namespace mgraph;

namespace paradocks
{

Paradocks::Paradocks(ParadocksConfig* cfg) : _cfg(cfg)
{
	// prepare protein input information
	try
	{
		_as=Vec3_t(boost::lexical_cast<Float_t>(cfg->X()),
					boost::lexical_cast<Float_t>(cfg->Y()),
					boost::lexical_cast<Float_t>(cfg->Z()));
		_rad=boost::lexical_cast<Float_t>(cfg->Rad());
	}
	catch(boost::bad_lexical_cast e)
	{
		string emsg("Error in input conversion! <protein> line\n"+string(e.what()));
		throw ParadocksError(emsg);
	}
	
	// prepare ligand input information
	try
	{
		for (vector<string>::iterator i=_cfg->LigandIdx().begin(); i!=_cfg->LigandIdx().end(); ++i)
		{
			_ligandIdx.push_back(boost::lexical_cast<unsigned int>(*i));
		}
	}
	catch(boost::bad_lexical_cast e)
	{
		string emsg("Error in input conversion! <ligand> line\n"+string(e.what()));
		throw ParadocksError(emsg);
	}
	
	// prepare number of runs
	try
	{
		_runs=boost::lexical_cast<unsigned int>(cfg->Runs());
	}
	catch(boost::bad_lexical_cast e)
	{
		string emsg("Error in input conversion! <runs> line\n"+string(e.what()));
		throw ParadocksError(emsg);
	}
	
	// prepare seed
	try
	{
		unsigned int seed=boost::lexical_cast<unsigned int>(cfg->Seed());
		if (!seed)
		{
			timeval tv;
			gettimeofday(&tv, NULL);
			seed=(unsigned int)tv.tv_usec;
			cout<<"Seed random from time with "<<seed<<"."<<endl;
		}
		_rng.seed(RNG::result_type(seed));
	}
	catch(boost::bad_lexical_cast e)
	{
		string emsg("Error in input conversion! <seed> line\n"+string(e.what()));
		throw ParadocksError(emsg);
	}
	
	// prepare fitness function
	if (cfg->FitnessType()=="pscore") _fit.reset(new PScore(cfg->FitnessPar()));
	if (cfg->FitnessType()=="pmf04") _fit.reset(new PMF04(cfg->FitnessPar()));
	if (!_fit)
	{
		string emsg("No valid fitness function selected!");
		throw ParadocksError(emsg);
	}
	
	// prepare optimizer type
	if (cfg->OptimizerType()=="pso") _opt.reset(new PSO(cfg->OptimizerPar(), _fit.get(), _rng, _as, _rad));
	if (cfg->OptimizerType()=="random") _opt.reset(new RANDOM(cfg->OptimizerPar(), _fit.get(), _rng, _as, _rad));                                                                   
    if (cfg->OptimizerType()=="rhc") _opt.reset(new RANDOM_HILL_CLIMBING(cfg->OptimizerPar(), _fit.get(), _rng, _as, _rad));   
	if (cfg->OptimizerType()=="hbee") _opt.reset(new HoneyBee(cfg->OptimizerPar(), _fit.get(), _rng, _as, _rad));
	if (!_opt)
	{
		string emsg("No valid optimizer selected!");
		throw ParadocksError(emsg);
	}
	
	// prepare output
	if (cfg->OutIteration()=="none") _outiter=0;
	else if (cfg->OutIteration()=="best") _outiter=1;
	else if (cfg->OutIteration()=="all") _outiter=2;
	else
	{
		string emsg(cfg->OutIteration()+" is not a valid type for <iteration> line");
		throw ParadocksError(emsg);
	}
	if (cfg->OutEnd()=="best") _outend=1;
	else if (cfg->OutEnd()=="all") _outend=2;
	else
	{
		string emsg(cfg->OutEnd()+" is not a valid type for <end> line");
		throw ParadocksError(emsg);
	}
	_outprefix=cfg->OutPrefix();	
}

void Paradocks::run()
{
	// read protein and broadcast data
	ChainList prot;
	prot=ReadMOL2(_cfg->Protein());
	
	// initialize fitnessfunction with the protein
	_fit->initProtein(prot, _as, _rad);
	
	// loop over the ligands
	// errors within this loop should only affect the current ligand
	for (unsigned int i=0; i!=_cfg->Ligand().size(); ++i)
	{
		try
		{
			// read ligand and broadcast data
			ChainList lig;
			lig=ReadMOL2(_cfg->Ligand()[i]);
			
			// get the indexed ligand
			ChainList selectedligand;
			if (_ligandIdx[i]>lig.size())
			{
				string emsg("No ligand with index " + boost::lexical_cast<string>(_ligandIdx[i]) + " in file \"" + _cfg->Ligand()[i] + "\".");
				throw ParadocksError(emsg);
			}
			
			selectedligand.push_back(lig[_ligandIdx[i]-1]);
			
			// initialize fitnessfunction with the ligand
			_fit->initLigand(selectedligand);
			
			for (unsigned int j=0; j<_runs; ++j)
			{
				boost::filesystem::path protpath(_cfg->Protein());
				boost::filesystem::path ligpath(_cfg->Ligand()[i]);
				boost::filesystem::path filename(_outprefix+string("_")+
					boost::filesystem::basename(protpath)+
					string("_")+
					boost::filesystem::basename(ligpath)+
					string("_")+
					boost::lexical_cast<string>(_ligandIdx[i])+
					string("_")+
					boost::lexical_cast<string>(j+1)+string(".mol2"));
				
				cout<<"Open outputfile \""<<filename.string()<<"\"."<<endl;
				boost::filesystem::ofstream outfile(filename);
				if (!outfile.good())
				{
					string emsg("Can't open outputfile \""+filename.string()+"\".");
					throw ParadocksConfigError(emsg);
				}
					
				outfile<<"#protein file = \""<<_cfg->Protein()<<"\""<<endl;
				outfile<<"#ligand file = \""<<_cfg->Ligand()[i]<<"\""<<endl;
				_opt->run(outfile, _outiter, _outend);
				outfile.close();
			}
			
			clDelete(lig);
		}
		// catch the error and proceed with the next ligand
		catch(std::runtime_error e)
		{
			cout<<e.what()<<endl;
			cout<<"Proceed with the next ligand!"<<endl;
		}
	}
}

}
