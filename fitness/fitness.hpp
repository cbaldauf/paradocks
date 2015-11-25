#ifndef FFITNESS_H
#define FFITNESS_H

#include <sstream>
#include "mgraph/graph.hpp"
#include <boost/function.hpp>

namespace paradocks
{
using namespace std;
using namespace mgraph;

/**
\defgroup fitness ParaDockS Fitness function Interface 
**/

/**
\brief Fitness error class.
This class should be thrown by all implemented Fitness classes in case of an error.
Derived from std::runtime_error.
\ingroup fitness
\ingroup framework
**/
class FitnessError : public runtime_error
{
public:
	FitnessError(string err) : runtime_error("FitnessError: "+err) {}
};

/**
\brief Fitness class.
Virtual base class interface definition for all Fitness class implementations.
\ingroup fitness
\ingroup framework
**/
class Fitness
{
public:
	/// Constructor
	Fitness(vector<string>& parameter) : _parameter(parameter) {}
	/// Destructor
	virtual ~Fitness() {}
	
	/**
	\brief Initialise the fitness function with the protein.
	\param prot The protein as ChainList.
	\param as Center of the active site.
	\param r Radius of the active site.
	**/
	virtual void initProtein(const ChainList& prot, const Vec3_t& as, Float_t r);
	
	/**
	\brief Return the number of terms implemented in the fitness function.
	Pure virtual function!
	\return Number of fitness (scoring) terms returned by eval.
	**/
	virtual unsigned int TermsCount() = 0;
	
	/**
	\brief Initialise the Fitness function with the ligand.
	Pure virtual function!
	\param lig The ligand as ChainList.
	\return Number of rotatable bonds in the ligand.
	**/
	virtual void initLigand(const ChainList& lig) = 0;
	
	/**
	\brief Return the number rotatable bonds.
	\return Number of rotatable bonds in the current ligand.
	**/
	virtual unsigned int RotCount() { return _rotcount; }
	
	/**
	\brief Evaluate the current configuration.
	This method evaluates the configuration of the ligand specified in the parameter
	and returns a field of scoring values. Pure virtual function!
	\param pos Position of the centerfragment.
	\param ori Orientation of the centerfragment.
	\param angle Vector of the ligand's torsion angles. The size needs to be equal to TermsCount().
	\param result A pointer to the vector with the results.
	**/
	virtual void eval(const Vec3_t& pos, const Quat_t& ori,
		const vector<Float_t>& angle, vector<Float_t>*& result) = 0;
	
	/**
	\brief Write out a ligands configuration as Tripos MOL2 data.
	This method writes out a ligands configuration as Tripos MOL2 data to a stream.
	\param pos Position of the centerfragment.
	\param ori Orientation of the centerfragment.
	\param angle Vector of the ligand's torsion angles. The size needs to be equal to TermsCount().
	\param outfile The Tripos MOL2 data will be writen to this stream.
	**/
    void writeOptimum(const Vec3_t& pos, const Quat_t& ori,
    	const vector<Float_t>& angle, ostream& outfile);
	
	/// write initial coordinates back
	void restoreLigand();
	/// return ligand offset
	Vec3_t Offset() { return _offset; }
	/// All output from the fitness class should go into this stringstream.
	stringstream log;
	
protected:
	/**
	\brief Provides some code for initLigand(const ChainList& lig).
	This method finds rotatable bonds in the ligand and does some preparation.
	\param lig The ligand as ChainList.
	\param h_handling Handling of hydrogens:\n
			0 - no hydrogens are rotated \n
			1 - polar hydrogens are rotated \n
			2 - all hydrogens are rotated \n
	**/
	void initLigand(const ChainList& lig, unsigned int h_handling);
	
	// parameter for the fitness function
	vector<string>& _parameter;
	
	// protein
	// center of the active site
	Vec3_t _center;
	// radius of the active site
	Float_t _radius;
	// protein atoms
	AtomList _patoms;
	
	// ligand
	// ligand atoms
	AtomList _latoms;
	// ligand offset, for reconstruction of original coordinates
	Vec3_t _offset;
	// initial ligand coordinates
	vector<Vec3_t> _lcoord;
	// rigid fragments are saved in this struct
	struct rotfragment
	{
		unsigned int from, to;
		int dist;
		vector<unsigned int> atomindex;
	};
	struct less_rotfragment : public binary_function<rotfragment, rotfragment, bool> {
	bool operator()(rotfragment x, rotfragment y) { return x.dist > y.dist; }
	};
	// by rotatable bonds connected fragments
	vector<rotfragment> _allrotfragments;
	// the number of rotatable bonds
	unsigned int _rotcount;
	
	// extends to rigid fragment
	AtomList RigidFragment(AtomKey a, unsigned int h_handling);
	// returns true if this atom ist the last heavy atom of a branch
	bool BranchEnd(AtomKey a, unsigned int h_handling);
	// returns every atom af the branch
	AtomList CollectBranch(AtomKey from, AtomKey to);
	// update the ligand coordinats
	vector<Vec3_t> updateLigandCoord(const Vec3_t& pos, const Quat_t& ori,
		const vector<Float_t>& angle);
	
};

/**
\brief Grid class.
This class can be used to aproximate the energie value within a grid by linear triangulation.
\ingroup fitness
\ingroup framework
**/
class Grid
{
public:
	/**
	\brief Construct an empty grid.
	**/
	Grid() {}
	virtual ~Grid() {}
	/**
	\brief setup the position and size of the grid.
	\param p1 point defining one corner of the grid box.
	\param p2 point defining another corner of the grid box.
	\param spacing lattice spacing.
	\param typecount number of grids.
	\param outside function value for all points outside the grid.
	**/
	void setupGrid(const Vec3_t& p1, const Vec3_t& p2, Float_t spacing, unsigned int typecount, Float_t outside);
	/**
	\brief fill a grid with precalculated values.
	\param type the type which should be calculated.
	\param calc callback to the evaluation function.
	**/
	void calculate(unsigned int type, boost::function<Float_t(Vec3_t)>& calc);
	/**
	\brief check if there are allready values in the grid for this type.
	\param type type.
	\return true if the grid of this type is not empty.
	**/
	bool iscalculated(unsigned int type);
	/**
	\brief trilinear interpolation.
	\param p position.
	\param type type.
	\return interpolated value.
	**/
	Float_t interpolate(const Vec3_t& p, unsigned int type);
protected:
	// offset corner off the grid
	Vec3_t _offset;
	// other end of the grid
	Vec3_t _end;
	// lattice spacing
	Float_t _spacing;
	Float_t _outside;
	unsigned int _xsize, _ysize, _zsize;
	// the grid
	vector<vector<vector<vector<Float_t> > > > _grid;
};

/**
\brief StaticNeighborhood class.
This class provides a static neighborhood list for atoms.
\ingroup fitness
\ingroup framework
**/
class StaticNeighborhood
{
public:
	/**
	\brief Construct an empty neighborhood list.
	**/
	StaticNeighborhood() {}
	virtual ~StaticNeighborhood() {}
	/**
	\brief Initialize the neighborhood list.
	\param al atoms, which should be put in the neighborhood list.
	\param cutoff maximum distance for atoms which belong to the neighborhood.
	\param spacing lattice spacing.
	**/
	void setupList(const AtomList& al, Float_t cutoff, Float_t spacing);
	/**
	\brief Return the indices of the atoms in the neighborhood.
	\param p the returned neighborhood list is arround position p.
	**/
	const vector<unsigned int> getNeighborhood(const Vec3_t& p);
protected:
	// offset corner off the neighborhood lattice
	Vec3_t _offset;
	// end of the lattice
	Vec3_t _end;
	// lattice spacing
	Float_t _spacing;
	unsigned int _xsize, _ysize, _zsize;
	// the neighborhood list
	vector<vector<vector<vector<unsigned int> > > > _neighborhoodlist;
	vector<unsigned int> _empty;
};

}

#endif /*FFITNESS_H*/
