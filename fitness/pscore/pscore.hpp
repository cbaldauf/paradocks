#ifndef PSCORE_H
#define PSCORE_H

#include "fitness/fitness.hpp"
#include "mgraph/graph.hpp"

#define VDWTYPECOUNT 45

namespace paradocks
{

class PScore : public Fitness
{
public:
	/// Constructor
	PScore(vector<string>& parameter);
	/// Destructor
	~PScore() {};
	
	/// Initialise the Fitness function with the protein.
	void initProtein(const ChainList& prot, const Vec3_t& as, Float_t r);
	
	/// Return the number of terms implemented in the fitness function.
	unsigned int TermsCount() { return _termscount; }
	
	/// Initialise the Fitness function with the ligand.
	void initLigand(const ChainList& lig);
	
	///calculate the vdw interaction for a single type at a position pos
	Float_t calcvdw(const Vec3_t& pos, unsigned int type);
	
	/// Evaluate the current configuration.
	void eval(const Vec3_t& pos, const Quat_t& ori,
		const vector<Float_t>& angle, vector<Float_t>*& result);
	
	/// vdw radius
	static Float_t vdW_type[VDWTYPECOUNT];
	
private:
	// general fitness function informations
	// the vector with the scoring terms (results)
	enum {_termscount=2};
	vector<Float_t> _result;
	
	struct DA
	{
		// can be frozen(0), rot(1) or free(2)
		unsigned int state;
		unsigned int root; // index of root atom
		unsigned int d;
		vector<unsigned int> h; // indices of hydrogens
		unsigned int a;
		vector<unsigned int> lp; // indices of lone pairs
	};
	
	struct HBond
	{
		DA donor;
		DA acceptor;
		Float_t f1;
		Float_t f2;
		Float_t f3;
		Float_t f;
	};
	
	struct less_hbond : public binary_function<HBond, HBond, bool> {
	bool operator()(HBond x, HBond y) { return x.f < y.f; }
	};
	struct more_hbond : public binary_function<HBond, HBond, bool> {
	bool operator()(HBond x, HBond y) { return x.f > y.f; }
	};

	// protein informations
	// vdw protein atom neighbourhood list
	StaticNeighborhood _vdw_nh;
	// protein atoms for vdW interaction
	vector<unsigned int> _vdw_patom;
	// protein vdW parameter
	vector<Float_t> _pvdw;
	// vdw grid
	Grid _vdw_grid;
		
	// donor-acceptor atom neighbourhood list
	StaticNeighborhood _da_nh;
	// protein atoms for h-bonding
	vector<unsigned int> _da_patom;
	// protein donor-acceptor radius (same like vdw)
	vector<Float_t> _da_pradius;
	// protein donor-acceptor parameter
	vector<DA> _pda;

	
	// ligand informations
	// ligand atoms for vdW interaction
	vector<unsigned int> _vdw_latom;
	// ligand vdW parameter
	vector<unsigned int> _lvdw;
	// internal ligand vdw clashes penalty
	vector<vector<unsigned int> > _internal_vdw_neighborhood_list;
	vector<vector<Float_t> > _internal_vdw_neighborhood_list_d0;
	
	// ligand atoms for h-bonding
	vector<unsigned int> _da_latom;
	// ligand donor-acceptor radius (same like vdw)
	vector<Float_t> _da_lradius;
	// da ligand
	vector<DA> _lda;

	// return the vdw type of the atom
	unsigned int vdWType(AtomKey a);
	// return the h-bonding type of the atom, need the right al for indices
	bool daType(AtomKey a, DA& type, AtomList& al);
	

	// returns all atoms which are more than n-1 bonds away in the same subgraph
	AtomList One_N_Neighborhood(AtomKey a, unsigned int n);
	// add lp for h-bonding
	void AddLP(const ChainList& cl);
	
	Float_t _hbfaktor;
	Float_t _hbdist;
	Float_t _hbangle;
	Float_t _inv_hbdist;
	Float_t _inv_hbangle;
};

}
#endif /*XSCORE_H*/


