#ifndef PMF04_H
#define PMF04_H

#include "fitness/fitness.hpp"
#include "mgraph/graph.hpp"



namespace paradocks
{

class PMF04 : public Fitness
{
public:
	/// Constructor
 	PMF04(vector<string>& parameter);
/*	PMF04(vector<string>& parameter) : Fitness(parameter), _result(_termscount)
	{
		PARAMETER=parameter[0];
		_SF1=boost::lexical_cast<double>(parameter[1]);
		_SF2=boost::lexical_cast<double>(parameter[2]);
		_SF3=boost::lexical_cast<double>(parameter[3]);
		_vdw_range=boost::lexical_cast<double>(parameter[4]);
		_vdw_repulsion_cutoff=boost::lexical_cast<double>(parameter[5]);
		_grid_size=boost::lexical_cast<double>(parameter[6]);
	}*/
	/// Destructor
	~PMF04() {};
	
	/// Initialise the Fitness function with the protein.
	void initProtein(const ChainList& prot, const Vec3_t& as, Float_t r);
	
	/// Return the number of terms implemented in the fitness function.
	unsigned int TermsCount() { return _termscount; }
	
	/// Initialise the Fitness function with the ligand.
	void initLigand(const ChainList& lig);
	
	///calculate the vdw interaction for a single type at a position pos
	Float_t calcVDW(const Vec3_t& pos, unsigned int type);
	
	///calculate the statistical potential for a single type at a position pos
	Float_t calcSP(const Vec3_t& pos, unsigned int type);
	
	/// Evaluate the current configuration.
	void eval(const Vec3_t& pos, const Quat_t& ori, const vector<Float_t>& angle, 
	vector<Float_t>*& result);
	
	/// vdw radius
	static Float_t vdW_type[46];
		
private:
	
	// general fitness function informations
	// the vector with the scoring terms (results)
	enum {_termscount=1};
	vector<Float_t> _result;	
	
	static Float_t _parameter[18][31][60];
	static unsigned short _assignDistance[3601];	
	
	string PARAMETER;
	double _SF1, _SF2, _SF3, _vdw_repulsion_cutoff, _vdw_range,_grid_size;
	
	// protein informations
	// vdw protein atom neighbourhood list for vdW
	StaticNeighborhood _vdw_nh;
	// protein atoms for vdW interaction
	vector<unsigned int> _vdw_patom;
	// protein vdW parameter
	vector<Float_t> _pvdw;
	// vdw grid
	Grid _vdw_grid;
	
	// vdw protein atom neighbourhood list for sp (with hydrogen atoms)
	StaticNeighborhood _sp_nh;
	// protein atoms for statistical potential
	vector<unsigned int> _sp_patom;
	// statistical potential grid
	Grid _sp_grid;	
	
	//ligand informations
	// ligand atoms for vdW interaction
	vector<unsigned int> _vdw_latom;
	// ligand vdW parameter
	vector<unsigned int> _lvdw;
	// internal ligand vdw clashes penalty
	vector<vector<unsigned int> > _internal_vdw_neighborhood_list;
	vector<vector<Float_t> > _internal_vdw_neighborhood_list_d0;
	
	// return the vdw type of the atom
	unsigned int vdWType(AtomKey a);
	
	// returns all atoms which are more than n-1 bonds away in the same subgraph
	AtomList One_N_Neighborhood(AtomKey a, unsigned int n);	
	
};

}
#endif
