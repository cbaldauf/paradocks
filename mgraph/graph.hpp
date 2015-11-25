#ifndef GRAPH_H
#define GRAPH_H

#include "mgraphtypes.hpp"
#include "chain.hpp"
#include "residue.hpp"
#include "atom.hpp"
#include "bond.hpp"

//TODO create a stable_unique for removal of duplicates(alBond and blAtoms)
//implement WriteMOL2
//aLPPos make it right or delete it
//delete lonepairs in CorrectGraphByTopology

namespace mgraph {

/**
\defgroup GSMF Graph/Tree Structure Modification Functions
These functions can be used to modify the hierarchical tree and the molecular graph.
**/

// chain functions
/// Create a new chain.
/** \ingroup GSMF
Create a new chain with name n and return the new ChainKey.
**/
ChainKey cCreate(string n);
/// Delete a chain.
/** \ingroup GSMF
Delete chain c and all of its children.
**/
void cDelete(ChainKey c);
/// Delete several chains.
/** \ingroup GSMF
Delete all chains in ChainList c and all of its children.
**/
void clDelete(const ChainList& cl);
// residue functions
/// Create a new residue.
/** \ingroup GSMF
Create a new residue of type rt at the end of chain p and return the new ResidueKey.
**/
ResidueKey rCreate(ChainKey p, RT rt);
/// Delete a residue.
/** \ingroup GSMF
Delete residue r and all of its children.
**/
void rDelete(ResidueKey r);
/// Delete several residues.
/** \ingroup GSMF
Delete all residues in ResidueList rl and all of its children.
**/
void rlDelete(const ResidueList& rl);
// atom functions
/// Create a new atom.
/** \ingroup GSMF
Create a new atom with element e, coordinates c and geometrie g at the end of residue p and
return the new AtomKey.
**/
AtomKey aCreate(ResidueKey p, ELE e, Vec3_t c, GEO g);
/// Delete an atom.
/** \ingroup GSMF
Delete atom a and all of its bonds.
**/
void aDelete(AtomKey a);
/// Delete several atoms.
/** \ingroup GSMF
Delete all atoms in AtomList al and all of its bonds.
**/
void alDelete(const AtomList& al);
// bond functions
/// Create a new bond.
/** \ingroup GSMF
Create a new bond with bond type bt from atom a1 to atom a2 and return the new BondKey.
If a bond between this atoms already exists it will be overwriten.
**/
BondKey bCreate(AtomKey a1, AtomKey a2, BT bt);
/// Delete a bond.
/** \ingroup GSMF
Delete bond b.
**/
void bDelete(BondKey b);
/// Delete several bonds.
/** \ingroup GSMF
Delete all bonds in BondList bl.
**/
void blDelete(const BondList& bl);



/**
\defgroup GSQF Graph/Tree Structure Query Functions
These functions can be used to get elements/parents/children of the 
hierarchical tree and the edges/nodes of the molecular graph.
**/

// chain functions
/// Return all chains.
/** \ingroup GSQF
Return all chains in the order they are stored in the global container.
**/
ChainList Chains();
/// Return the residues of a chain.
/** \ingroup GSQF
Return the residues of chain c in the order they are stored in the chain.
**/
ResidueList cResidues(ChainKey c);
/// Return the residues of several chains.
/** \ingroup GSQF
Return the residues of all chains in cl. The residues are simply concatenated.
This list is not uniq if one chain appears twice.
**/
ResidueList clResidues(const ChainList& cl);
/// Return the atoms of a chain.
/** \ingroup GSQF
Return the atoms of chain c by concatenating the atoms of all residue in c.
**/
AtomList cAtoms(ChainKey c);
/// Return the atoms of several chains.
/** \ingroup GSQF
Return the atoms of all chains in cl. The atoms are simply concatenated.
This list is not uniq if one chain appears twice.
**/
AtomList clAtoms(const ChainList& cl);
// residue functions
/// Return all residues.
/** \ingroup GSQF
Return all residues in the order they are stored in the chains. This is done by
concatenating all residues of all chains.
**/
ResidueList Residues();
/// Return the parent chain.
/** \ingroup GSQF
Return the parent chain.
**/
ChainKey rChain(ResidueKey r);
/// Return the parent chain of several residues.
/** \ingroup GSQF
Return the parent chain of several residues. The result will have one chain for every residue in cl.
**/
ChainList rlChain(const ResidueList& rl);
/// Return the atoms of a residue.
/** \ingroup GSQF
Return the atoms of residue r in the order they are stored in the residue.
**/
AtomList rAtoms(ResidueKey r);
/// Return the atoms of several residues.
/** \ingroup GSQF
Return the atoms of all residues in rl. The atoms are simply concatenated.
This list is not uniq if one residue appears twice.
**/
AtomList rlAtoms(ResidueList rl);
// atom functions
/// Return all atoms.
/** \ingroup GSQF
Return all atoms.This is done by concatenating all atoms of all residues of all chains.
**/
AtomList Atoms();
/// Return the parent chain of the parent residue.
/** \ingroup GSQF
Return the parent chain of the parent residue.
**/
ChainKey aChain(AtomKey a);
/// Return the parent chain of several atoms.
/** \ingroup GSQF
Return the parent chain of all the parent residues of all atoms in al. 
The result will have one chain for every atom in al.
**/
ChainList alChains(const AtomList& al);
/// Return the parent residue.
/** \ingroup GSQF
Return the parent residue.
**/
ResidueKey aResidue(AtomKey a);
/// Return the parent residues of several atoms.
/** \ingroup GSQF
Return the parent residues of all atoms in al. The result will have one 
residue for every atom in al.
**/
ResidueList alResidues(const AtomList& al);
/// Return the bonds of an atom.
/** \ingroup GSQF
Return the bonds of atom a.
**/
BondList aBonds(AtomKey a);
/// Return the bonds of several atoms.
/** \ingroup GSQF
Return the bonds of all atoms in al. The result is uniq.
**/
BondList alBonds(const AtomList& al);
/// Return the adjacent atoms of a.
/** \ingroup GSQF
Return the adjacent atoms of a. The order of the result is the same like the result from aBonds.
**/
AtomList aAdjAtoms(AtomKey a);
/// Return the bond between a1 and a2.
/** \ingroup GSQF
If atom a1 and atom a2 form a bond it will be in the result. Otherwise it will be empty.
**/
BondList aaFormBond(AtomKey a1, AtomKey a2);
/// Return the adjacent atom of a connected by b.
/** \ingroup GSQF
Return the adjacent atom of a connected by b. If there is no atom connected to a by b
the result will be empty.
**/
AtomList abAdjAtom(AtomKey a1, BondKey b);
// bond functions
/// Return all bonds.
/** \ingroup GSQF
Return all bonds.
**/
BondList Bonds();
/// Return the atoms of the bond.
/** \ingroup GSQF
Return the atoms of the bond.
**/
AtomList bAtoms(BondKey b);
/// Return the atoms of several bonds.
/** \ingroup GSQF
Return the atoms of the several bonds. The result is uniq.
**/
AtomList blAtoms(const BondList& bl);
/// Return true if bond b is in a ring.
/** \ingroup GSQF
Return true if bond b is in a ring. The size of the ring doesn't matter.
**/
bool isRingBond(BondKey b);
/// Return true if bond b is in a ring with x atoms.
/** \ingroup GSQF
Return true if bond b is in a ring with x atoms. Note: b can also be in a ring of different size.
**/
bool isXRingBond(BondKey b, unsigned int x);
/// Return true if atom a is in a ring.
/** \ingroup GSQF
Return true if atom a is in a ring. The size of the ring doesn't matter.
**/
bool isRingAtom(AtomKey a);
/// Return true if atom a is in a ring with x atoms.
/** \ingroup GSQF
Return true if atom a is in a ring with x atoms. Note: a can also be in a ring of different size.
**/
bool isXRingAtom(AtomKey a, unsigned int x);
/// Return all atoms which are connected to atom a.
/** \ingroup GSQF
Return all atoms which are connected to atom a.
**/
AtomList aExtendConnected(AtomKey a);
/// Return all atoms which are connected to any atom in al.
/** \ingroup CGSQF
Return all atoms which are connected to any atom in al.
**/
AtomList alExtendConnected(AtomList al);



/**
\defgroup CKQF Chemical Knowledge Query Functions
This functions can be used to get chemical features.
**/

/// Return the number of bonded heavy atoms.
/** \ingroup CKQF
Return the number of bonded heavy atoms.
**/
unsigned int aHeavyCount(AtomKey a);
/// Return the number of bonded hydrogen.
/** \ingroup CKQF
Return the number of bonded hydrogen.
**/
unsigned int aHCount(AtomKey a);
/// Return the number of lonepairs(explicit).
/** \ingroup CKQF
Return the number of lonepairs.
**/
unsigned int aLPCount(AtomKey a);
/// Return true if not all adjacent atoms C or H.
/** \ingroup CKQF
Return true if at least one of the adjacent atoms is not C or H.
Otherwise return false.
**/
bool aPolar(AtomKey a);
/// Return true if a is a metall atom.
/** \ingroup CKQF
Return true if a is a metall. Matalls are:
Li, Na, Mg, Al, K, Ca, Mn, Fe, Co, Ni, Cu, Zn
**/
bool aMetall(AtomKey a);

/// Return a vector with the positions of missing lonepairs according to geometry.
/** \ingroup CKQF
Return a vector with Coordinates for missing lone pairs, so that the lonepairs are as
far as possible away from the atoms and each other.
**/
//vector<Vec3_t> aLPPos(AtomKey a);

/// Fix small problems in the graph, so that type perception works well.
/** \ingroup CKQF
This function only uses connection informations. It depends on correct hydrogens 
and knows nothing about anorganic compounds. \n
Procedure: \n
1. delete lone pairs \n
2. remove all bonds from metals (Li, Na, Mg, Al, K, Ca, Fe, Zn, Pb, U) \n
3. atoms with no bonds are interpreted as ions \n
4. simple atoms with 1 bond are examined \n
5. simple atoms with 4 bond are examined \n
6. hydroxy/thiol ether/thioether group \n
7. carboxylic acid and nitro group \n
8. amide group and peptide bonds(but not cyclic) \n
9. guanidino group
**/
void CorrectGraphByTopology(ChainList cl);


/**
\defgroup SAF Structural Analysis Functions
Functions to measure and compare molecules.
**/



/**
\brief MOL Class
Class MOL can be used to copy molecules.
\ingroup GSMF
**/
class MOL
{
public:
	/**
	\brief Make a copy of all chains and save the state in internal variables.
	\param cl the chains to copy.
	**/
	void extractMOL(const ChainList& cl);
	/**
	\brief Recreate all objects from internal variables.
	**/
	ChainList createMOL();
private:
	struct CHAIN
	{
		CHAIN() {}
		CHAIN(string name, string tag) : _name(name), _tag(tag) {}
		string _name, _tag;
	};
	
	struct RESIDUE
	{
		RESIDUE() {}
		RESIDUE(RT rt, string name, unsigned int id) : _rt(rt), _name(name), _id(id) {}
		RT _rt;
		string _name;
		unsigned int _id;
	};
	
	struct ATOM
	{
		ATOM() {}
		ATOM(ELE ele, GEO geo, Vec3_t coord, string name, unsigned int fftype, Float_t charge, Float_t fcharge, string mol2type) :
		_ele(ele), _geo(geo), _coord(coord), _name(name), _fftype(fftype), _charge(charge), _fcharge(fcharge), _mol2type(mol2type) {}
		ELE _ele;
		GEO _geo;
		Vec3_t _coord;
		string _name;
		unsigned int _fftype;
		Float_t _charge, _fcharge;
		string _mol2type;
	};
	
	struct BOND
	{
		BOND() {}
		BOND(unsigned int a1, unsigned int a2, BT bt, unsigned int fftype, unsigned int a1pos, unsigned int a2pos) : 
		_a1(a1), _a2(a2), _bt(bt), _fftype(fftype), _a1pos(a1pos), _a2pos(a2pos) {}
		unsigned int _a1, _a2;
		BT _bt;
		unsigned int _fftype, _a1pos, _a2pos;
	};
	
	vector<CHAIN> _c;
	vector<vector<unsigned int> > _cres;
	vector<RESIDUE> _r;
	vector<vector<unsigned int> > _ratoms;
	vector<ATOM> _a;
	vector<BOND> _b;
};



///////////////////
// io functions
/////////////////
/// read a Tripos mol2 file
ChainList ReadMOL2(string file);
///////////////////////
// end io functions
/////////////////////


////////////////////////////////
// graph propertie functions
//////////////////////////////
/// return the distance in bonds between atom a1 and a2
int GraphDistance(AtomKey a1, AtomKey a2);
/// search for subgraphs matching query
//void GraphMatch(GraphQuery query);
////////////////////////////////////
// end graph propertie functions
//////////////////////////////////

/// return the atom adjacent to a1, connected by b
AtomKey aAdjacentAtom(AtomKey a1, BondKey b);



}
#endif

