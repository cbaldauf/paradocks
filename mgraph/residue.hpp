#ifndef RESIDUE_H
#define RESIDUE_H

#include "mgraphtypes.hpp"
#include "residuetype.hpp"

namespace mgraph {

/**
\brief Residue Class
This is a container to group atoms in small organising units. It has a
parent Chain, a residue type (class RT) and a AtomList with all atoms
of the residue. It can also have a name and a id.
**/
class Residue {
friend ResidueKey rCreate(ChainKey p, RT rt);
friend void rDelete(ResidueKey k);
friend AtomKey aCreate(ResidueKey p, ELE e, Vec3_t c, GEO g);
friend void aDelete(AtomKey k);
public:
	/// returns the parent chain of the residue
	inline ChainKey Parent() const
	{
		#ifdef DEBUG2
		if (_parent.expired()) throw PhoenixError("Residue::Parent: can't get allready expired ChainKey in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
		#endif
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Residue::Parent: can't get parent of invalid residue.");
		ChainKey parent(_parent.lock());
		if (!parent->Valid()) throw PhoenixError("Residue::Parent: parent is invalid.");
		return parent;
		#else
		return _parent;
		#endif
	}
	/// returns the residue type of the residue
	inline RT Type() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Residue::RT: can't get type of invalid residue.");
		#endif
		return _rt;
	}
	/// sets the residue type of the residue
	inline void Type(RT rt)
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Residue::RT: can't assign type to invalid residue.");
		#endif
		_rt=rt;
	}
	/// returns the atoms of the residue
	inline AtomList Atoms() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Residue::Atoms: can't get atoms of invalid residue.");
		#endif
		return _atm;
	}
	/// returns the name of the residue
	inline string Name() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Residue::Name: can't get name of invalid residue.");
		#endif
		return _name;
	}
	/// sets the name of the residue
	inline void Name(string n)
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Residue::Name: can't assign name to invalid residue.");
		#endif
		_name=n;
	}
	/// returns the id of the residue
	inline unsigned int Id() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Residue::Id: can't get id of invalid residue.");
		#endif
		return _id;
	}
	/// sets the id of the residue
	inline void Id(unsigned int i)
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Residue::Id: can't assign id to invalid residue.");
		#endif
		_id=i;
	}
	#ifdef DEBUG
	/// returns a validation flag, only with -DDEBUG
	inline bool Valid() const
	{
		return _valid;
	}
	#endif
	// public destructor for smart pointer
	~Residue() {}
private:
	// only rCreate is allowed to construct residues
	Residue(ChainKey c, RT rt) : _parent(c), _rt(rt), _id(0)
	#ifdef DEBUG
	, _valid(false)
	#endif
	{}
	// the parent chain
	#ifdef DEBUG
	WeakChainKey _parent;
	#else
	ChainKey _parent;
	#endif
	// the residue type
	RT _rt;
	// the atomlist
	AtomList _atm;
	// the name	
	string _name;
	// the id
	unsigned int _id;
	#ifdef DEBUG
	// valid flag
	bool _valid;
	#endif
};

}
#endif
