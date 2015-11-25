#ifndef CHAIN_H
#define CHAIN_H

#include "mgraphtypes.hpp"

namespace mgraph
{

/**
\brief Chain Class
This class is the top level in the molecular tree structure. It has a name,
a ResidueList with all residues and a tag.
\ingroup graph
**/
class Chain {
friend ChainKey cCreate(string n);
friend void cDelete(ChainKey k);
friend ResidueKey rCreate(ChainKey p, RT rt);
friend void rDelete(ResidueKey k);
public:
	/// returns the name of the chain
	inline string Name() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Chain::Name: can't get ResidueKey of invalid atom.");
		#endif
		return _name;
	}
	/// sets the name of the chain
	inline void Name(string n)
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Chain::Name: can't assign name to invalid chain.");
		#endif
		_name=n;
	}
	/// return the residues of the chain
	inline ResidueList Residues() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Chain::Residues: can't get residues of invalid chain.");
		#endif
		return _res;
	}
	/// returns the tag of the chain
	inline string Tag() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Chain::Tag: can't get tag of invalid chain.");
		#endif
		return _tag;
	}
	/// sets the tag of the chain
	inline void Tag(string t)
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Chain::Tag: can't assign tag to invalid chain.");
		#endif
		_tag=t;
	}
	#ifdef DEBUG
	/// returns a validation flag, only with -DDEBUG
	inline bool Valid() const
	{
		return _valid;
	}
	#endif
private:
	// only cCreate is allowed to construct chains
	Chain(string n) : _name(n)
	#ifdef DEBUG
	, _valid(false)
	#endif
	{ }
	// the name
	string _name;
	// the residuelist
	ResidueList _res;
	// the tag
	string _tag;
	#ifdef DEBUG
	// valid flag
	bool _valid;
	#endif
};

}
#endif
