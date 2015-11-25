#ifndef BOND_H
#define BOND_H

#include "mgraphtypes.hpp"

// TODO remove MOL2Type

namespace mgraph {

/**
\brief Bond Type Class
This class has a nested enum for different bondtypes. Supported bondtypes are: \n
BT::s - single bond \n
BT::re - bond involved in resonance stabilized system eg. carboxylate, 	peptide but not aromatic \n
BT::ar - aromatic bond \n
BT::d - double bond \n
BT::t - triple bond \n
BT::UNK - unknown
\ingroup graph
**/
class BT {
public:
	/// the nested enum for bondtypes
	enum _BT
	{
		s,	// single bond
		re,	// bond involved in resonance stabilized system
		ar,	// aromatic bond
		d,	// double bond
		t,	// triple bond
		UNK
	};
	// default construct to UNK
	BT() : _bt(UNK) {}
	/// constructor from enum
	BT(_BT b) : _bt(b) {}
	/// constructor from string
	BT(string s)
	{
		map<string, BT::_BT>::const_iterator i(_convertstr2bt.find(s));
		#ifdef ERROR
		if (i == _convertstr2bt.end()) throw PhoenixError("BT::BT: can't create BT from string: "+s);
		#endif
		_bt=i->second;
	}
	/// destructor
	~BT() {}
	/// relational operator==
	inline bool operator== (BT t) const { return (_bt==t._bt); }
	/// relational operator==
	inline bool operator== (_BT t) const {	return (_bt==t); }
	/// relational operator!=
	inline bool operator!= (BT t) const { return !operator==(t); }
	/// relational operator!=
	inline bool operator!= (_BT t) const { return !operator==(t); }
	/// cast from BT to unsigned int
	inline operator unsigned int() const {	return (unsigned int)_bt; }
	/// cast from BT to string
	inline operator string() const { return _convertbt2str[_bt]; }

private:
	_BT _bt;
	static const string _convertbt2str[];
	static const map<string, _BT> _convertstr2bt;
};

/// operator << for class BT
inline ostream& operator << (ostream& os, const BT& bt)
{
	return os<<string(bt);
}

/**
\brief Bond Class
Class Bond represents the edges in a molecular graph. A Bond has two AtomKeys for the atoms,
which are connected by the bond, and a bondtype (class BT). The atoms are named Atom1 and Atom2.
Its also possible to store a forcefield type (unsigned int) in the bond.
\ingroup graph
**/
class Bond {
friend BondKey bCreate(AtomKey a1, AtomKey a2, BT bt);
friend void bDelete(BondKey k);
public:
	/// returns Atom1
	inline AtomKey Atom1() const
	{
		#ifdef DEBUG2
		if (_a1.expired()) throw PhoenixError("Bond::Atom1: can't get allready expired AtomKey in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
		#endif
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Bond::Atom1: can't get Atom1 of invalid bond.");
		return _a1.lock();
		#else
		return _a1;
		#endif
	}
	/// returns Atom2
	inline AtomKey Atom2() const
	{
		#ifdef DEBUG2
		if (_a2.expired()) throw PhoenixError("Bond::Atom2: can't get allready expired AtomKey in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
		#endif
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Bond::Atom2: can't get Atom2 of invalid bond.");
		return _a2.lock();
		#else
		return _a2;
		#endif
	}
	/// returns the type of the bond
	inline BT Type() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("BT::Type: can't get type of invalid bond.");
		#endif
		return _bt;
	}
	/// sets the type of the bond
	inline void Type(BT bt)
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("BT::Type: can't assign type to invalid bond.");
		#endif
		_bt=bt;
	}
	/// returns the forcefield type
	inline unsigned int FFType() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("BT::FFType: can't get FF-type of invalid bond.");
		#endif
		return _fftype;
	}
	/// sets the forcefield type
	inline void FFType(unsigned int n)
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("BT::FFType: can't assign FF-type to invalid bond.");
		#endif
		_fftype=n;
	}
	#ifdef DEBUG
	/// returns a validation flag, only with -DDEBUG
	inline bool Valid() const
	{
		return _valid;
	}
	#endif
	
	// this should be removed
	inline string MOL2Type() const { return _mol2type; }
	inline void MOL2Type(string t) { _mol2type=t; }
	// public destructor for smart pointer
	~Bond() {}
private:
	// only bCreate is allowed to construct bonds
	Bond(AtomKey a1, AtomKey a2, BT bt) : _bt(bt)
	#ifdef DEBUG
	, _valid(false)
	#endif
	{
		#ifdef ERROR
		if (a1==a2) throw PhoenixError("Bond::Bond: can't construct selfedge");
		#endif
		_a1=a1; _a2=a2;
	}
	// Atom1 and Atom2
	#ifdef DEBUG
	WeakAtomKey _a1;
	WeakAtomKey _a2;
	#else
	AtomKey _a1;
	AtomKey _a2;
	#endif
	// bondtype
	BT _bt;
	// forcefield type
	unsigned int _fftype;
	#ifdef DEBUG
	// valid flag
	bool _valid;
	#endif
	
	// mol2 type
	string _mol2type;
};

}
#endif
