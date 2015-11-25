#ifndef ATOM_H
#define ATOM_H

#include "mgraphtypes.hpp"

// TODO remove MOL2Type

namespace mgraph {

/**
\brief Element Class
This class has a nested enum for all elements.
\ingroup graph
**/
class ELE {
public:
	/// the nested enum for elements
	enum _ELE
	{
		LP,
		H, He,
		Li, Be, B, C, N, O, F, Ne,
		Na, Mg, Al, Si, P, S, Cl, Ar,
		K, Ca,
			Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn,
				Ga, Ge, As, Se, Br, Kr,
		Rb, Sr,
			Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd,
				In, Sn, Sb, Te, I, Xe,
		Cs, Ba,
			La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb,
			Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg,
				Tl, Pb, Bi, Po, At, Rn,
		Fr, Ra,
			Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No,
			Lr, Rf,
		DU, Any
	};
	
	// default construct to DU
	ELE() : _e(DU) {}
	/// constructor from enum
	ELE(_ELE e) : _e(e) {}
	/// constructor from atomic number
	ELE(unsigned int i)
	{
		#ifdef ERROR
		if (i>DU) throw PhoenixError("ELE::ELE: can't create ELE from number: "+boost::lexical_cast<string>(i));
		#endif
		_e=_ELE(i);
	}
	/// constructor from string
	ELE(string s)
	{
		map<string, _ELE>::const_iterator i(_convertstr2ele.find(s));
		#ifdef ERROR
		if (i == _convertstr2ele.end()) throw PhoenixError("ELE::ELE: can't create ELE from string: "+s);
		#endif
		_e=i->second;
	}
	/// destructor
	~ELE() {}
	/// operator==
	inline bool operator== (ELE e) const { return (_e==e._e); }
	/// operator==
	inline bool operator== (_ELE e) const { return (_e==e); }
	/// operator!=
	inline bool operator!= (ELE e) const { return !operator==(e); }
	/// operator!=
	inline bool operator!= (_ELE e) const { return !operator==(e); }
	/// cast from ELE to unsigned int
	inline operator unsigned int() const { return (unsigned int)_e; }
	/// cast from ELE to string
	inline operator string() const { return _convertele2str[_e]; }
	
private:
	_ELE _e;
	static const string _convertele2str[];
	static const map<string, ELE::_ELE> _convertstr2ele;
};

/// operator << for class ELE
inline ostream& operator << (ostream& os, const ELE& e)
{
	return os<<(string)e;
}

/**
\brief Geometry Class
This class has a nested enum for geometries. Supported geometries are: \n
GEO::none - atoms with no or one bond eg. H, Fl, Cl \n
GEO::lin - linear geometry for atoms with 2 bonds eg. sp carbon \n
GEO::tri - trigonal planar geometry eg. sp2 carbon/nitrogen, carboxylate oxygen(lonepairs also count) \n
GEO::tet - tetrahedral geometry eg. sp3 carbon but also hydroxy oxygen \n
GEO::bip - trigonal bipyramidal \n
GEO::oct - octahedral \n
GEO::UNK - unknown geometry
\ingroup graph
**/
class GEO {
public:
	/// the nested enum for geometries
	enum _GEO
	{
		none,
		lin,			//linear
		tri,			//trigonal
		tet,			//tetrahedral
		bip,			//trigonal bipyramidal
		oct,			//octahedral
		UNK,
	};
	
	// default construct to UNK
	GEO() : _g(UNK) {}
	/// constructor from enum
	GEO(_GEO g) : _g(g) {}
	/// constructor from string
	GEO(string s)
	{
		map<string, GEO::_GEO>::const_iterator i(_convertstr2g.find(s));
		if (i != _convertstr2g.end()) _g=i->second;
		else throw PhoenixError("GEO::GEO: can't create GEO from string: "+s);
	}
	/// relational operator==
	inline bool operator== (GEO g) const { return (_g==g._g); }
	/// relational operator==
	inline bool operator== (_GEO g) const { return (_g==g); }
	/// relational operator!=
	inline bool operator!= (GEO g) const { return !operator==(g); }
	/// relational operator!=
	inline bool operator!= (_GEO g) const { return !operator==(g); }
	/// cast from GEO to unsigned int
	inline operator unsigned int() const { return (unsigned int)_g; }
	/// cast from GEO to string
	inline operator string() const { return _convertg2str[_g]; }
	
private:
	_GEO _g;
	static const string _convertg2str[];
	static const map<string, GEO::_GEO> _convertstr2g;
};

/// operator << for class GEO
inline ostream& operator << (ostream& os, const GEO& g)
{
	return os<<(string)g;
}

/**
\brief Atom Class
Class Atom represents the nodes in the molecular graph and the last layer in the hirachical
tree. An Atom has an parent Residue, an Element (class ELE), a Geometry (class GEO), 3d
coordinates (Vec3_t) and a list of bonds. It can also have a name, a forcefield type, a
charge and a formal charge.
\ingroup graph
**/
class Atom {
friend AtomKey aCreate(ResidueKey p, ELE e, Vec3_t c, GEO g);
friend void aDelete(AtomKey k);
friend BondKey bCreate(AtomKey a1, AtomKey a2, BT bt);
friend void bDelete(BondKey k);
public:
	/// returns the parent ResidueKey
	inline ResidueKey Parent() const
	{
		#ifdef DEBUG2
		if (_parent.expired()) throw PhoenixError("Atom::Parent: can't get allready expired ResidueKey in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
		#endif
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::Parent: can't get ResidueKey of invalid atom.");
		return _parent.lock();
		#else
		return _parent;
		#endif
	}	
	/// returns the element of the atom
	inline ELE Element() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::Element: can't get element of invalid atom.");
		#endif
		return _ele;
	}
	/// sets the element of the atom
	inline void Element(ELE e)
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::Element: can't assign element to invalid atom.");
		#endif
		_ele=e;
	}
	/// returns the geometry of the atom
	inline GEO Geometry() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::Geometry: can't get geometry of invalid atom.");
		#endif
		return _geo;
	}
	/// sets the geometry of the atom
	inline void Geometry(GEO g)
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::Geometry: can't assign geometry to invalid atom.");
		#endif
		_geo=g;
	}
	/// returns the coordinates of the atom
	inline Vec3_t Coord() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::Coord: can't get coordinates of invalid atom.");
		#endif
		return _coord;
	}
	/// sets the coordinates of the atom
	inline void Coord(Vec3_t c)
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::Coord: can't assign coordinates to invalid atom.");
		#endif
		_coord=c;
	}
	inline BondList Bonds() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::Bonds: can't get bonds of invalid atom.");
		#endif
		return _bnd;
	}
	/// returns the name of the atom
	inline string Name() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::Name: can't get name of invalid atom.");
		#endif
		return _name;
	}
	/// sets the name of the atom
	inline void Name(string n)
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::Name: can't assign name to invalid atom.");
		#endif
		_name=n;
	}
	/// returns the forcefield type of the atom
	inline unsigned int FFType() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::FFType: can't get forcefield type of invalid atom.");
		#endif
		return _fftype;
	}
	/// sets the forcefield type of the atom
	inline void FFType(unsigned int t)
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::FFType: can't assign forcefield type to invalid atom.");
		#endif
		_fftype=t;
	}
	/// returns the charge of the atom	
	inline Float_t Charge() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::Charge: can't get charge of invalid atom.");
		#endif
		return _charge;
	}
	/// sets the charge of the atom
	inline void Charge(Float_t c)
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::Charge: can't assign charge to invalid atom.");
		#endif
		_charge=c;
	}
	/// returns the formal charge of the atom
	inline Float_t FormalCharge() const
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::FormalCharge: can't get formal charge of invalid atom.");
		#endif
		return _fcharge;
	}
	/// sets the formal charge of the atom
	inline void FormalCharge(Float_t c)
	{
		#ifdef DEBUG
		if (!_valid) throw PhoenixError("Atom::FormalCharge: can't assign formal charge to invalid atom.");
		#endif
		_fcharge=c;
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
	~Atom() {}
private:
	// only aCreate is allowed to construct atoms
	Atom(ResidueKey r, ELE e, Vec3_t c, GEO g) : _parent(r), _ele(e),
	_geo(g), _coord(c), _fftype(0), _charge(0), _fcharge(99)
	#ifdef DEBUG
	, _valid(false)
	#endif
	{}
	// the parent residue
	#ifdef DEBUG
	WeakResidueKey _parent;
	#else
	ResidueKey _parent;
	#endif
	// the element
	ELE _ele;
	// the geometry
	GEO _geo;
	// the coordinates
	Vec3_t _coord;
	// the bondlist
	BondList _bnd;
	// the atom name
	string _name;
	// atom ff-type
	unsigned int _fftype;
	// the charge of the atom
	Float_t _charge;
	// the formal charge of the atom
	Float_t _fcharge;
	#ifdef DEBUG
	// valid flag
	bool _valid;
	#endif
	
	// mol2 type
	string _mol2type;
};

}
#endif
