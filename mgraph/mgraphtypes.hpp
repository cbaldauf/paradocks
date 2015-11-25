#ifndef MGRAPHTYPES_H
#define MGRAPHTYPES_H

#include <stdexcept>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <iostream>
#include <numeric>
#include <math/osg/Vec3f>
#include <math/osg/Vec3d>
#include <math/osg/Quat>
#include <math/osg/Matrixf>
#include <math/osg/Matrixd>
#include <boost/lexical_cast.hpp>
#include <boost/assign.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/tuple/tuple.hpp>

#if ( defined(DEBUG) || defined(DEBUG2))
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#endif

/**
\defgroup graph Graph/Tree Types and Classes
**/

/**
\brief All functionality of the molecular graph is in namespace mgraph.
There are Chains, Residues, Atoms and Bonds. Chains have Residues, Residues have Atoms and Atoms
are connected by Bonds. Atoms and Bonds form a graph, while Chains, Residues and Atoms are 
organized in a tree structure. The objects itself are stored in hidden structures and manipulation
is possible by functions. Access to the members of the objects is provided through pointer. 
They are called ChainKey, ResidueKey, AtomKey and BondKey.
**/
namespace mgraph {
using namespace osg;
using namespace std;

// class forward declarations
class Chain;
class Residue;
class RT;
class Atom;
class ELE;
class GEO;
class Bond;
class BT;

/**
\brief errors class for mgraph
\ingroup graph
All errors thrown in namespace mgraph are of this type.
**/
class PhoenixError : public runtime_error
{
public:
	PhoenixError(string err) : runtime_error("PhoenixError: "+err) { }
};

/**
\brief floating point number type.
\ingroup graph
**/
typedef float Float_t;
/**
\brief general purpose 3d vector (osg::Vec3d Class).
\ingroup graph
**/
typedef Vec3f Vec3_t;
/**
\brief general purpose quaternion (osg::Quat Class).
\ingroup graph
**/
typedef Quat Quat_t;
/**
\brief general purpose 4x4 matrix (osg::Matrixd Class).
\ingroup graph
**/
typedef Matrixf Matrix_t;


#if !( defined(DEBUG) || defined(DEBUG2))
/**
\brief ChainKey is used to access the members of Chain and as result/input for functions.
\ingroup graph
**/
typedef Chain* ChainKey;
/**
\brief ResidueKey is used to access the members of Residue and as result/input for functions.
\ingroup graph
**/
typedef Residue* ResidueKey;
/**
\brief AtomKey is used to access the members of Atom and as result/input for functions.
\ingroup graph
**/
typedef Atom* AtomKey;
/**
\brief BondKey is used to access the members of Bond and as result/input for functions.
\ingroup graph
**/
typedef Bond* BondKey;
#else
//ChainKey is used to access the members of Chain and as result/input for functions.
//With -DDEBUG its an shared pointer.
typedef boost::shared_ptr<Chain> ChainKey;
typedef boost::weak_ptr<Chain> WeakChainKey;
//ResidueKey is used to access the members of Residue and as result/input for functions.
//With -DDEBUG its an shared pointer.
typedef boost::shared_ptr<Residue> ResidueKey;
typedef boost::weak_ptr<Residue> WeakResidueKey;
//AtomKey is used to access the members of Atom and as result/input for functions.
//With -DDEBUG its an shared pointer.
typedef boost::shared_ptr<Atom> AtomKey;
typedef boost::weak_ptr<Atom> WeakAtomKey;
//BondKey is used to access the members of Bond and as result/input for functions.
//With -DDEBUG its an shared pointer.
typedef boost::shared_ptr<Bond> BondKey;
typedef boost::weak_ptr<Bond> WeakBondKey;
#endif
/**
\brief container to store "ChainKey"s. Some functions expect this container for in/output.
\ingroup graph
**/
typedef vector<ChainKey> ChainList;
/**
\brief container to store "ResidueKey"s. Some functions expect this container for in/output.
\ingroup graph
**/
typedef vector<ResidueKey> ResidueList;
/**
\brief container to store "AtomKey"s. Some functions expect this container for in/output.
\ingroup graph
**/
typedef vector<AtomKey> AtomList;
/**
\brief container to store "BondKey"s. Some functions expect this container for in/output.
\ingroup graph
**/
typedef vector<BondKey> BondList;

/**
\brief operator<< for Vec3_t
\ingroup graph
**/
inline ostream& operator<<(ostream& os, const Vec3_t& v)
{
	return os<<v.x()<<" "<<v.y()<<" "<<v.z();
}

/**
\brief operator<< for Quat_t
\ingroup graph
**/
inline ostream& operator<<(ostream& os, const Quat_t& v)
{
	return os<<v.x()<<" "<<v.y()<<" "<<v.z()<<" "<<v.w();
}

/**
\brief normalize Quat_t
\ingroup graph
**/
inline void normalize(Quat_t& r)
{
	Quat_t::value_type norm = r.length();
	if (norm>0.0)
	{
		double inv = 1.0/norm;
		r._v[0] *= inv;
		r._v[1] *= inv;
		r._v[2] *= inv;
		r._v[3] *= inv;
	}
}

}

#endif /* PHOENIXTYPES_H */
