#ifndef RESIDUETYPE_H
#define RESIDUETYPE_H

#include "mgraphtypes.hpp"

namespace mgraph {

/**
\brief Residue Type Class
This class has a nested enum for different residues types. Supported residues types are: \n
standard as \n
ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, \n
PRO, SER, THR, TRP, TYR, VAL \n
modified standard as \n
ASZ, ASH - neutral ASP \n
GLZ, GLH - neutral GLU \n
HID - neutral HIS, H at ND1, converted to HIS \n
HIE - neutral HIS, H at NE2 \n
HIP - positive HIS \n
LYZ, LYN - neutral LYS \n
CYM - negative CYS \n
CYX - half cystine, CYS in disulfide bridges \n
TYM - negative TYR \n
N and C terminus \n
AMN - N-terminus ammonium \n
AMI - N-terminus amin \n
CXL - C-terminus carboxylate \n
CXC - C-terminus carbon acid \n
UNK - unspecified type
\ingroup graph
**/
class RT {
public:
	/// the nested enum for residues types
	enum _RT
	{
		UNK,
		ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE,
		PRO, SER, THR, TRP, TYR, VAL,

		ASZ, GLZ, HIE, HIP, LYZ, CYM, CYX, TYM,
		
		AMN, AMI, CXL, CXC
	};	
	
	// default construct to UNK
	RT() : _rt(UNK) {}
	/// constructor from enum
	RT(_RT rt) : _rt(rt) {}
	/// constructor from string
	RT(string s)
	{
		map<string, RT::_RT>::const_iterator i(_convertstr2rt.find(s));
		#ifdef ERROR
		if (i == _convertstr2rt.end()) cerr<<"RT::RT: can't create RT from string: "<<s<<endl;
		#endif
		_rt=i->second;
	}
	/// destructor
	~RT() {}
	/// relational operator==
	inline bool operator== (RT rt) const {	return (_rt==rt._rt); }
	/// relational operator==
	inline bool operator== (_RT rt) const { return (_rt==rt); }
	/// relational operator==
	inline bool operator!= (RT rt) const { return !operator==(rt); }
	/// relational operator==
	inline bool operator!= (_RT rt) const { return !operator==(rt); }
	/// cast from RT to unsigned int
	inline operator unsigned int() const {	return (unsigned int)_rt; }
	/// cast from RT to string
	inline operator string() const { return _convertrt2str[_rt]; }
	
private:
	_RT _rt;
	static const string _convertrt2str[];
	static const map<string, RT::_RT> _convertstr2rt;
};

/// operator << for class RT
inline ostream& operator << (ostream& os, const RT& rt)
{
	return os<<(string)rt;
}

}
#endif
