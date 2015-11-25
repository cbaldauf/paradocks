#include "graph.hpp"

#include <boost/spirit.hpp>
#include <boost/spirit/iterator.hpp>
 
namespace mgraph {

using namespace boost;
using namespace boost::spirit;

typedef char char_t;
typedef file_iterator<char_t> fiterator_t;
typedef scanner<fiterator_t> fscanner_t;
typedef rule<fscanner_t> frule_t;

//! sybyl atom types parser
//! http://www.tripos.com/mol2/atom_types.html#11465
struct sybylatomtype : symbols<boost::tuple<ELE, GEO, string> >
{
sybylatomtype()
{
	add
	("LP", boost::tuple<ELE, GEO, string>(ELE(ELE::LP), GEO(GEO::none), "LP"))
	("H.spc", boost::tuple<ELE, GEO, string>(ELE(ELE::H), GEO(GEO::none), "H.spc"))
	("H.t3p", boost::tuple<ELE, GEO, string>(ELE(ELE::H), GEO(GEO::none), "H.t3p"))
	("H", boost::tuple<ELE, GEO, string>(ELE(ELE::H), GEO(GEO::none), "H"))
	("C.3", boost::tuple<ELE, GEO, string>(ELE(ELE::C), GEO(GEO::tet), "C.3"))
	("C.ar", boost::tuple<ELE, GEO, string>(ELE(ELE::C), GEO(GEO::tri), "C.ar"))		
	("C.2", boost::tuple<ELE, GEO, string>(ELE(ELE::C), GEO(GEO::tri), "C.2"))
	("C.1", boost::tuple<ELE, GEO, string>(ELE(ELE::C), GEO(GEO::lin), "C.1"))
	("C.cat", boost::tuple<ELE, GEO, string>(ELE(ELE::C), GEO(GEO::tri), "C.cat"))
	("Du.C", boost::tuple<ELE, GEO, string>(ELE(ELE::C), GEO(GEO::none), "Du.C"))
	("N.3", boost::tuple<ELE, GEO, string>(ELE(ELE::N), GEO(GEO::tet), "N.3"))
	("N.4", boost::tuple<ELE, GEO, string>(ELE(ELE::N), GEO(GEO::tet), "N.4"))
	("N.pl3", boost::tuple<ELE, GEO, string>(ELE(ELE::N), GEO(GEO::tri), "N.pl3"))
	("N.ar", boost::tuple<ELE, GEO, string>(ELE(ELE::N), GEO(GEO::tri), "N.ar"))
	("N.am", boost::tuple<ELE, GEO, string>(ELE(ELE::N), GEO(GEO::tri), "N.am"))
	("N.2", boost::tuple<ELE, GEO, string>(ELE(ELE::N), GEO(GEO::tri), "N.2"))
	("N.1", boost::tuple<ELE, GEO, string>(ELE(ELE::N), GEO(GEO::lin), "N.1"))
	("O.3", boost::tuple<ELE, GEO, string>(ELE(ELE::O), GEO(GEO::tet), "O.3"))
	("O.co2", boost::tuple<ELE, GEO, string>(ELE(ELE::O), GEO(GEO::tri), "O.co2"))
	("O.2", boost::tuple<ELE, GEO, string>(ELE(ELE::O), GEO(GEO::tri), "O.2"))
	("O.spc", boost::tuple<ELE, GEO, string>(ELE(ELE::O), GEO(GEO::tri), "O.spc"))
	("O.t3p", boost::tuple<ELE, GEO, string>(ELE(ELE::O), GEO(GEO::tri), "O.t3p"))
	("S.3", boost::tuple<ELE, GEO, string>(ELE(ELE::S), GEO(GEO::tet), "S.3"))
	("S.2", boost::tuple<ELE, GEO, string>(ELE(ELE::S), GEO(GEO::tri), "S.2"))
	("S.o", boost::tuple<ELE, GEO, string>(ELE(ELE::S), GEO(GEO::tet), "S.o"))
	("S.o2", boost::tuple<ELE, GEO, string>(ELE(ELE::S), GEO(GEO::tet), "S.o2"))
	("P.3", boost::tuple<ELE, GEO, string>(ELE(ELE::P), GEO(GEO::tet), "P.3"))	
	("Du", boost::tuple<ELE, GEO, string>(ELE(ELE::DU), GEO(GEO::none), "Du"))
	("Any", boost::tuple<ELE, GEO, string>(ELE(ELE::DU), GEO(GEO::none), "Any"))
	// for correct parsing
	("He", boost::tuple<ELE, GEO, string>(ELE(ELE::He), GEO(GEO::UNK), "He"))
	("Hg", boost::tuple<ELE, GEO, string>(ELE(ELE::Hg), GEO(GEO::UNK), "Hg"))
	("Ho", boost::tuple<ELE, GEO, string>(ELE(ELE::Ho), GEO(GEO::UNK), "Ho"))
	("F", boost::tuple<ELE, GEO, string>(ELE(ELE::F), GEO(GEO::none), "F"))
	("Cl", boost::tuple<ELE, GEO, string>(ELE(ELE::Cl), GEO(GEO::none), "Cl"))
	("Br", boost::tuple<ELE, GEO, string>(ELE(ELE::Br), GEO(GEO::none), "Br"))
	("I", boost::tuple<ELE, GEO, string>(ELE(ELE::I), GEO(GEO::none), "I"))
	("Fe", boost::tuple<ELE, GEO, string>(ELE(ELE::Fe), GEO(GEO::none), "Fe"))
	;
}
} sybylatom_p;


struct atomicsymbol : symbols<ELE>
{
atomicsymbol()
{
	add
	("Li", ELE(ELE::Li))
	("Be", ELE(ELE::Be))
	("B",  ELE(ELE::B))
	("C",  ELE(ELE::C))
	("N",  ELE(ELE::N))
	("O",  ELE(ELE::O))
	("F",  ELE(ELE::F))
	("Ne", ELE(ELE::Ne))
	("Na", ELE(ELE::Na))
	("Mg", ELE(ELE::Mg))
	("Al", ELE(ELE::Al))
	("Si", ELE(ELE::Si))
	("P",  ELE(ELE::P))
	("S",  ELE(ELE::S))
	("Cl", ELE(ELE::Cl))
	("Ar", ELE(ELE::Ar))
	("K",  ELE(ELE::K))
	("Ca", ELE(ELE::Ca))
	("Sc", ELE(ELE::Sc))
	("Ti", ELE(ELE::Ti))
	("V",  ELE(ELE::V))
	("Cr", ELE(ELE::Cr))
	("Mn", ELE(ELE::Mn))
	("Fe", ELE(ELE::Fe))
	("Co", ELE(ELE::Co))
	("Ni", ELE(ELE::Ni))
	("Cu", ELE(ELE::Cu))
	("Zn", ELE(ELE::Zn))
	("Ga", ELE(ELE::Ga))
	("Ge", ELE(ELE::Ge))
	("As", ELE(ELE::As))
	("Se", ELE(ELE::Se))
	("Br", ELE(ELE::Br))
	("Kr", ELE(ELE::Kr))
	("Rb", ELE(ELE::Rb))
	("Sr", ELE(ELE::Sr))
	("Y",  ELE(ELE::Y))
	("Zr", ELE(ELE::Zr))
	("Nb", ELE(ELE::Nb))
	("Mo", ELE(ELE::Mo))
	("Tc", ELE(ELE::Tc))
	("Ru", ELE(ELE::Ru))
	("Rh", ELE(ELE::Rh))
	("Pd", ELE(ELE::Pd))
	("Ag", ELE(ELE::Ag))
	("Cd", ELE(ELE::Cd))
	("In", ELE(ELE::In))
	("Sn", ELE(ELE::Sn))
	("Sb", ELE(ELE::Sb))
	("Te", ELE(ELE::Te))
	("I",  ELE(ELE::I))
	("Xe", ELE(ELE::Xe))
	("Cs", ELE(ELE::Cs))
	("Ba", ELE(ELE::Ba))
	("La", ELE(ELE::La))
	("Ce", ELE(ELE::Ce))
	("Pr", ELE(ELE::Pr))
	("Nd", ELE(ELE::Nd))
	("Pm", ELE(ELE::Pm))
	("Sm", ELE(ELE::Sm))
	("Eu", ELE(ELE::Eu))
	("Gd", ELE(ELE::Gd))
	("Tb", ELE(ELE::Tb))
	("Dy", ELE(ELE::Dy))
	("Er", ELE(ELE::Er))
	("Tm", ELE(ELE::Tm))
	("Yb", ELE(ELE::Yb))
	("Lu", ELE(ELE::Lu))
	("Hf", ELE(ELE::Hf))
	("Ta", ELE(ELE::Ta))
	("W",  ELE(ELE::W))
	("Re", ELE(ELE::Re))
	("Os", ELE(ELE::Os))
	("Ir", ELE(ELE::Ir))
	("Pt", ELE(ELE::Pt))
	("Au", ELE(ELE::Au))
	("Tl", ELE(ELE::Tl))
	("Pb", ELE(ELE::Pb))
	("Bi", ELE(ELE::Bi))
	("Po", ELE(ELE::Po))
	("At", ELE(ELE::At))
	("Rn", ELE(ELE::Rn))
	("Fr", ELE(ELE::Fr))
	("Ra", ELE(ELE::Ra))
	("Ac", ELE(ELE::Ac))
	("Th", ELE(ELE::Th))
	("Pa", ELE(ELE::Pa))
	("U",  ELE(ELE::U))
	("Np", ELE(ELE::Np))
	("Pu", ELE(ELE::Pu))
	("Am", ELE(ELE::Am))
	("Cm", ELE(ELE::Cm))
	("Bk", ELE(ELE::Bk))
	("Cf", ELE(ELE::Cf))
	("Es", ELE(ELE::Es))
	("Fm", ELE(ELE::Fm))
	("Md", ELE(ELE::Md))
	("No", ELE(ELE::No))
	("Lr", ELE(ELE::Lr))
	("Rf", ELE(ELE::Rf))
	;
}
} element_p;

struct sybylbondtype : symbols<boost::tuple<BT, string> >
{
sybylbondtype()
{
	add
	("1",  boost::tuple<BT, string> (BT(BT::s), "1") )
	("2",  boost::tuple<BT, string> (BT(BT::d), "2") )
	("3",  boost::tuple<BT, string> (BT(BT::t), "3") )
	("am", boost::tuple<BT, string> (BT(BT::re), "am") )
	("ar", boost::tuple<BT, string> (BT(BT::ar), "ar") )
	("du", boost::tuple<BT, string> (BT(BT::UNK), "du") )
	("un", boost::tuple<BT, string> (BT(BT::UNK), "un") )
        ;
}
} sybylbondtype_p;


struct sybylresidue : symbols<RT>
{
	sybylresidue()
	{
	add
	("ALA", RT(RT::ALA))
	("ARG", RT(RT::ARG))
	("ASN", RT(RT::ASN))
	("ASP", RT(RT::ASP))
	("CYS", RT(RT::CYS))
	("GLN", RT(RT::GLN))
	("GLU", RT(RT::GLU))
	("GLY", RT(RT::GLY))
	("HIS", RT(RT::HIS))
	("ILE", RT(RT::ILE))
	("LEU", RT(RT::LEU))
	("LYS", RT(RT::LYS))
	("MET", RT(RT::MET))
	("PHE", RT(RT::PHE))
	("PRO", RT(RT::PRO))
	("SER", RT(RT::SER))
	("THR", RT(RT::THR))
	("TRP", RT(RT::TRP))
	("TYR", RT(RT::TYR))
	("VAL", RT(RT::VAL))
	// modufied as
	("ASZ", RT(RT::ASZ))
	("GLZ", RT(RT::GLZ))
	("HID", RT(RT::HIS))
	("HIE", RT(RT::HIE))
	("HIP", RT(RT::HIP))
	("LYZ", RT(RT::LYZ))
	("CYM", RT(RT::CYM))
	("CYX", RT(RT::CYX))
	("TYM", RT(RT::TYM))
	// terminus
	("AMN", RT(RT::AMN))
	("AMI", RT(RT::AMI))
	("CXL", RT(RT::CXL))
	("CXC", RT(RT::CXC));
    }
} sybylresidue_p;


struct bondinfo_t
{
	bondinfo_t() : t(BT(BT::UNK)) {}
	unsigned int origin_atom_id;
	unsigned int target_atom_id;
	BT t;
	string mol2t;
};

struct atominfo_t
{
	atominfo_t() : t(ELE(ELE::DU)), g(GEO(GEO::UNK)) {}
	string name;
	Float_t x;
	Float_t y;
	Float_t z;
	ELE t;
	GEO g;
	unsigned int subst_id;
	bool hascharge;
	Float_t charge;
	string mol2t;
};


struct substinfo_t
{
	substinfo_t() : t(RT(RT::UNK)) {}
	string name;
	string chain;
	RT t;
};

struct parsed_info
{
	parsed_info (string fn) : filename(fn), mcount(0) {};
	unsigned int id;
	string mname;
	string filename;
	unsigned int mcount;
	ChainList chains;
	bondinfo_t binfo;
	map<unsigned int, bondinfo_t> bdata;
	atominfo_t ainfo;
	map<unsigned int, atominfo_t> adata;
	substinfo_t sinfo;
	map<unsigned int, substinfo_t> sdata;
};

class DuplicateBondKeyError : public runtime_error
{
public:
	DuplicateBondKeyError(unsigned int i) : 
		runtime_error("Error: bond_id " + lexical_cast<string>(i) + " occours more than one time.") { }
};

class DuplicateAtomKeyError : public runtime_error
{
public:
	DuplicateAtomKeyError(unsigned int i) : 
		runtime_error("Error: atom_id " + lexical_cast<string>(i) + " occours more than one time.") { }
};

class DuplicateSubstKeyError : public runtime_error
{
public:
	DuplicateSubstKeyError(unsigned int i) : 
		runtime_error("Error: subst_id " + lexical_cast<string>(i) + " occours more than one time.") { }
};


struct fMOL2 : public grammar<fMOL2>
{
	struct atom_charge_actor
	{
		atom_charge_actor(atominfo_t& ainfo) : _ainfo(ainfo) {}
		void operator()(Float_t  c) const
		{
			_ainfo.hascharge=true;
			_ainfo.charge=c;
		}
		atominfo_t& _ainfo;
	};
	
	struct sybylatom_actor
	{
		sybylatom_actor(atominfo_t& ainfo) : _ainfo(ainfo) {}			
		void operator()(boost::tuple<ELE, GEO, string> sa) const
		{
			_ainfo.t = sa.get<0>();
			_ainfo.g = sa.get<1>();
			_ainfo.mol2t = sa.get<2>();
		}
		atominfo_t& _ainfo;
	};
	
	struct atomline_actor
	{
		atomline_actor(parsed_info& p) : _p(p) {}			
		void operator()(fiterator_t, fiterator_t) const
		{
			map<unsigned int, atominfo_t>::iterator i(_p.adata.find(_p.id));
			if (i!=_p.adata.end()) throw DuplicateAtomKeyError(_p.id);
			else _p.adata[_p.id]=_p.ainfo;
			_p.ainfo=atominfo_t();
		}
		parsed_info& _p;
	};
	
	struct sybylbond_actor
	{
		sybylbond_actor(bondinfo_t& binfo) : _binfo(binfo) {}			
		void operator()(boost::tuple<BT, string> sb) const
		{
			_binfo.t = sb.get<0>();
			_binfo.mol2t = sb.get<1>();
		}
		bondinfo_t& _binfo;
	};
	
	struct bondline_actor
	{
		bondline_actor(parsed_info& p) : _p(p) {}			
		void operator()(fiterator_t, fiterator_t) const
		{
			map<unsigned int, bondinfo_t>::iterator i(_p.bdata.find(_p.id));
			if (i!=_p.bdata.end()) throw DuplicateBondKeyError(_p.id);
			else _p.bdata[_p.id]=_p.binfo;
			_p.binfo=bondinfo_t();
		}
		parsed_info& _p;
	};
	
	struct substline_actor
	{
		substline_actor(parsed_info& p) : _p(p) {}			
		void operator()(fiterator_t, fiterator_t) const
		{
			map<unsigned int, substinfo_t>::iterator i(_p.sdata.find(_p.id));
			if (i!=_p.sdata.end()) throw DuplicateSubstKeyError(_p.id);
			else _p.sdata[_p.id]=_p.sinfo;
			_p.sinfo=substinfo_t();
		}
		parsed_info& _p;
	};
	
	struct molecule_actor
	{
		molecule_actor(parsed_info& p) : _p(p) {}			
		void operator()(fiterator_t, fiterator_t) const
		{
			_p.mcount++;
			if ((_p.mname=="") || (_p.mname=="****"))
				_p.mname=_p.filename + " molecule " + lexical_cast<string>(_p.mcount);
			// a default fallback subst
			_p.sdata[0]=substinfo_t();
			// all new residue keys will be saved here
			map<unsigned int, ResidueKey> cr;
			// all new chain keys will be saved here, indexed by chain tag
			map<string, ChainKey> cc;
			// all new chain keys will be saved here
			map<unsigned int, AtomKey> ca;
			// loop over every atom and create residues and chains
			for (map<unsigned int, atominfo_t>::iterator i=_p.adata.begin();
				i!=_p.adata.end(); ++i)
			{
				// has the parsed subst_id an entry in _p.sdata?
				map<unsigned int, substinfo_t>::iterator 
					si(_p.sdata.find(i->second.subst_id));
				// if no substructure with this id was parsed, use a default substructure
				if (si==_p.sdata.end()) si=_p.sdata.find(0);
					
				// find the residue for the atom
				map<unsigned int, ResidueKey>::iterator ri(cr.find(si->first));
				if (ri==cr.end())
				{
					// residue needs to be created
					// check if the substructure belongs to a new chain
					if (si->second.chain=="****") si->second.chain="";
					map<string, ChainKey>::iterator ci(cc.find(si->second.chain));
					if (ci==cc.end())
					{
						// this chain is new
						cc[si->second.chain]=cCreate(_p.mname);
						cc[si->second.chain]->Tag(si->second.chain);
					}
					// now we can create the residue
					//cr[si->first]=rCreate(cc[si->second.chain],si->second.name, si->second.t);
					cr[si->first]=rCreate(cc[si->second.chain], si->second.t);
					cr[si->first]->Name(si->second.name);
					ri=cr.find(si->first);					
				}
				AtomKey ak;
				ak=aCreate(ri->second, i->second.t, Vec3_t(i->second.x, i->second.y, i->second.z), i->second.g);
				//ak->Geometry(i->second.g);
				ak->Name(i->second.name);
				ak->MOL2Type(i->second.mol2t);
				if (i->second.hascharge) ak->Charge(i->second.charge);
				ca[i->first]=ak;				
			}
			// loop over every bond
			map<unsigned int, bondinfo_t>::iterator i;
			for (map<unsigned int, bondinfo_t>::iterator i=_p.bdata.begin();
				i!=_p.bdata.end(); ++i)
			{
				map<unsigned int, AtomKey>::iterator key1, key2;
				key1=ca.find(i->second.origin_atom_id);
				key2=ca.find(i->second.target_atom_id);
				if ((key1==ca.end()) || (key2==ca.end()))
				{
					throw PhoenixError("Found incomplete bond.");
				}
				BondKey b=bCreate(key1->second, key2->second, i->second.t);
				b->MOL2Type(i->second.mol2t);
			}
			for (map<string, ChainKey>::iterator i=cc.begin();
				i!=cc.end(); ++i)
			{
				_p.chains.push_back(i->second);
			}
			_p.bdata.clear();
			_p.adata.clear();
			_p.sdata.clear();
			
		}
		parsed_info& _p;
	};

	template <typename ScannerT>
	struct definition
	{
		definition(fMOL2 const& self)
		{
			ignored_line = comment_p("#") | *blank_p >> eol_p;
			line = (*(anychar_p - eol_p)) >> eol_p;
			moleculertiline = str_p("@<TRIPOS>MOLECULE") >> eol_p;
			moleculerti = 
				moleculertiline >>
				// the name of the molecule
				*blank_p >> (*(anychar_p - eol_p))[assign_a(self._p.mname)] >> eol_p >>
				// next 5 lines are ignored
				repeat_p(5)[((*(anychar_p - eol_p)) >> eol_p)];
			atomline = 
				// atom_id
				*blank_p >> uint_p[assign_a(self._p.id)] >>
				// atom_name
				+blank_p >> (+(alnum_p|punct_p))[assign_a(self._p.ainfo.name)] >>
				// x y z
				+blank_p >> real_p[assign_a(self._p.ainfo.x)] >>
				+blank_p >> real_p[assign_a(self._p.ainfo.y)] >>
				+blank_p >> real_p[assign_a(self._p.ainfo.z)] >>
				// atom_type
				+blank_p >> (sybylatom_p[sybylatom_actor(self._p.ainfo)] 
					| element_p[assign_a(self._p.ainfo.t)]) >>
				// subst_id
				!(+blank_p >> uint_p[assign_a(self._p.ainfo.subst_id)] >>
				// subst_name(not needed)				
				!(+blank_p >> +(alnum_p|punct_p) >>
				// charge				
				!(+blank_p >> real_p[atom_charge_actor(self._p.ainfo)] >>
				// status_bit(not needed)
				!(+blank_p >> +(alnum_p|punct_p) )))) >>
				*(blank_p - eol_p) >> eol_p;
			atomrti =
				str_p("@<TRIPOS>ATOM") >> eol_p >>
				*(atomline[atomline_actor(self._p)] | ignored_line);
			bondline =
				// bond_id
				*blank_p >> uint_p[assign_a(self._p.id)] >>
				// origin_atom_id
				+blank_p >>	uint_p[assign_a(self._p.binfo.origin_atom_id)] >>
				// target_atom_id
				+blank_p >>	uint_p[assign_a(self._p.binfo.target_atom_id)] >>
				// bond_type
				+blank_p >>	sybylbondtype_p[sybylbond_actor(self._p.binfo)] >>
				// status_bit(not needed)
				!(+blank_p >> +(alnum_p|punct_p) ) >>
				*(blank_p - eol_p) >> eol_p;
			bondrti = 
				str_p("@<TRIPOS>BOND") >> eol_p >>
				*(bondline[bondline_actor(self._p)] | ignored_line);
			substline =
				// subst_id
				*blank_p >> uint_p[assign_a(self._p.id)] >>
				// subst_name
				+blank_p >> (+(alnum_p|punct_p))[assign_a(self._p.sinfo.name)] >>
				// root_atom(not needed)
				+blank_p >>	uint_p >>
				// subst_type(not needed)
				!(+blank_p >> +(alnum_p|punct_p) >>
				// dict_type(not needed)
				!(+blank_p >> uint_p >>
				// chain
				!(+blank_p >> (+(alnum_p|punct_p))[assign_a(self._p.sinfo.chain)] >>
				// subst_type
				!(+blank_p >> (sybylresidue_p[assign_a(self._p.sinfo.t)] | (+(alnum_p|punct_p))) >>
				// inter_bonds(not needed)
				!(+blank_p >> uint_p >>
				// status_bit(not needed)
				!(+blank_p >> +(alnum_p|punct_p) >>
				// comment
				!(+blank_p >> +(alnum_p|punct_p|blank_p) ))))))) >>
				*(blank_p - eol_p) >> eol_p;
			substrti = 
				str_p("@<TRIPOS>SUBSTRUCTURE") >> eol_p >>
				*(substline[substline_actor(self._p)] | ignored_line);
			rtiline = 
				str_p("@<TRIPOS>") >> +(alnum_p|punct_p) >> eol_p;
			otherrti = 
				(rtiline - moleculertiline) >>
				*(line -rtiline);
			molecule =
				moleculerti >>
				*(atomrti | bondrti | substrti | otherrti);
			first = *ignored_line >>
				+molecule[molecule_actor(self._p)];
		}
		frule_t ignored_line, line, moleculertiline, moleculerti, atomline;
		frule_t atomrti, bondline, bondrti, substline, substrti, rtiline;
		frule_t otherrti, molecule, first;
		frule_t const& start() const { return first; }
	};
	fMOL2(parsed_info& p) : _p(p) {}
	parsed_info& _p;
};

ChainList ReadMOL2(string s)
{
	file_iterator<> first(s);
	if (!first)
	{
		string emsg("Could not open file "+s+" for reading!");
		throw PhoenixError(emsg);
	}
	file_iterator<> last = first.make_end();
	cout << "Read file " << s << " as Tripos mol2 file." << endl;
	parsed_info p(s);
	fMOL2 MOL2_p(p);
	try
	{
		parse_info<fiterator_t> i=parse(first, last, MOL2_p);
		if (i.full)
		{
			#ifdef DEBUG
			cout << "Finished with no errors." << endl;
			#endif
		}
		else
		{
			cout << "parsing failed!" << endl;
			first=i.stop;
			cout << "Error in this text:-------" << endl;
			while ((first != last) && (first != i.stop+10))
			cout << *first++;
			cout << "end-----------------------" << endl;
			
		}
	}
	catch (DuplicateAtomKeyError e)
	{
		cout<< e.what() << endl;
	}
	catch (DuplicateBondKeyError e)
	{
		cout<< e.what() << endl;
	}
	catch (DuplicateSubstKeyError e)
	{
		cout<< e.what() << endl;
	}
	#ifdef DEBUG
	cout << "Created " << p.mcount << " molecule/s with " << p.chains.size() << " chains." << endl;
	#endif
	return p.chains;
}


}
