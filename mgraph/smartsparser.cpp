#include "graph.hpp"

#include "smartsparser.hpp"
#include <boost/spirit/include/classic.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <iostream>


using namespace boost;
using namespace boost::spirit::classic;		// entspricht dem alten boost::spirit 1.8.x
using namespace mgraph;


/// create a storage structure for SMARTS using BGL ----------------------
///

// property of vertices - in this case atom properties such as name, charged, number of hydrogens, negation, return sign
struct smart_at
{
	smart_at() : t(ELE(ELE::DU)), charge(99), num_of_hydrogen(99), negation(0), return_sign(0) {}
	string name;
	ELE t;
	unsigned int at_id;
	int charge;		/// set to default value --> i.e. if user do not set up a charge, then the charge never minds
	unsigned int num_of_hydrogen;
	bool negation;		/// set default value --> FALSE
	bool return_sign;	/// set default value --> FALSE
};
// property of edges - in this case bond properties such as name --> will be used in the BGL presentaion for edge properties

struct parsed_info
{
 	parsed_info(string fn) : smartstring(fn) {};
	string smartstring;
	smart_at ainfo;
// 	smart_bt binfo;
};
	
/// hier werden alle Atome abgespeichert
vector <smart_at> myAtoms;

/// save ring bonds: atomkey, bondtype (-, =, #, ~), id of bond (0-9)
vector <boost::tuple<unsigned, unsigned, unsigned> > ringbonds;

/// save all other bonds: atomkey, bondtype (-, =, #, ~), atomkey
vector <boost::tuple<unsigned, unsigned, unsigned> > bond;

/// save information about the number of atoms within different branches 
struct branch_content
{
	branch_content(int a=0, int b=0) :a(a), b(b) {} 
	unsigned int a;
	unsigned int b;
};

vector <branch_content> branch_info;
int test=0; /// to rememeber the atom before the branch opening (therefore the atoms within each branch are counted)
int atom_counter=-1; /// to set up the atom id

////////////////////////////////////////////////////////////////////////////
//
//  Semantic actions
//
////////////////////////////////////////////////////////////////////////////

/// create symbol table
struct atomicsymbol : symbols<ELE>
{
	atomicsymbol()
	{
		add
			("*", ELE(ELE::Any))	("He", ELE(ELE::He))	("Li", ELE(ELE::Li))	("Be", ELE(ELE::Be))	("B",  ELE(ELE::B))	("C",  ELE(ELE::C))	
			("N",  ELE(ELE::N))	("O",  ELE(ELE::O))	("F",  ELE(ELE::F))	("Ne", ELE(ELE::Ne))	("Na", ELE(ELE::Na))	("Mg", ELE(ELE::Mg))	
			("Al", ELE(ELE::Al))	("Si", ELE(ELE::Si))	("P",  ELE(ELE::P))	("S",  ELE(ELE::S))	("Cl", ELE(ELE::Cl))	("Ar", ELE(ELE::Ar))
			("K",  ELE(ELE::K))	("Ca", ELE(ELE::Ca))	("Sc", ELE(ELE::Sc))	("Ti", ELE(ELE::Ti))	("V",  ELE(ELE::V))	("Cr", ELE(ELE::Cr))	
			("Mn", ELE(ELE::Mn))	("Fe", ELE(ELE::Fe))	("Co", ELE(ELE::Co))	("Ni", ELE(ELE::Ni))	("Cu", ELE(ELE::Cu))	("Zn", ELE(ELE::Zn))	
			("Ga", ELE(ELE::Ga))	("Ge", ELE(ELE::Ge))	("As", ELE(ELE::As))	("Se", ELE(ELE::Se))	("Br", ELE(ELE::Br))	("Kr", ELE(ELE::Kr))
			("Rb", ELE(ELE::Rb))	("Sr", ELE(ELE::Sr))	("Y",  ELE(ELE::Y))	("Zr", ELE(ELE::Zr))	("Nb", ELE(ELE::Nb))	("Mo", ELE(ELE::Mo))
			("Tc", ELE(ELE::Tc))	("Ru", ELE(ELE::Ru))	("Rh", ELE(ELE::Rh))	("Pd", ELE(ELE::Pd))	("Ag", ELE(ELE::Ag))	("Cd", ELE(ELE::Cd))
			("In", ELE(ELE::In))	("Sn", ELE(ELE::Sn))	("Sb", ELE(ELE::Sb))	("Te", ELE(ELE::Te))	("I",  ELE(ELE::I))	("Xe", ELE(ELE::Xe))
			("Cs", ELE(ELE::Cs))	("Ba", ELE(ELE::Ba))	("La", ELE(ELE::La))	("Ce", ELE(ELE::Ce))	("Pr", ELE(ELE::Pr))	("Nd", ELE(ELE::Nd))
			("Pm", ELE(ELE::Pm))	("Sm", ELE(ELE::Sm))	("Eu", ELE(ELE::Eu))	("Gd", ELE(ELE::Gd))	("Tb", ELE(ELE::Tb))	("Dy", ELE(ELE::Dy))
			("Er", ELE(ELE::Er))	("Tm", ELE(ELE::Tm))	("Yb", ELE(ELE::Yb))	("Lu", ELE(ELE::Lu))	("Hf", ELE(ELE::Hf))	("Ta", ELE(ELE::Ta))
			("W",  ELE(ELE::W))	("Re", ELE(ELE::Re))	("Os", ELE(ELE::Os))	("Ir", ELE(ELE::Ir))	("Pt", ELE(ELE::Pt))	("Au", ELE(ELE::Au))
			("Tl", ELE(ELE::Tl))	("Pb", ELE(ELE::Pb))	("Bi", ELE(ELE::Bi))	("Po", ELE(ELE::Po))	("At", ELE(ELE::At))	("Rn", ELE(ELE::Rn))
			("Fr", ELE(ELE::Fr))	("Ra", ELE(ELE::Ra))	("Ac", ELE(ELE::Ac))	("Th", ELE(ELE::Th))	("Pa", ELE(ELE::Pa))	("U",  ELE(ELE::U))
			("Np", ELE(ELE::Np))	("Pu", ELE(ELE::Pu))	("Am", ELE(ELE::Am))	("Cm", ELE(ELE::Cm))	("Bk", ELE(ELE::Bk))	("Cf", ELE(ELE::Cf))
			("Es", ELE(ELE::Es))	("Fm", ELE(ELE::Fm))	("Md", ELE(ELE::Md))	("No", ELE(ELE::No))	("Lr", ELE(ELE::Lr))	("Rf", ELE(ELE::Rf));
	}
} element_p;

struct hydrogen_property : symbols<unsigned>
{
	hydrogen_property()
	{
		add
		("H1", 1)
		("H2", 2)
		("H3", 3)
		("H4", 4);
	}
} h_property_p;

struct charge_property : symbols<unsigned>
{
	charge_property()
	{
		add
		("+1", 1)
		("+2", 2)
		("+3", 3)
		("-1", -1)
		("-2", -2)
		("-3", -3);
	}
	
} ch_property_p;

struct bond_symbol : symbols<unsigned>
{
	bond_symbol()
	{
		add
		("-", 1)
		("=", 2)
		("#", 3)
		(":", 4)
		("~", 5);
	}
} bond_p;

struct atom_p
{
	atom_p(smart_at& ainfo) : _ainfo(ainfo) {}			
	void operator()(ELE c) const
	{
		atom_counter++;
		_ainfo.t=c;
		_ainfo.name=string(c);
	}
	smart_at& _ainfo;
};
struct charge_p
{
	charge_p(smart_at& ainfo) : _ainfo(ainfo) {}			
	void operator()(unsigned u) const
	{
		_ainfo.charge= u;
	}
	smart_at& _ainfo;
};
	
struct num_of_hydrogen_p
{
	num_of_hydrogen_p(smart_at& ainfo) : _ainfo(ainfo) {}			
	void operator()(unsigned u) const
	{
		_ainfo.num_of_hydrogen= u;
	}
	smart_at& _ainfo;
};
	
struct at_end_p
{
	at_end_p(smart_at& ainfo) : _ainfo(ainfo) {}			
	void operator()(char const* first, char const* last) const
	{
		_ainfo.at_id=atom_counter;
		myAtoms.push_back(_ainfo);
		_ainfo=smart_at();
	}
	smart_at& _ainfo;
};

struct at_negation_p
{
	at_negation_p(smart_at& ainfo) : _ainfo(ainfo) {}			
	void operator()(char const* first, char const* last) const
	{
		_ainfo.negation=1;
	}
	smart_at& _ainfo;
};

struct at_return_p
{
	at_return_p(smart_at& ainfo) : _ainfo(ainfo) {}			
	void operator()(char const* first, char const* last) const
	{
		_ainfo.return_sign=1;
	}
	smart_at& _ainfo;
};

struct ringbond_p
{
	void operator()(char const *first, char const *last) const
	{
		int bond(*find(bond_p, string(first,first+1).c_str()));
		int  number(atoi(string(first+1,last).c_str()));
		ringbonds.push_back(boost::tuple<unsigned, unsigned, unsigned> (atom_counter, bond, number));
	}
};

struct bt_smart_p
{
	void operator()(unsigned bond_type) const
	{
  		if (branch_info.size())
		{
		    if(branch_info.back().a)
			branch_info.pop_back();
		}
		for(vector< branch_content>::iterator i=branch_info.begin(); i!=branch_info.end();++i)
 			if(!(*i).a) 
				(*i).b++;
		bond.push_back(boost::tuple<unsigned, unsigned, unsigned> (atom_counter-test, bond_type, atom_counter+1));
		test=0;
	}
};

struct open_branch_action
{
	void operator()(char const *first, char const *last) const
	{
		if (branch_info.size())
		{
		    if(branch_info.back().a)
		    {
			unsigned int num = branch_info.back().b;
			branch_info.pop_back();
			branch_info.push_back(branch_content(0,num));
		    }
		    else
			branch_info.push_back(branch_content());
		}
		else
			branch_info.push_back(branch_content());
	}
};

struct close_branch_action
{
	void operator()(char const *first, char const *last) const
	{
		for(vector< branch_content>::reverse_iterator ri=branch_info.rbegin(); ri!=branch_info.rend();++ri)
			if(!(*ri).a) 
			{	
				(*ri).a++;
				break;
			}
		for(vector< branch_content>::iterator i=branch_info.begin(); i!=branch_info.end();++i)	
 			cout << (*i).b << " ";
 		cout << endl;
		test=branch_info.back().b;
	}
};

////////////////////////////////////////////////////////////////////////////
//
//  My SMARTS grammar
//
////////////////////////////////////////////////////////////////////////////
	
struct mySMARTS : public grammar<mySMARTS>
{
	template <typename ScannerT>
	struct definition
	{
		rule<ScannerT> any_atom, atom_symbol, atom_number, atom, open_square_bracket, close_square_bracket, negation, open_branch, close_branch, digit, atom_property, branch, molecule, my_operator, smarts_query, return_sign, ring_bond;
		definition(mySMARTS const& self)
		{
			open_square_bracket = ch_p('[')
					;
			close_square_bracket = ch_p(']')
					;
			/// additional signs
			negation = ch_p('!')
					;
			return_sign = ch_p('^')
					;
			/// atom definition
			atom = open_square_bracket >> 
					(	
						!(return_sign[at_return_p(self._p.ainfo)]) >> !(negation[at_negation_p(self._p.ainfo)]) >> element_p[atom_p(self._p.ainfo)]
						|
						!(negation[at_negation_p(self._p.ainfo)]) >> !(return_sign[at_return_p(self._p.ainfo)]) >> element_p[atom_p(self._p.ainfo)] 
					) >> 
					*atom_property >> 
					close_square_bracket[at_end_p(self._p.ainfo)]
					;
			/// numbers between 0 and 9
			uint_parser<unsigned, 10, 1, 1> uint1_p;
			digit = limit_d(0u, 9u)[uint1_p]
					;
			/// atom properties
			ring_bond = bond_p >> digit
				;
			atom_property = ch_property_p[charge_p(self._p.ainfo)] | 
					h_property_p[num_of_hydrogen_p(self._p.ainfo)]
					;
			/// branches
			open_branch = ch_p('(')
					;
			close_branch = ch_p(')')
					;
			branch = open_branch[open_branch_action()] >> 
				    bond_p[bt_smart_p()] >> 
				    atom >> 
						*ring_bond[ringbond_p()] >> 
						*(
						  	branch | 
						  	(
							 	bond_p[bt_smart_p()] >> 
								atom >> 
								*ring_bond[ringbond_p()]
							)
					 	) >> 
					close_branch[close_branch_action()]
					;
			/// whole molecule
			molecule = 
					atom >> 
					*ring_bond[ringbond_p()] >> 
					*( 
						branch | 
						(
						 	bond_p[bt_smart_p()] >> 
							atom >> 
							*ring_bond[ringbond_p()]
						)
					 )
					;
		}
	rule<ScannerT> const& start() const {return molecule; }			/// member function start returns a reference to the start rule
	};
	mySMARTS(parsed_info& p) : _p(p) {}
	parsed_info& _p;
};

// create a typedef for the SMARTS Graph type

struct smart_bt_t
{
	typedef edge_property_tag kind;
// 	unsigned int b_type;
};
typedef property<smart_bt_t, unsigned> EdgeProperty;
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, smart_at, EdgeProperty> MyGraph;
typedef boost::graph_traits<MyGraph>::vertex_descriptor MyVertex;

class MyVisitor : public boost::default_dfs_visitor
{
public:
  void discover_vertex(MyVertex v, const MyGraph& g) const
  {
    cout << v << "\tELE: "<< g[v].t << "\tID: "<< g[v].at_id << "\tcha: "<< g[v].charge << "\tnHy: "
    << g[v].num_of_hydrogen << "\tret: "<< g[v].return_sign << "\tneg: "<< g[v].negation << endl;

    return;
  }
};
/*
string name;
	ELE t;
	unsigned int at_id;
	int charge;		
	unsigned int num_of_hydrogen;
	bool negation;		
	bool return_sign;	*/

// ChainList
void mySmartsQuery(char const* str)
{
	parsed_info p(str);
	mySMARTS mysm(p);
	parse_info<> i=parse(str, mysm, space_p);
	if (i.full)
	{
		cout << "Finished with no errors." << endl;
		/// #################################################
		/// try to create Boost Graph Library form parsed data
		/// #################################################
		MyGraph g;
		/// insert vertices (with properties)
		for(unsigned i=0; i<myAtoms.size(); i++)
		{
		    MyGraph::vertex_descriptor v = add_vertex(g);
		    g[v] = myAtoms[i];
		}
		/// insert edges (with properties)
		for(unsigned i=0; i<bond.size(); i++)
		{
		    boost::add_edge(bond[i].get<0>(),bond[i].get<2>(), EdgeProperty(bond[i].get<1>()), g);		    
		}
		MyVisitor vis;
		boost::depth_first_search(g, boost::visitor(vis));
	
	   
// 	    typedef graph_traits<MyGraph>::edge_iterator edge_iter;
// 	    std::pair<edge_iter, edge_iter> ep;
// 	    edge_iter ei, ei_end;
// 	    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
// 		std::cout << edge_property_tag[*ei] << endl; 

	}
	else
	{
		cout << "parsing failed!" << endl;
	}
	cout << "number of atoms: " << myAtoms.size() << endl;
	cout << "###### general atom information ######" << endl;
	for(unsigned int i=0; i<myAtoms.size();i++)
		cout  <<  myAtoms[i].name << " " << myAtoms[i].t << " " << myAtoms[i].charge << " " << myAtoms[i].num_of_hydrogen << " " << myAtoms[i].at_id<< endl;
	cout << "####### ring bond information #######" << endl;
	for(unsigned int i=0; i<ringbonds.size();i++)
		cout  <<  ringbonds[i].get<0>() << " " <<  ringbonds[i].get<1>() << " " <<  ringbonds[i].get<2>() << " "<< endl;
	cout << "###### general bond information ######" << endl;
	for(unsigned int i=0; i<bond.size();i++)
		cout  <<  bond[i].get<0>() << " " <<  bond[i].get<1>() << " " <<  bond[i].get<2>() << " "<< endl;
	cout << "###### branch_content information ######" << endl;
	/*struct smart_at a = {"C",1,2,0,0};
	a_vertices.push_back(a);
	struct smart_at b;
	b.name="N";
	b.atkey=2;
	b.num_of_hydrogen=1;
	b.return_sign=true;
	b.negation=false;
	a_vertices.push_back(b);*/
	

// 	parse_info<> info = parse(str, ms, space_p);
/*	if(info.full) cout << "parsing succeed" << endl;
	else cout << "parsing failed" << endl;*/
}

