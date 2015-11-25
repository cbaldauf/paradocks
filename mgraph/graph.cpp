#include "graph.hpp"

//#include "combination.hpp"

namespace mgraph {

using namespace boost::lambda;

// global chain container
ChainList _chains;

//////////
// GSMF
////////

// chain functions
ChainKey cCreate(string n)
{
	// create a new chain and push it at the end of _chains
	ChainKey c(new Chain(n));
	_chains.push_back(c);
	#ifdef DEBUG
	// make valid
	c->_valid=true;
	#endif
	return c;
}

void cDelete(ChainKey c)
{
	// delete all residues
	rlDelete(cResidues(c));
	// remove chain from global chain container
	ChainList::iterator ki(find(_chains.begin(), _chains.end(), c));
	#ifdef DEBUG
	if (!c->Valid()) throw PhoenixError("cDelete: invalid ChainKey.");
	// invalidate chain
	c->_valid=false;
	#endif
	#ifdef DEBUG2
	if (ki==_chains.end()) throw PhoenixError("cDelete: can't delete non existing chain in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
	#endif
	_chains.erase(ki);
	#if !( defined(DEBUG) || defined(DEBUG2))
	// free memory
	delete c;
	#endif
}

void clDelete(const ChainList& cl)
{
	for_each(cl.begin(), cl.end(), bind(cDelete, _1));
}

// residue functions
ResidueKey rCreate(ChainKey p, RT rt)
{
	#ifdef DEBUG
	if (!p->Valid()) throw PhoenixError("rCreate: invalid ChainKey.");
	#endif
	// create a new residue
	ResidueKey r(new Residue(p, rt));
	// and push it at the end of the residue list of the parent chain
	p->_res.push_back(r);
	#ifdef DEBUG
	// make valid
	r->_valid=true;
	#endif
	return r;
}

void rDelete(ResidueKey r)
{
	// delete all atoms
	alDelete(rAtoms(r));
	// get parent chain
	ChainKey p=r->Parent();
	#ifdef DEBUG
	if (!p->Valid()) throw PhoenixError("rDelete: invalid ChainKey.");
	if (!r->Valid()) throw PhoenixError("rDelete: invalid ResidueKey.");
	// invalidate chain
	r->_valid=false;
	#endif
	// remove from parent chain
	ResidueList::iterator ki(find(p->_res.begin(), p->_res.end(), r));
	#ifdef DEBUG2
	if (ki==p->_res.end()) throw PhoenixError("rDelete: can't delete non existing residue in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
	#endif
	p->_res.erase(ki);
	#if !( defined(DEBUG) || defined(DEBUG2))
	// free memory
	delete r;
	#endif
}

void rlDelete(const ResidueList& rl)
{
	for_each(rl.begin(), rl.end(), bind(rDelete, _1));
}

// atom functions
AtomKey aCreate(ResidueKey p, ELE e, Vec3_t c, GEO g)
{
	#ifdef DEBUG
	if (!p->Valid()) throw PhoenixError("aCreate: invalid ResidueKey.");
	#endif
	// create a new atom
	AtomKey a(new Atom(p, e, c, g));
	// push it in residue p
	p->_atm.push_back(a);
	#ifdef DEBUG
	// make valid
	a->_valid=true;
	#endif
	return a;
}

void aDelete(AtomKey a)
{
	// delete all bonds of this atom
	blDelete(a->Bonds());
	// get parent residue
	ResidueKey p=a->Parent();
	#ifdef DEBUG
	if (!p->Valid()) throw PhoenixError("aDelete: invalid ResidueKey.");
	if (!a->Valid()) throw PhoenixError("aDelete: invalid AtomKey.");
	// invalidate atom
	a->_valid=false;
	#endif
	// remove from parent residue
	AtomList::iterator ki(find(p->_atm.begin(), p->_atm.end(), a));
	#ifdef DEBUG2
	if (ki==p->_atm.end()) throw PhoenixError("aDelete: can't delete non existing residue in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
	#endif
	p->_atm.erase(ki);
	#if !( defined(DEBUG) || defined(DEBUG2))
	// free memory
	delete a;
	#endif
}

void alDelete(const AtomList& al)
{
	for_each(al.begin(), al.end(), bind(aDelete, _1));
}

// bond functions
BondKey bCreate(AtomKey a1, AtomKey a2, BT bt)
{
	#ifdef DEBUG
	if (!a1->Valid() || !a2->Valid()) throw PhoenixError("bCreate: invalid AtomKey.");
	#endif
	// check for existence
	BondList bond(aaFormBond(a1, a2));
	if (!bond.empty())
	{
		#ifdef ERROR
		cout<<"bCreate: WARNING! Reasignment of BondType." << endl;
		#endif
		(*bond.begin())->Type(bt);
		return (*bond.begin());
	}
	// bond is new - create it
	BondKey b(new Bond(a1, a2, bt));
	// store the bond in the parent atoms
	a1->_bnd.push_back(b);
	a2->_bnd.push_back(b);
	#ifdef DEBUG
	// make valid
	b->_valid=true;
	#endif
	return b;
}

void bDelete(BondKey b)
{
	AtomKey a1=b->Atom1();
	AtomKey a2=b->Atom2();
	#ifdef DEBUG
	if (!a1->Valid() || !a2->Valid()) throw PhoenixError("bDelete: invalid AtomKey.");
	if (!b->Valid()) throw PhoenixError("bDelete: invalid BondKey.");
	b->_valid=false;
	#endif
	// remove from edge list
	BondList::iterator ki1(find(a1->_bnd.begin(), a1->_bnd.end(), b));
	BondList::iterator ki2(find(a2->_bnd.begin(), a2->_bnd.end(), b));
	#ifdef DEBUG2
	if (ki1==a1->_bnd.end() || ki2==a2->_bnd.end()) throw PhoenixError("bDelete: can't delete non existing bond in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
	#endif
	a1->_bnd.erase(ki1);
	a2->_bnd.erase(ki2);
	#if !( defined(DEBUG) || defined(DEBUG2))
	// free memory
	delete b;
	#endif
}

void blDelete(const BondList& bl)
{
	for_each(bl.begin(), bl.end(), bind(bDelete, _1));
}


//////////
// GSQF
////////

// chain functions
ChainList Chains()
{
	return _chains;
}

ResidueList cResidues(ChainKey c)
{
	return c->Residues();
}

ResidueList clResidues(const ChainList& cl)
{
	ResidueList result;
	for (ChainList::const_iterator i=cl.begin(); i!=cl.end(); ++i)
	{
		ResidueList tmp((*i)->Residues());
		copy(tmp.begin(), tmp.end(), back_inserter(result));
	}
	return result;
}

AtomList cAtoms(ChainKey c)
{
	return rlAtoms(c->Residues());
}

AtomList clAtoms(const ChainList& cl)
{
	return rlAtoms(clResidues(cl));
}

// residue functions
ResidueList Residues()
{
	return clResidues(_chains);
}

ChainKey rChain(ResidueKey r)
{
	return r->Parent();
}

ChainList rlChain(const ResidueList& rl)
{
	ChainList result;
	for (ResidueList::const_iterator i=rl.begin(); i!=rl.end(); ++i)
	{
		result.push_back((*i)->Parent());
	}
	return result;
}

AtomList rAtoms(ResidueKey r)
{
	return r->Atoms();
}

AtomList rlAtoms(ResidueList rl)
{
	AtomList result;
	for (ResidueList::iterator i=rl.begin(); i!=rl.end(); ++i)
	{
		AtomList tmp((*i)->Atoms());
		copy(tmp.begin(), tmp.end(), back_inserter(result));
	}
	return result;
}

// atom functions
AtomList Atoms()
{
	return clAtoms(_chains);
}

ChainKey aChain(AtomKey a)
{
	return a->Parent()->Parent();
}

ChainList alChains(const AtomList& al)
{
	ChainList result;
	for (AtomList::const_iterator i=al.begin(); i!=al.end(); ++i)
	{
		result.push_back((*i)->Parent()->Parent());
	}
	return result;
}

ResidueKey aResidue(AtomKey a)
{
	return a->Parent();
}

ResidueList alResidues(const AtomList& al)
{
	ResidueList result;
	for (AtomList::const_iterator i=al.begin(); i!=al.end(); ++i)
	{
		result.push_back((*i)->Parent());
	}
	return result;

}

BondList aBonds(AtomKey a)
{
	return a->Bonds();
}

BondList alBonds(const AtomList& al)
{
	set<BondKey> sorted_bonds;
	BondList result;
	for (AtomList::const_iterator i=al.begin(); i!=al.end(); ++i)
	{
		BondList tmp(aBonds(*i));
		for (BondList::iterator j=tmp.begin(); j!=tmp.end(); ++j)
		{
			if (sorted_bonds.find(*j)==sorted_bonds.end())
			{
				sorted_bonds.insert(*j);
				result.push_back(*j);
			}
		}
	}
	return result;
}

AtomList aAdjAtoms(AtomKey a)
{
	AtomList result;
	BondList bonds(aBonds(a));
	for (BondList::iterator i=bonds.begin(); i!=bonds.end(); ++i)
	{
		result.push_back( (*i)->Atom1()==a ? (*i)->Atom2() : (*i)->Atom1() );
	}
	return result;
}

BondList aaFormBond(AtomKey a1, AtomKey a2)
{
	BondList result;
	BondList bonds_a1(aBonds(a1));
	for (BondList::iterator i=bonds_a1.begin(); i!=bonds_a1.end(); ++i)
	{
		AtomKey a=((*i)->Atom1()==a1) ? (*i)->Atom2() : (*i)->Atom1();
		if (a==a2)
		{
			result.push_back(*i);
			return result;
		}
	}
	return result;
}

AtomList abAdjAtom(AtomKey a, BondKey b)
{
	AtomList result;
	BondList bonds(aBonds(a));
	for (BondList::iterator i=bonds.begin(); i!=bonds.end(); ++i)
	{
		if (b==(*i))
		{
			result.push_back(((*i)->Atom1()==a) ? (*i)->Atom2() : (*i)->Atom1());
			return result;
		}
	}
	return result;
}

// bond functions
BondList Bonds()
{
	return alBonds(Atoms());
}

AtomList bAtoms(BondKey b)
{
	AtomList result;
	result.push_back(b->Atom1());
	result.push_back(b->Atom2());
	return result;
}

AtomList blAtoms(const BondList& bl)
{
	// all atoms are put into a set - so the result is uniq
	set<AtomKey> sorted_atoms;
	AtomList result;
	for (BondList::const_iterator i=bl.begin(); i!=bl.end(); ++i)
	{
		if (sorted_atoms.find((*i)->Atom1())==sorted_atoms.end())
		{
			sorted_atoms.insert((*i)->Atom1());
			result.push_back((*i)->Atom1());
		}
		if (sorted_atoms.find((*i)->Atom2())==sorted_atoms.end())
		{
			sorted_atoms.insert((*i)->Atom2());
			result.push_back((*i)->Atom2());
		}
	}
	return result;
}

// ring check
bool isRingBond(BondKey b)
{
	map<AtomKey, bool> visited;
	// bfs search by using a queue
	queue<AtomKey> disc;
	// atom1 is visited
	visited[b->Atom1()]=1;
	AtomKey a2=b->Atom2();
	
	// start the bfs at all adjacent atoms of a1, exept a2
	AtomList a1_adj(aAdjAtoms(b->Atom1()));
	for (AtomList::iterator i=a1_adj.begin(); i!=a1_adj.end(); ++i)
	{
		if ((*i)!=a2)
		{
			disc.push(*i);
			visited[(*i)]=1;
		}
	}
	// bfs search
	while (!disc.empty())
	{
		AtomKey front=disc.front();
		disc.pop();
		AtomList front_adj(aAdjAtoms(front));
		for (AtomList::iterator i=front_adj.begin(); i!=front_adj.end(); ++i)
		{
			if (!visited[*i])
			{
				if ((*i)==a2) return true;
				disc.push(*i);
				visited[(*i)]=1;
			}
		}
	}
	return false;
}

bool isXRingBond(BondKey b, unsigned int x)
{
	// at least 3 members to form a ring
	if (x<3) return false;
	
	AtomKey a2=b->Atom2();
	AtomKey a1=b->Atom1();
	// start the bfs at all adjacent atoms of a1, exept a2
	AtomList a1_adj(aAdjAtoms(a1));
	for (AtomList::iterator i=a1_adj.begin(); i!=a1_adj.end(); ++i)
	{
		if ((*i)==a2) continue;
		map<AtomKey, bool> visited;
		// bfs search by using a queue
		// maintain a counter for bfs depth
		queue<pair<AtomKey, unsigned int> > disc;
		// atom1 is visited
		visited[a1]=1;
		
		// 2 bonds are allready discovered
		disc.push(pair<AtomKey, unsigned int>(*i,2));
		visited[(*i)]=1;
				
		// bfs search
		while (!disc.empty())
		{
			pair<AtomKey, unsigned int> front=disc.front();
			disc.pop();
			if (front.second>x) continue;
			AtomList front_adj(aAdjAtoms(front.first));
			for (AtomList::iterator i=front_adj.begin(); i!=front_adj.end(); ++i)
			{
				if (!visited[*i])
				{
					// if this bond goes to a2 and there are x-1 bonds allready
					// discovered it is a x boned ring
					if ((*i)==a2 && front.second==x-1)
					{
						return true;
					}
					disc.push(pair<AtomKey, unsigned int>(*i,front.second+1));
					visited[(*i)]=1;
				}
			}
		}
	}
	return false;
}

bool isRingAtom(AtomKey a)
{
	map<AtomKey, bool> visited;
	// bfs search by using a queue
	queue<AtomKey> disc;
	
	// start the search at at all adjacent atoms of the adjacent atoms of a
	AtomList a_adj(aAdjAtoms(a));
	// at least two neighbour to form a ring
	if (a_adj.size()<2) return false;
	for (AtomList::iterator i=a_adj.begin(); i!=(a_adj.end()--); ++i)
	{
		visited[*i]=1;
		// initialize disc
		AtomList aa_adj(aAdjAtoms(*i));
		for (AtomList::iterator j=aa_adj.begin(); j!=aa_adj.end(); ++j)
		{
			
			if ((*j)!=a)
			{
				disc.push(*j);
				visited[*j]=1;
			}
		}
		while (!disc.empty())
		{
			AtomKey front=disc.front();
			disc.pop();
			AtomList front_adj(aAdjAtoms(front));
			for (AtomList::iterator k=front_adj.begin(); k!=front_adj.end(); ++k)
			{
				if (!visited[*k])
				{
					if ((*k)==a) return true;
					disc.push(*k);
					visited[*k]=1;
				}
			}
		}
	}
	return false;
}

bool isXRingAtom(AtomKey a, unsigned int x)
{
	// at least 3 members to form a ring
	if (x<3) return false;
	
	// start the search at at all adjacent atoms of the adjacent atoms of a
	AtomList a_adj(aAdjAtoms(a));
	// at least two neighbour to form a ring
	if (a_adj.size()<2) return false;
	for (AtomList::iterator i=a_adj.begin(); i!=(a_adj.end()--); ++i)
	{
		map<AtomKey, bool> visited;
		// bfs search by using a queue
		queue<pair<AtomKey, unsigned int> > disc;
		visited[*i]=1;
		// initialize disc
		AtomList aa_adj(aAdjAtoms(*i));
		for (AtomList::iterator j=aa_adj.begin(); j!=aa_adj.end(); ++j)
		{
			
			if ((*j)!=a)
			{
				disc.push(pair<AtomKey, unsigned int>(*j,2));
				visited[*j]=1;
			}
		}
		while (!disc.empty())
		{
			pair<AtomKey, unsigned int> front=disc.front();
			disc.pop();
			if (front.second>x) continue;
			AtomList front_adj(aAdjAtoms(front.first));
			for (AtomList::iterator k=front_adj.begin(); k!=front_adj.end(); ++k)
			{
				if (!visited[*k])
				{
					if ((*k)==a && front.second==x-1) return true;
					disc.push(pair<AtomKey, unsigned int>(*k,front.second+1));
					visited[*k]=1;
				}
			}
		}
	}
	return false;
}

AtomList aExtendConnected(AtomKey a)
{
	// simple bfs to find all connected atoms
	map<AtomKey, bool> visited;
	visited[a]=1;
	queue<AtomKey> disc;
	disc.push(a);
	AtomList result;
	result.push_back(a);
	
	while (!disc.empty())
	{
		AtomKey front=disc.front();
		disc.pop();
		AtomList front_adj(aAdjAtoms(front));
		for (AtomList::iterator k=front_adj.begin(); k!=front_adj.end(); ++k)
		{
			if (!visited[*k])
			{
				result.push_back(*k);
				disc.push(*k);
				visited[*k]=1;
			}
		}
	}
	return result;
}

AtomList alExtendConnected(AtomList al)
{
	// simple bfs to find all connected atoms
	map<AtomKey, bool> visited;
	queue<AtomKey> disc;
	AtomList result;
	for (AtomList::iterator i=al.begin(); i!=al.end(); i++)
	{
		visited[*i]=1;
		disc.push(*i);
		result.push_back(*i);
	}
	
	while (!disc.empty())
	{
		AtomKey front=disc.front();
		disc.pop();
		AtomList front_adj(aAdjAtoms(front));
		for (AtomList::iterator k=front_adj.begin(); k!=front_adj.end(); ++k)
		{
			if (!visited[*k])
			{
				result.push_back(*k);
				disc.push(*k);
				visited[*k]=1;
			}
		}
	}
	return result;
}

//////////
// CKQF
////////
unsigned int aHeavyCount(AtomKey a)
{
	AtomList adj_a(aAdjAtoms(a));
	unsigned int heavy=0;
	for (AtomList::iterator i=adj_a.begin(); i!=adj_a.end(); ++i)
	{
		if ((*i)->Element()>ELE::H) ++heavy;
	}
	return heavy;
}

unsigned int aHCount(AtomKey a)
{
	AtomList adj_a(aAdjAtoms(a));
	unsigned int h=0;
	for (AtomList::iterator i=adj_a.begin(); i!=adj_a.end(); ++i)
	{
		if ((*i)->Element()==ELE::H) ++h;
	}
	return h;
}

unsigned int aLPCount(AtomKey a)
{
	AtomList adj_a(aAdjAtoms(a));
	unsigned int lp=0;
	for (AtomList::iterator i=adj_a.begin(); i!=adj_a.end(); ++i)
	{
		if ((*i)->Element()==ELE::LP) ++lp;
	}
	return lp;
}

bool aPolar(AtomKey a)
{
	AtomList adj_a(aAdjAtoms(a));
	for (AtomList::iterator j=adj_a.begin(); j!=adj_a.end(); ++j)
	{
		if (((*j)->Element()>ELE::H) && ((*j)->Element()!=ELE::C)) return true;
	}
	return false;
}

bool aMetall(AtomKey a)
{
	if (a->Element()==ELE::Li ||
		a->Element()==ELE::Na ||
		a->Element()==ELE::Mg ||
		a->Element()==ELE::Al ||
		a->Element()==ELE::K  ||
		a->Element()==ELE::Ca ||
		a->Element()==ELE::Mn ||
		a->Element()==ELE::Fe ||
		a->Element()==ELE::Co ||
		a->Element()==ELE::Ni ||
		a->Element()==ELE::Cu ||
		a->Element()==ELE::Zn ||
		a->Element()==ELE::Pb ||
		a->Element()==ELE::U) return true;
	else return false;
}

/*vector<Vec3_t> aLPPos(AtomKey a)
{
	AtomList a_adj(aAdjAtoms(a));
	vector<Vec3_t> result;
	// tetraeder
	if (a->Geometry()==GEO::tet)
	{
		// add one lp
		if (a_adj.size()==3)
		{
			// let A, B, C be the vectors from atom a to the atoms in a_adj
			// we want to have angles A-LP, B-LP, C-LP be equal
			// A*LP+B*LP = 2* C*LP
			// A*LP+B*LP-2*C*LP=0
			// (A+B-2C)*LP=0
			// (A*B-2C) is orthogonal to LP
			// (A*C-2B) is orthogonal to LP
			// => (A*B-2C)*(A*C-2B) ist linear to LP
			vector<Vec3_t> c;
			for (AtomList::iterator i=a_adj.begin(); i!=a_adj.end(); ++i)
			{
				Vec3_t tmp=(*i)->Coord()-a->Coord();
				tmp.normalize();
				c.push_back(tmp);
			}
			Vec3_t tmp=(c[0]+c[1]-(c[2]*2))^(c[0]+c[2]-(c[1]*2));
			tmp.normalize();
			if (acos(c[0]*tmp)>acos(c[0]*(tmp*-1)))
			{
				tmp+=a->Coord();
				result.push_back(tmp);
			}
			else
			{
				tmp=a->Coord()-tmp;
				result.push_back(tmp);
			}
		}
		// add two lp
		if (a_adj.size()==2)
		{
			// mathematical not correct, but its close
			Vec3_t v1, v2;
			v1=a_adj[0]->Coord()-a->Coord();
			v2=a_adj[1]->Coord()-a->Coord();
			v1.normalize(); v2.normalize();
			// normal of lp plane
			Vec3_t v3=(v1-v2);
			// vector between lp
			Vec3_t v4=(v1+v2)*-1;
			Quat_t rot(osg::DegreesToRadians(109.5)/2, v3);
			Vec3_t tmp=rot*v4;
			tmp+=a->Coord();
			result.push_back(tmp);
			tmp=(rot.inverse())*v4;
			tmp+=a->Coord();
			result.push_back(tmp);
		}
		// add three lp
		if (a_adj.size()==1)
		{
			// we assign a perfect template
			// root is (1,0,0)
			Vec3_t root=a_adj[0]->Coord()-a->Coord();
			Vec3_t v1(-0.333807, 0, -0.942641);
			Vec3_t v2(-0.333807, 0.816351, 0.471321);
			Vec3_t v3(-0.333807, -0.816351, 0.471321);
			Quat_t rot;
			rot.makeRotate (Vec3_t(1,0,0), root);
			v1=rot*v1+a->Coord();
			v2=rot*v2+a->Coord();
			v3=rot*v3+a->Coord();
			result.push_back(v1);
			result.push_back(v2);
			result.push_back(v3);
		}
	}
	if (a->Geometry()==GEO::tri)
	{
		// add one lp
		if (a_adj.size()==2)
		{
			Vec3_t v1, v2;
			v1=a_adj[0]->Coord()-a->Coord();
			v2=a_adj[1]->Coord()-a->Coord();
			v1.normalize(); v2.normalize();
			Vec3_t tmp=(v1-v2)^(v1^v2);
			tmp.normalize();
			if (acos(v1*tmp)>acos(v1*(tmp*-1)))
			{
				tmp+=a->Coord();
				result.push_back(tmp);
			}
			else
			{
				tmp=a->Coord()-tmp;
				result.push_back(tmp);
			}
		}
		// add two lp
		if (a_adj.size()==1)
		{
			Vec3_t root=a_adj[0]->Coord()-a->Coord();
			root.normalize();
			AtomList adj_a_adj(aAdjAtoms(a_adj[0]));
			Vec3_t v1(1,0,0);
			for (AtomList::iterator i=adj_a_adj.begin(); i!=adj_a_adj.end(); ++i)
			{
				if ((*i)!=a)
				{
					v1=(*i)->Coord()-a_adj[0]->Coord();
					v1.normalize();
					break;
				}
			}
			Vec3_t v2=root^v1;
			Quat_t rot(osg::DegreesToRadians(120.0), v2);
			Vec3_t tmp=rot*root;
			tmp+=a->Coord();
			result.push_back(tmp);
			tmp=(rot.inverse())*root;
			tmp+=a->Coord();
			result.push_back(tmp);
		}
	}
	if (a->Geometry()==GEO::lin)
	{
		// add one lp
		if (a_adj.size()==1)
		{
			Vec3_t root=a_adj[0]->Coord()-a->Coord();
			root.normalize();
			Vec3_t tmp=root*-1+a->Coord();
			result.push_back(tmp);
		}
	}
	return result;
}*/



void CorrectGraphByTopology(ChainList cl)
{
	// delete lonepairs
	AtomList a(clAtoms(cl));
	for (mgraph::AtomList::iterator i=a.begin(); i!=a.end(); ++i)
	{
		if ((*i)->Element()==ELE::LP)
		{
			aDelete(*i);
		}
	}
	
	// no bonds on metall ions
	a=clAtoms(cl);
	for (mgraph::AtomList::iterator i=a.begin(); i!=a.end(); ++i)
	{
		if (aMetall(*i)) blDelete(aBonds(*i));
	}

	a=clAtoms(cl);
	BondList b(alBonds(a));
	// untested atoms and bonds 
	set<AtomKey> ua(a.begin(), a.end());
	set<BondKey> ub(b.begin(), b.end());
	#ifdef DEBUG
	cout<<ua.size()<<" atoms and "<<ub.size()<<" bonds to check..."<<endl;
	#endif
	
	set<AtomKey>::iterator i=ua.begin();
	while (i!=ua.end())
	{
		
		BondList bonds=aBonds(*i);
		AtomList a_adj(aAdjAtoms(*i));
		AtomList heavy;
		AtomList hydrogen;
		BondList heavy_b;
		BondList hydrogen_b;
		
		for (AtomList::iterator j=a_adj.begin(); j!=a_adj.end(); ++j)
		{
			if ((*j)->Element()!=ELE::H)
			{
				heavy.push_back(*j);
				BondList b=aaFormBond(*i, *j);
				#ifdef DEBUG
				if (b.empty()) throw PhoenixError("CorrectGraphByTopology: inconsistency in graph structure in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
				#endif
				heavy_b.push_back(*(b.begin()));
			}
			if ((*j)->Element()==ELE::H)
			{
				hydrogen.push_back(*j);
				BondList b=aaFormBond(*i, *j);
				#ifdef DEBUG
				if (b.empty()) throw PhoenixError("CorrectGraphByTopology: inconsistency in graph structure in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
				#endif
				hydrogen_b.push_back(*(b.begin()));
			}
		}
		
		// hydrogens
		if ((*i)->Element()==ELE::H)
		{
			if (bonds.size()>1) throw PhoenixError("CorrectGraphByTopology: H with # of bonds > 1.");
			(*i)->Geometry(GEO::none);
			if (bonds.size()==1)
			{
				(*i)->FormalCharge(0);
				if (ub.find(*bonds.begin())!=ub.end())
				{
					if (ub.find(*bonds.begin())!=ub.end())
					{
						(*bonds.begin())->Type(BT::s);
						ub.erase(*bonds.begin());
					}
				}
				ua.erase(i++);
				continue;
			}
			if (bonds.size()==0)
			{
				(*i)->FormalCharge(1);
				ua.erase(i++);
				continue;
			}
		}
		
		// atoms with no bonds are normaly ions
		// eg. halogene ions and metall ions
		if (bonds.empty())
		{
			// +1 ions
			if ((*i)->Element()==ELE::Li ||
			    (*i)->Element()==ELE::Na ||
			    (*i)->Element()==ELE::K)
			{
				(*i)->Geometry(GEO::none);
				(*i)->FormalCharge(1);
				ua.erase(i++);
				continue;
			}
			// +2 ions
			if ((*i)->Element()==ELE::Mg ||
			    (*i)->Element()==ELE::Ca ||
			    (*i)->Element()==ELE::Zn ||
			    (*i)->Element()==ELE::Fe ||
			    (*i)->Element()==ELE::Pb ||
			    (*i)->Element()==ELE::U)
			{
				(*i)->Geometry(GEO::none);
				(*i)->FormalCharge(2);
				ua.erase(i++);
				continue;
			}
			if ((*i)->Element()==ELE::Al)
			{
				(*i)->Geometry(GEO::none);
				(*i)->FormalCharge(3);
				ua.erase(i++);
				continue;
			}
			// -1 ions
			if ((*i)->Element()==ELE::F ||
			    (*i)->Element()==ELE::Cl ||
			    (*i)->Element()==ELE::Br ||
			    (*i)->Element()==ELE::I)
			{
				(*i)->Geometry(GEO::none);
				(*i)->FormalCharge(-1);
				ua.erase(i++);
				continue;
			}
		}
		
		// halogenes with only one bond
		if (bonds.size()==1)
		{
			if ((*i)->Element()==ELE::F  ||
			    (*i)->Element()==ELE::Cl ||
			    (*i)->Element()==ELE::Br ||
			    (*i)->Element()==ELE::I)
			{
				(*i)->Geometry(GEO::tet);
				(*i)->FormalCharge(0);
				
				BondList b;
				b=aaFormBond(*i, *a_adj.begin());
				#ifdef DEBUG
				if (b.empty()) throw PhoenixError("CorrectGraphByTopology: inconsistency in graph structure in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
				#endif
				if (ub.find(*b.begin())!=ub.end())
				{
					(*bonds.begin())->Type(BT::s);
					ub.erase(*bonds.begin());
				}
				ua.erase(i++);
				continue;
			}
		}
		
		// atoms with 4 bonds are sp3 carbon and positively charged nitrogens
		if (bonds.size()==4)
		{
			if ((*i)->Element()==ELE::C)
			{
				(*i)->Geometry(GEO::tet);
				(*i)->FormalCharge(0);
				for (BondList::iterator j=bonds.begin(); j!=bonds.end(); ++j)
				{
					if (ub.find(*j)!=ub.end())
					{
						(*j)->Type(BT::s);
						ub.erase(*j);
					}
				}
				ua.erase(i++);
				continue;
			}
			if ((*i)->Element()==ELE(ELE::N))
			{
				(*i)->Geometry(GEO::tet);
				(*i)->FormalCharge(1);
				for (BondList::iterator j=bonds.begin(); j!=bonds.end(); ++j)
				{
					if (ub.find(*j)!=ub.end())
					{
						(*j)->Type(BT::s);
						ub.erase(*j);
					}
				}
				ua.erase(i++);
				continue;
			}
		}
		
		// easy functional groups
		// hydroxy/thiol ether/thioether group
		if (((*i)->Element()==ELE::O || (*i)->Element()==ELE::S) && bonds.size()==2)
		{
			if (!isRingAtom(*i))
			{
				AtomList hc;
				for (AtomList::iterator j=a_adj.begin(); j!=a_adj.end(); ++j)
				{
					if ((*j)->Element()==ELE::C || (*j)->Element()==ELE::H)
					{
						hc.push_back(*j);
					}
				}
				if (hc.size()==2)
				{
					(*i)->Geometry(GEO::tet);
					(*i)->FormalCharge(0);
					for (BondList::iterator j=bonds.begin(); j!=bonds.end(); ++j)
					{
						// assign bonds
						if (ub.find(*j)!=ub.end())
						{
							(*j)->Type(BT::s);
							ub.erase(*j);
						}
					}
					ua.erase(i++);
					continue;
				}
			}
		}
		
		// carboxylic acid
		if ((*i)->Element()==ELE::C && bonds.size()==3)
		{
			AtomList oxygens;
			BondList re, s;
			for (AtomList::iterator j=a_adj.begin(); j!=a_adj.end(); ++j)
			{
				if ((*j)->Element()==ELE::O && (aBonds(*j).size()==1))
				{
					oxygens.push_back(*j);
					BondList b;
					b=aaFormBond(*i, *j);
					#ifdef DEBUG
					if (b.empty()) throw PhoenixError("CorrectGraphByTopology: inconsistency in graph structure in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
					#endif
					re.push_back(*b.begin());
				}
				else if ((*j)->Element()==ELE::C)
				{
					BondList b;
					b=aaFormBond(*i, *j);
					#ifdef DEBUG
					if (b.empty()) throw PhoenixError("CorrectGraphByTopology: inconsistency in graph structure in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
					#endif
					s.push_back(*b.begin());
				}
			}
			if (oxygens.size()==2 && s.size()==1)
			{
				// assign central carbon
				(*i)->Geometry(GEO::tri);
				(*i)->FormalCharge(0);
				// assign oxigens
				for (AtomList::iterator j=oxygens.begin(); j!=oxygens.end(); ++j)
				{
					if (ua.find(*j)!=ua.end())
					{
						(*j)->Geometry(GEO::tri);
						(*j)->FormalCharge(-0.5);
						ua.erase(*j);
					}
				}
				// assign bonds to the oxygens
				for (BondList::iterator j=re.begin(); j!=re.end(); ++j)
				{
					if (ub.find(*j)!=ub.end())
					{
						(*j)->Type(BT::re);
						ub.erase(*j);
					}
				}
				// assign root bond
				for (BondList::iterator j=s.begin(); j!=s.end(); ++j)
				{
					if (ub.find(*j)!=ub.end())
					{
						(*j)->Type(BT::s);
						ub.erase(*j);
					}
				}
				ua.erase(i++);
				continue;
			}
		}
		
		// nitro group
		if ((*i)->Element()==ELE::N && bonds.size()==3)
		{
			AtomList oxygens;
			BondList re, s;
			for (AtomList::iterator j=a_adj.begin(); j!=a_adj.end(); ++j)
			{
				if ((*j)->Element()==ELE::O && (aBonds(*j).size()==1))
				{
					oxygens.push_back(*j);
					BondList b;
					b=aaFormBond(*i, *j);
					#ifdef DEBUG
					if (b.empty()) throw PhoenixError("CorrectGraphByTopology: inconsistency in graph structure in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
					#endif
					re.push_back(*b.begin());
				}
				else if ((*j)->Element()==ELE::C)
				{
					BondList b;
					b=aaFormBond(*i, *j);
					#ifdef DEBUG
					if (b.empty()) throw PhoenixError("CorrectGraphByTopology: inconsistency in graph structure in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
					#endif
					s.push_back(*b.begin());
				}
			}
			if (oxygens.size()==2 && s.size()==1)
			{
				// assign central nitrogen
				(*i)->Geometry(GEO::tri);
				(*i)->FormalCharge(1);
				// assign oxigens
				for (AtomList::iterator j=oxygens.begin(); j!=oxygens.end(); ++j)
				{
					if (ua.find(*j)!=ua.end())
					{
						(*j)->Geometry(GEO::tri);
						(*j)->FormalCharge(-0.5);
						ua.erase(*j);
					}
				}
				// assign bonds to the oxygens
				for (BondList::iterator j=re.begin(); j!=re.end(); ++j)
				{
					if (ub.find(*j)!=ub.end())
					{
						(*j)->Type(BT::re);
						ub.erase(*j);
					}
				}
				// assign root bond
				for (BondList::iterator j=s.begin(); j!=s.end(); ++j)
				{
					if (ub.find(*j)!=ub.end())
					{
						(*j)->Type(BT::s);
						ub.erase(*j);
					}
				}
				ua.erase(i++);
				continue;
			}
		}
		
		// amide group and peptide bonds
		if ((*i)->Element()==ELE::C && bonds.size()==3 && !isXRingAtom(*i, 5) && !isXRingAtom(*i, 6))
		{
			AtomList oxygen, nitrogen, carbon;
			for (AtomList::iterator j=a_adj.begin(); j!=a_adj.end(); ++j)
			{
				if ((*j)->Element()==ELE::O && (aBonds(*j).size()==1))
				{
					oxygen.push_back(*j);
				}
				else if ((*j)->Element()==ELE::N && (aBonds(*j).size()==3))
				{
					AtomList bn=aAdjAtoms(*j);
					unsigned int hccount=0;
					for (AtomList::iterator k=bn.begin(); k!=bn.end(); ++k)
					{
						if ((*k)->Element()==ELE::C || (*k)->Element()==ELE::H)
						{
							++hccount;
						}
					}
					if (hccount==3) nitrogen.push_back(*j);
				}
				else if ((*j)->Element()==ELE::C || (*j)->Element()==ELE::H)
				{
					carbon.push_back(*j);
				}
			}
			if (oxygen.size()==1 && nitrogen.size()==1 && carbon.size()==1)
			{
				// assign central carbon
				(*i)->Geometry(GEO::tri);
				(*i)->FormalCharge(0);
				// assign oxygen
				if (ua.find(*oxygen.begin())!=ua.end())
				{
					(*oxygen.begin())->Geometry(GEO::tri);
					(*oxygen.begin())->FormalCharge(-0.5);
					ua.erase(*oxygen.begin());
				}
				// assign the bond to the oxygen
				BondList bl=aaFormBond(*i, *oxygen.begin());
				#ifdef DEBUG
				if (bl.empty()) throw PhoenixError("CorrectGraphByTopology: inconsistency in graph structure in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
				#endif
				BondKey b=*bl.begin();
				if (ub.find(b)!=ub.end())
				{
					b->Type(BT::re);
					ub.erase(b);
				}
				// assign nitrogen
				if (ua.find(*nitrogen.begin())!=ua.end())
				{
					(*nitrogen.begin())->Geometry(GEO::tri);
					(*nitrogen.begin())->FormalCharge(0.5);
					ua.erase(*nitrogen.begin());
				}
				// assign the bond to the nitrogen
				bl=aaFormBond(*i, *nitrogen.begin());
				#ifdef DEBUG
				if (bl.empty()) throw PhoenixError("CorrectGraphByTopology: inconsistency in graph structure in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
				#endif
				b=*bl.begin();
				if (ub.find(b)!=ub.end())
				{
					b->Type(BT::re);
					ub.erase(b);
				}
				// remaining bonds on the nitrogen are single bonds
				BondList nbnd=aBonds(*nitrogen.begin());
				for (BondList::iterator j=nbnd.begin(); j!=nbnd.end(); ++j)
				{
					if (ub.find(*j)!=ub.end())
					{
						(*j)->Type(BT::s);
						ub.erase(*j);
					}
				}
				// assign the bond to the carbon
				bl=aaFormBond(*i, *carbon.begin());
				#ifdef DEBUG
				if (bl.empty()) throw PhoenixError("CorrectGraphByTopology: inconsistency in graph structure in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
				#endif
				b=*bl.begin();
				if (ub.find(b)!=ub.end())
				{
					b->Type(BT::s);
					ub.erase(b);
				}
				ua.erase(i++);
				continue;
			}
		}
		
		// hier weiter machen
		
		if (bonds.size()==1)
		{
			if ((*i)->Element()==ELE::N)
			{
				(*i)->Geometry(GEO::lin);
				(*i)->FormalCharge(0);
				//assign the bond
				if (ub.find(*bonds.begin())!=ub.end())
				{
					(*bonds.begin())->Type(BT::t);
					ub.erase(*bonds.begin());
				}
				ua.erase(i++);
				continue;
			}
			// isocyanide
			if ((*i)->Element()==ELE::C)
			{
				if ((*a_adj.begin())->Element()==ELE::N)
				{
					(*i)->Geometry(GEO::lin);
					(*i)->FormalCharge(-1);
					// bond to the nitrogen
					if (ub.find(*bonds.begin())!=ub.end())
					{
						(*bonds.begin())->Type(BT::t);
						ub.erase(*bonds.begin());
					}
					
					BondList n_bonds(aBonds(*a_adj.begin()));
					if (ua.find(*a_adj.begin())!=ua.end())
					{
						
						(*a_adj.begin())->Geometry(GEO::lin);
						if (n_bonds.size()==2) (*a_adj.begin())->FormalCharge(+1);
						if (n_bonds.size()==1) (*a_adj.begin())->FormalCharge(0);
						ua.erase(*a_adj.begin());
						// remaining bond in ub is the single bond on N
						for (BondList::iterator j=n_bonds.begin(); j!=n_bonds.end(); ++j)
						{
							if (ub.find(*j)!=ub.end())
							{
								(*j)->Type(BT::s);
								ub.erase(*j);
							}	
						}
					}
				}
				ua.erase(i++);
				continue;
			}
		}
		
		// search guanidino group
		if ((*i)->Element()==ELE::C)
		{
			if (a_adj.size()==3)
			{
				if (!isRingAtom(*i))
				{
					AtomList nitrogen;
					for (AtomList::iterator j=a_adj.begin(); j!=a_adj.end(); ++j)
					{
						if ((*j)->Element()==ELE::N &&
						    (aBonds(*j).size()==3) &&
						    !isRingAtom(*j))
						{
							nitrogen.push_back(*j);
						}
					}
					if (nitrogen.size()==3)
					{
						// assign central carbon
						(*i)->Geometry(GEO::tri);
						(*i)->FormalCharge(0);
						for (AtomList::iterator j=nitrogen.begin(); j!=nitrogen.end(); ++j)
						{
							if (ua.find(*j)!=ua.end())
							{
								(*j)->Geometry(GEO(GEO::tri));
								(*j)->FormalCharge(+0.333);
								ua.erase(*j);
							}
							BondKey cnbond=*(aaFormBond((*i), (*j))).begin();
							if (ub.find(cnbond)!=ub.end())
							{
								cnbond->Type(BT(BT::re));
								ub.erase(cnbond);
							}
							// remaining bonds on the nitrogen are single bonds
							BondList nbnd=aBonds(*j);
							for (BondList::iterator k=nbnd.begin(); k!=nbnd.end(); ++k)
							{
								if (ub.find(*k)!=ub.end())
								{
									(*k)->Type(BT::s);
									ub.erase(*k);
								}
							}
						}
						ua.erase(i++);
						continue;
					}
				}
			}
		}
		
		
		
		
		
		// remove aromatic atoms
		if (bonds.size()==3)
		{
			BondList ar, s;
			for (unsigned int j=0 ; j!=bonds.size(); ++j)
			{
				if (bonds[j]->Type()==BT::ar)
				{
					ar.push_back(bonds[j]);
				}
				if (bonds[j]->Type()==BT::s)
				{
					s.push_back(bonds[j]);
				}
			}
			if ((ar.size()==2 && s.size()==1) || ar.size()==3)
			{
				(*i)->Geometry(GEO::tri);
				(*i)->FormalCharge(0);
				for (BondList::iterator j=bonds.begin(); j!=bonds.end(); ++j)
				{
					if (ub.find(*j)!=ub.end())
					{
						ub.erase(*j);
					}
				}
				ua.erase(i++);
				continue;
			}
		}
		// remove sp2 atoms
		if (bonds.size()==3)
		{
			BondList d, s;
			for (unsigned int j=0 ; j!=bonds.size(); ++j)
			{
				if (bonds[j]->Type()==BT::d)
				{
					d.push_back(bonds[j]);
				}
				if (bonds[j]->Type()==BT::s)
				{
					s.push_back(bonds[j]);
				}
			}
			if (d.size()==1 && s.size()==2)
			{
				(*i)->Geometry(GEO::tri);
				(*i)->FormalCharge(0);
				for (BondList::iterator j=bonds.begin(); j!=bonds.end(); ++j)
				{
					if (ub.find(*j)!=ub.end())
					{
						ub.erase(*j);
					}
				}
				ua.erase(i++);
				continue;
			}
		}
		++i;
	}
		
//	#ifdef DEBUG
//	cout<<ua.size()<<" atoms and "<<ub.size()<<" bonds not succesfully checked!"<<endl;
//	cout<<"this is not necessarily an error."<<endl;
//	i=ua.begin();
/*	while (i!=ua.end())
	{
		if ((*i)->Parent()->Type()!= RT::HIS && (*i)->Parent()->Type()!= RT::TRP )
		{
			cout<<"atom: "<<(*i)->Name()<<" res: "<<(*i)->Parent()->Name()<< " geo: "<<(*i)->Geometry() <<endl ;
		}
		++i;
	}*/
//	#endif	
}

void MOL::extractMOL(const ChainList& cl)
{
	AtomList atoms(clAtoms(cl));
	AtomList atoms_ext(alExtendConnected(atoms));
	#ifdef ERROR
	if (atoms.size()!=atoms_ext.size()) throw 
	PhoenixError("MOL::extractMOL: at least one of the atoms belonging to cl\nare connected to another chain.");
	#endif
	
	unsigned int resid=0;
	unsigned int atomid=0;
	_c.clear();
	_cres.clear();
	_r.clear();
	_ratoms.clear();
	_a.clear();
	
	map<AtomKey, unsigned int> key2index;
	// for each chain
	for (ChainList::const_iterator i=cl.begin(); i!=cl.end(); ++i)
	{
		// get chain information
		CHAIN tmpc((*i)->Name(), (*i)->Tag());
		// save chain information
		_c.push_back(tmpc);
		// for each residue of the chain
		ResidueList r(cResidues(*i));
		vector<unsigned int> cres;
		for (ResidueList::iterator j=r.begin(); j!=r.end(); ++j)
		{
			// get residue information
			RESIDUE tmpr((*j)->Type(), (*j)->Name(), (*j)->Id());
			// save residue information
			_r.push_back(tmpr);
			// save id of the residue; for mapping of residue to chain
			cres.push_back(resid);
			// for each atom of the residue
			AtomList a(rAtoms(*j));
			vector<unsigned int> ratoms;
			for (AtomList::iterator k=a.begin(); k!=a.end(); ++k)
			{
				// get atom information
				ATOM tmpa((*k)->Element(), (*k)->Geometry(), (*k)->Coord(), (*k)->Name(),
					(*k)->FFType(), (*k)->Charge(), (*k)->FormalCharge(), (*k)->MOL2Type());
				// save atom information
				_a.push_back(tmpa);
				// save id of the atom; for mapping of atom to residue
				ratoms.push_back(atomid);
				// map AtomKey to id for bonds
				key2index[*k]=atomid;
				atomid++;
			}
			_ratoms.push_back(ratoms);
			resid++;
		}
		_cres.push_back(cres);
	}
	// save bonds with the position in the atoms bond vector
	BondList bl(alBonds(atoms));
	for (BondList::iterator i=bl.begin(); i!=bl.end(); ++i)
	{
		BondList a1bnd((*i)->Atom1()->Bonds());
		BondList a2bnd((*i)->Atom2()->Bonds());
		BondList::iterator i1(find(a1bnd.begin(), a1bnd.end(), (*i)));
		BondList::iterator i2(find(a2bnd.begin(), a2bnd.end(), (*i)));
		#ifdef DEBUG2
		if (i1==a1bnd.end() || i2==a2bnd.end()) throw 
		PhoenixError("MOL::extractMOL: a1bnd or a2bnd boundary error in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
		#endif
		BOND tmpb(key2index[(*i)->Atom1()], key2index[(*i)->Atom2()], (*i)->Type(), (*i)->FFType(), distance(a1bnd.begin(), i1), distance(a2bnd.begin(), i2));
		_b.push_back(tmpb);
	}
}

ChainList MOL::createMOL()
{
	ChainList result;
	map<unsigned int, AtomKey> index2key;
	for (unsigned int i=0; i!=_c.size(); ++i)
	{
		// create chain
		ChainKey newchain=cCreate(_c[i]._name);
		newchain->Tag(_c[i]._tag);
		result.push_back(newchain);
		#ifdef DEBUG2
		if (i>=_cres.size()) throw 
		PhoenixError("MOL::createMOL: _cres boundary error in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
		#endif
		for (unsigned int j=0; j!=_cres[i].size(); ++j)
		{
			// create residue
			unsigned int rindex=_cres[i][j];
			#ifdef DEBUG2
			if (rindex>=_r.size()) throw 
			PhoenixError("MOL::createMOL: _r boundary error in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
			#endif
			ResidueKey newres=rCreate(newchain, _r[rindex]._rt);
			newres->Name(_r[rindex]._name);
			newres->Id(_r[rindex]._id);
			#ifdef DEBUG2
			if (rindex>=_ratoms.size()) throw 
			PhoenixError("MOL::createMOL: _ratoms boundary error in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
			#endif
			for (unsigned int k=0; k!=_ratoms[rindex].size(); ++k)
			{
				// create atom and
				unsigned int aindex=_ratoms[rindex][k];
				#ifdef DEBUG2
				if (aindex>=_a.size()) throw 
				PhoenixError("MOL::createMOL: _a boundary error in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
				#endif
				AtomKey newatom=aCreate(newres, _a[aindex]._ele, _a[aindex]._coord, _a[aindex]._geo);
				newatom->Name(_a[aindex]._name);
				newatom->FFType(_a[aindex]._fftype);
				newatom->Charge(_a[aindex]._charge);
				newatom->FormalCharge(_a[aindex]._fcharge);
				newatom->MOL2Type(_a[aindex]._mol2type);
				// map index to new AtomKey to id for bonds
				index2key[aindex]=newatom;
			}
		}
	}
	// create bond if all bonds, which need to be created previously, have been created
	queue<BOND> bonds;
	for (unsigned int i=0; i!=_b.size(); ++i)
	{
		#ifdef DEBUG2
		if (index2key.find(_b[i]._a1)==index2key.end() || index2key.find(_b[i]._a2)==index2key.end()) throw 
		PhoenixError("MOL::createMOL: index2key index error in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
		#endif
		AtomKey a1=index2key[_b[i]._a1];
		AtomKey a2=index2key[_b[i]._a2];
		if (a1->Bonds().size()==_b[i]._a1pos && a2->Bonds().size()==_b[i]._a2pos)
		{
			BondKey newbond=bCreate(a1, a2, _b[i]._bt);
			newbond->FFType(_b[i]._fftype);
		}
		else
		{
			bonds.push(_b[i]);
		}
	}
	//repeat until all bonds have bean created
	#ifdef DEBUG2
	unsigned int size=bonds.size();
	unsigned int counter=size;
	#endif
	while (!bonds.empty())
	{
		BOND front=bonds.front(); bonds.pop();
		AtomKey a1=index2key[front._a1];
		AtomKey a2=index2key[front._a2];
		if (a1->Bonds().size()==front._a1pos && a2->Bonds().size()==front._a2pos)
		{
			BondKey newbond=bCreate(a1, a2, front._bt);
			newbond->FFType(front._fftype);
		}
		else
		{
			bonds.push(front);
		}
		#ifdef DEBUG2
		if (size!=bonds.size())
		{
			size=bonds.size();
			counter=size;
		}
		else
		{
			if (!counter) throw 
			PhoenixError("MOL::createMOL: traped in endless loop in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
			--counter;
		}
		#endif
	}	
	return result;
}


/////////
// SAF
//////



////////////////////////////////////
// end graph structure functions
//////////////////////////////////


////////////////////////////////
// graph propertie functions
//////////////////////////////




int GraphDistance(AtomKey a1, AtomKey a2)
{
	map<AtomKey, bool> visited;
	int result=0;
	visited[a1]=1;
	
	AtomList disc;
	disc.push_back(a1);
	while (!disc.empty())
	{
		++result;
		AtomList newdisc;
		for (AtomList::iterator i=disc.begin(); i!=disc.end(); ++i)
		{
			// all adjent atoms
			AtomList tmp=aAdjAtoms(*i);
			for (AtomList::iterator j=tmp.begin(); j!=tmp.end(); ++j)
			{
				// if the atom is new
				if (!visited[*j])
				{
					if ((*j)==a2) return result;
					// set graph discovered to 1
					visited[*j]=1;
					// and store in newdisc
					newdisc.push_back(*j);
				}
			}
		}
		disc=newdisc;
	}
	return 0;
	
}

struct Sol
{
	// maps query atom and bonds to graph atoms and bonds
	map<unsigned int, AtomKey> atommatch;
	map<unsigned int, BondKey> bondmatch;
	// visited graph atoms and bonds
	map<AtomKey, bool> atom_v;
	map<BondKey, bool> bond_v;
	// visited query atoms and bonds
	vector<bool> q_atom_v;
	vector<bool> q_bond_v;
	// query atoms for the next iteration
	queue<unsigned int> disc;
};

class Increment
{
public:
	Increment(): value(0) {}
	unsigned int operator()()
	{
		return value++;
	}
private:
	int value;
};

void printsol(vector<Sol>& sol)
{
	for (vector<Sol>::iterator i=sol.begin(); i!=sol.end(); ++i)
	{
		cout<<aResidue(i->atommatch[0])->Name()<<" ";
		for (map<unsigned int, AtomKey>::iterator j=i->atommatch.begin();
		j!=i->atommatch.end(); ++j) cout <<j->second->Name()<<" ";
		cout<<endl; 
	}	
	cout<<endl;
	cout<<endl;	
}

/*void GraphMatch(GraphQuery query)
{
	// all finished solutions
	vector<Sol> result;
	// find all possible starting points for the query
	vector<Sol> allsol;
	{
	vector<bool> atom_v(query._a.size(), false);
	vector<bool> bond_v(query._b.size(), false);
	// start at first atom of the query
	queue<unsigned int> disc;
	disc.push(0);
	atom_v[0]=true;
	// check all atoms of the graph
	AtomList a=Atoms();
	for (AtomList::iterator i=a.begin(); i!=a.end(); ++i)
	{
		if (query._a[0].Match(**i))
		{
			Sol sol;
			sol.atommatch[0]=*i;
			sol.atom_v[*i]=true;
			sol.q_atom_v=atom_v;
			sol.q_bond_v=bond_v;
			sol.disc=disc;
			allsol.push_back(sol);
		}
	}
	}
	
	while (allsol.size())
	{
		vector<Sol> newsol;
		for (vector<Sol>::iterator i=allsol.begin(); i!=allsol.end(); ++i)
		{
			Sol sol=*i;
			// if disc is empty the query traversal is finished
			// and if the solution is still in allsol it is a valid solution
			if (sol.disc.empty())
			{
				result.push_back(sol);
				continue; 
			}
			
			// prepare query
			// get first atom of disc
			unsigned int query_a=sol.disc.front();
			sol.disc.pop();
			// find unvisited bonds of the query_a
			vector <unsigned int> q_edges;
			vector <unsigned int> q_backedges;
			// and find new atoms in the query
			vector <unsigned int> q_atom_on_edge;
			// and backedge atoms
			vector <unsigned int> q_atom_on_backedge;
			for (unsigned int j=0; j!=query._a_l_bond[query_a].size(); ++j)
			{
				// unvisited?
				if (!sol.q_bond_v[query._a_l_bond[query_a][j]])
				{
					// check if the bond is a backedge
					// bond is a backedge if the connected atom is visited
					if (sol.q_atom_v[query._a_l_atom[query_a][j]])
					{
						q_backedges.push_back(query._a_l_bond[query_a][j]);
						sol.q_bond_v[query._a_l_bond[query_a][j]]=true;
						q_atom_on_backedge.push_back(query._a_l_atom[query_a][j]);
					}
					else
					{
						q_edges.push_back(query._a_l_bond[query_a][j]);
						sol.q_bond_v[query._a_l_bond[query_a][j]]=true;
						q_atom_on_edge.push_back(query._a_l_atom[query_a][j]);
						sol.q_atom_v[query._a_l_atom[query_a][j]]=true;
					}
				}
			}
			
			// this is the graph atom which is in progress now
			AtomKey graph_a=sol.atommatch[query_a];
			
			// check backedges first
			vector <unsigned int> q_backedges_notmatched;
			for (unsigned int j=0; j!=q_backedges.size(); ++j)
			{
				// the backedges goes from query_a to q_atom_on_backedge[j]
				// check if there is allready a matched atom for q_atom_on_backedge[j]
				bool match=false;
				if (sol.atommatch.find(q_atom_on_backedge[j])!=sol.atommatch.end())
				{
					// is there a bond in the graph?
					BondList b=aFormBond(graph_a, sol.atommatch[q_atom_on_backedge[j]]);
					// if there is a bond check it
					if (b.size())
					{
						if (query._b[q_backedges[j]].Match(**b.begin()))
						{
							// sucessfull backedge match
							match=true;
							sol.bondmatch[q_backedges[j]]=*b.begin();
							sol.bond_v[*b.begin()]=true;
						}
					}
				}
				if (!match && !query._b[q_backedges[j]].Opt())
				{
					q_backedges_notmatched.push_back(q_backedges[j]);
				}
			}
			
			// if there are unmatched not optional backedges
			// the sollution can be ommited
			if (q_backedges_notmatched.size())
			{
				continue;
			}
			
			// prepare graph
			// find unvisited bonds in the graph
			BondList gb=aBonds(graph_a);
			vector<BondKey> g_edges;
			// and find new atoms in the graph
			vector<AtomKey> g_atom_on_edge;
			bool unmatched_backedges=false;
			for (BondList::iterator j=gb.begin(); j!=gb.end(); ++j)
			{
				if (!sol.bond_v[*j])
				SAF{
					AtomKey connected_atom=((*j)->Atom1()==graph_a ? (*j)->Atom2() : (*j)->Atom1() );
					if (!sol.atom_v[connected_atom])
					{
						g_edges.push_back(*j);
						g_atom_on_edge.push_back(connected_atom);
					}
					else
					{
						unsigned int q_connected_atom=query_a;
						for (map<unsigned int, AtomKey>::iterator k=sol.atommatch.begin();
							k!=sol.atommatch.end(); ++k)
						{
							if (k->second==connected_atom)
							{
								q_connected_atom=k->first;
								break;
							}
						}
						// there are unmatched backedges on this atom
						if (query._a[query_a].Strict() || query._a[q_connected_atom].Strict())
						{
							// this sollution can be ommited
							unmatched_backedges=true;
							break;
						}
					}
				}
			}
			//cout<<g_edges.size()<<endl;

			if (unmatched_backedges)
			{
				continue;
			}
			
			// find matches
			// in the case of strict checking all query bonds must be matched by graph bonds
			if (query._a[query_a].Strict())
			{
				// if there are more bonds in the graph the solution can be ommited
				if (g_edges.size() > q_edges.size())
				{
					continue;
				}
				// test permutations
				vector<unsigned int> perm;
				generate_n(back_inserter(perm), q_edges.size(), Increment());
				do
				{
					Sol thissol=sol;
					bool matches=true;
					unsigned int j=0;
					// match the graph edges with the query edge permutation
					for (;j!=g_edges.size(); ++j)
					{
						if (query._b[q_edges[perm[j]]].Match(*g_edges[j]) && 
						    query._a[q_atom_on_edge[perm[j]]].Match(*g_atom_on_edge[j]))
						{
							thissol.bondmatch[q_edges[perm[j]]]=g_edges[j];
							thissol.bond_v[g_edges[j]]=true;
							thissol.atommatch[q_atom_on_edge[perm[j]]]=g_atom_on_edge[j];
							thissol.atom_v[g_atom_on_edge[j]]=true;
							thissol.disc.push(q_atom_on_edge[j]);
						}
						else
						{
							matches=false;
							break;
						}
					}
					// remaining query edges need to be optional
					for (;j!=q_edges.size(); ++j)
					{
						if (!query._b[q_edges[perm[j]]].Opt())
						{
							matches=false;
							break;
						}
					}
					if (matches)
					{
						newsol.push_back(thissol);
					}
				}
				while (next_permutation(perm.begin(), perm.end()));
			}
			// in the case of non strict checking at least 
			// all existing not optional query edges must be matched
			else
			{
				if (q_edges.empty())
				{
					newsol.push_back(sol);
				}
				else
				{
					// find not optional query edges
					vector <unsigned int> nopt_q_edges;
					vector <unsigned int> nopt_q_atom_on_edge;
					vector <unsigned int> opt_q_edges;
					vector <unsigned int> opt_q_atom_on_edge;
					for (unsigned int i=0; i!=q_edges.size(); ++i)
					{
						if (query._b[q_edges[i]].Opt())
						{
							opt_q_edges.push_back(q_edges[i]);
							opt_q_atom_on_edge.push_back(q_atom_on_edge[i]);
						}
						else
						{
							nopt_q_edges.push_back(q_edges[i]);
							nopt_q_atom_on_edge.push_back(q_atom_on_edge[i]);
						}
					}
					
					
					// if there are less graph edges than not optional query edges
					// this solution can be ommited
					if (g_edges.size() < nopt_q_edges.size())
					{
						continue;
					}
					
					// test permutations
					vector<unsigned int> perm;
					generate_n(back_inserter(perm), g_edges.size(), Increment());
					do
					do
					{
						Sol thissol=sol;
						bool matches=true;
						unsigned int j=0;
						// match the non optional query edges with the g_edges permutation
						for (;j!=nopt_q_edges.size(); ++j)
						{
							if (query._b[nopt_q_edges[j]].Match(*g_edges[perm[j]]) && 
							query._a[nopt_q_atom_on_edge[j]].Match(*g_atom_on_edge[perm[j]]))
							{
								thissol.bondmatch[nopt_q_edges[j]]=g_edges[perm[j]];
								thissol.bond_v[g_edges[perm[j]]]=true;
								thissol.atommatch[nopt_q_atom_on_edge[j]]=g_atom_on_edge[perm[j]];
								thissol.atom_v[g_atom_on_edge[perm[j]]]=true;
								thissol.disc.push(nopt_q_atom_on_edge[j]);
							}
							else
							{
								matches=false;
								break;
							}
						}
						
						if (matches)
						{
							// from here everything is a valid match
							unsigned int k=0;
							for (;(j!=(nopt_q_edges.size()+opt_q_edges.size()) &&
							       j!=g_edges.size()); ++j)
							{
								if (query._b[opt_q_edges[k]].Match(*g_edges[perm[j]]) && 
								    query._a[opt_q_atom_on_edge[k]].Match(*g_atom_on_edge[perm[j]]))
								{
									thissol.bondmatch[opt_q_edges[k]]=g_edges[perm[j]];
									thissol.bond_v[g_edges[perm[j]]]=true;
									thissol.atommatch[opt_q_atom_on_edge[k]]=g_atom_on_edge[perm[j]];
									thissol.atom_v[g_atom_on_edge[perm[j]]]=true;
									thissol.disc.push(opt_q_atom_on_edge[k]);
								}																
								++k;
							}
						}
						if (matches)
						{
							newsol.push_back(thissol);
						}
					}
					while (next_permutation(perm.begin(), perm.begin()+nopt_q_edges.size()+opt_q_edges.size()));
					while (btb::next_combination(perm.begin(), perm.begin()+nopt_q_edges.size()+opt_q_edges.size(), perm.end()));
				}
			}
		}
		allsol=newsol;
	}
	cout<<result.size()<< " ergebnisse"<<endl;
	printsol(result);
}*/
////////////////////////////////////
// end graph propertie functions
//////////////////////////////////


//////////////////////////////////////////
// automatic type assignment functions
////////////////////////////////////////

//////////////////////////////////////////////
// end automatic type assignment functions
////////////////////////////////////////////



// to be removed

AtomKey aAdjacentAtom(AtomKey a1, BondKey b)
{
        return b->Atom1()==a1 ? b->Atom2() : b->Atom1();
}



}
