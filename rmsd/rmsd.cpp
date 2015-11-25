#include "mgraph/graph.hpp"


using namespace std;
using namespace mgraph;

struct Priority
{
	Priority() {}
	Priority(AtomKey a);
	bool operator== (Priority p) const;
	bool operator> (Priority p) const;
	vector<vector<unsigned int> > _plist;
};

Priority::Priority(AtomKey a)
{
	map<AtomKey, bool> visited;
	AtomList disc;
	
	vector<unsigned int> tmp(1, (unsigned int)(a)->Element());
	_plist.push_back(tmp);
	visited[a]=1;
	disc.push_back(a);
	while (!disc.empty())
	{
		vector<unsigned int> tmp2;
		AtomList newdisc;
		for (AtomList::iterator i=disc.begin(); i!=disc.end(); ++i)
		{
			AtomList adj(aAdjAtoms(*i));
			for (AtomList::iterator j=adj.begin(); j!=adj.end(); ++j)
			{
				if (!visited[*j])
				{
					newdisc.push_back(*j);
					tmp2.push_back((unsigned int)(*j)->Element());
					visited[*j]=1;
				}
			}
		}
		sort(tmp2.begin(), tmp2.end(), greater<unsigned int>());
		_plist.push_back(tmp2);
		disc=newdisc;
	}
}

bool Priority::operator==(Priority p) const
{
	return _plist==p._plist;
}

bool Priority::operator> (Priority p) const
{
	return _plist>p._plist;
}

struct Path
{
	map<AtomKey, bool> visited;
	queue<AtomKey> disc;
	AtomList path;
};


class RMSD
{
public:
	RMSD(map<AtomKey, Priority>& p) : refmap(p) {}
	void initRef(const ChainList& ref);
	void calcRMSD(const ChainList& mol);

	vector<AtomList> refpath;
	map<AtomKey, Priority> refmap;
};

struct equal_priority: public binary_function<AtomKey, AtomKey, bool> 
{
	equal_priority(map<AtomKey, Priority> &map) : _map(map) {}
	bool operator()(const AtomKey x, const AtomKey y) { return _map[x] == _map[y]; }
	map<AtomKey, Priority> &_map;
};

struct greater_priority: public binary_function<AtomKey, AtomKey, bool> 
{
	greater_priority(map<AtomKey, Priority> &map) : _map(map) {}
	bool operator()(const AtomKey x, const AtomKey y) { return _map[x] > _map[y]; }
	map<AtomKey, Priority> &_map;
};

void RMSD::initRef(const ChainList& ref)
{
	AtomList refatoms(clAtoms(ref));
	if (!refatoms.size()) throw PhoenixError("RMSD::initRef: no atoms.");
	
	// priority less compare
	equal_priority ep(refmap);
	greater_priority gp(refmap);
	
	sort(refatoms.begin(), refatoms.end(), gp);
	
	AtomList start;
	{
		AtomList::iterator i=refatoms.begin();
		start.push_back(*i);
		++i;
		for (; i!=refatoms.end(); ++i)
		{
			if (ep((*i),*start.begin()))
			{
				start.push_back(*i);
			}
			else
			{
				break;
			}
		}
	}
	
	vector<Path> allpath;
	for (AtomList::iterator i=start.begin(); i!=start.end(); ++i)
	{
		Path tmp;
		tmp.path.push_back(*i);
		tmp.visited[*i]=1;
		tmp.disc.push(*i);
		allpath.push_back(tmp);
	}
	
	for (unsigned int i=0; i!=allpath.size(); ++i)
	{
		Path tmp=allpath[i];
		// start a bfs search
		while (!tmp.disc.empty())
		{
			AtomKey front=tmp.disc.front(); tmp.disc.pop();
			AtomList adj(aAdjAtoms(front));
			AtomList newatoms;
			for (AtomList::iterator j=adj.begin(); j!=adj.end(); ++j)
			{
				if (!tmp.visited[*j])
				{
					tmp.visited[*j]=1;
					newatoms.push_back(*j);
				}
			}
			sort(newatoms.begin(), newatoms.end(), gp);
			list<AtomKey> newatomslist(newatoms.begin(), newatoms.end());
			vector<AtomList> minipath(1);
			
			while (!newatomslist.empty())
			{
				list<AtomKey>::iterator j=newatomslist.begin();
				AtomList equalatoms;
				for (list<AtomKey>::iterator k=newatomslist.begin(); k!=newatomslist.end(); ++k)
				{
					if (ep(*j, *k))
					{
						equalatoms.push_back(*k);
					}
				}
				
				for (AtomList::iterator k=equalatoms.begin(); k!=equalatoms.end(); ++k)
				{
					tmp.disc.push(*k);
					newatomslist.remove(*k);
				}
				
				if (equalatoms.size()==1)
				{
					for (vector<AtomList>::iterator k=minipath.begin(); k!=minipath.end(); ++k)
					{
						k->push_back(*equalatoms.begin());
					}
				}
				else
				{
					// need to build up permutation
					sort(equalatoms.begin(), equalatoms.end());
					vector<AtomList> tmpminipath;
					for (vector<AtomList>::iterator k=minipath.begin(); k!=minipath.end(); ++k)
					{
						AtomList seq=*k;
						copy(equalatoms.begin(), equalatoms.end(), back_inserter(seq));
						tmpminipath.push_back(seq);
						while (next_permutation(equalatoms.begin(), equalatoms.end()))
						{
							AtomList seq2=*k;
							copy(equalatoms.begin(), equalatoms.end(), back_inserter(seq2));
							tmpminipath.push_back(seq2);
						}
					}
					minipath=tmpminipath;
				}
			}
			
			if (minipath.size()==1)
			{
				for (AtomList::iterator k=minipath[0].begin(); k!=minipath[0].end(); ++k)
				{
					tmp.path.push_back(*k);
				}
			}
			else
			{
				Path backup=tmp;
				for (AtomList::iterator k=minipath[0].begin(); k!=minipath[0].end(); ++k)
				{
					tmp.path.push_back(*k);
				}
				for (unsigned int l=1; l!=minipath.size(); ++l)
				{
					Path tmp2=backup;
					for (AtomList::iterator k=minipath[l].begin(); k!=minipath[l].end(); ++k)
					{
						tmp2.path.push_back(*k);
					}
					allpath.push_back(tmp2);
				}
			}
		}
		allpath[i]=tmp;
	}
	
	for (unsigned int i=0; i!=allpath.size(); ++i)
	{
		refpath.push_back(allpath[i].path);
	}

}

void RMSD::calcRMSD(const ChainList& mol)
{
	AtomList refatoms(clAtoms(mol));
	if (!refatoms.size()) throw PhoenixError("RMSD::initRef: no atoms.");
	
	// priority less compare
	equal_priority ep(refmap);
	greater_priority gp(refmap);
	
	sort(refatoms.begin(), refatoms.end(), gp);
	
	AtomList start;
	start.push_back(*refatoms.begin());
	
	Path molpath;
	molpath.path.push_back(*refatoms.begin());
	molpath.visited[*refatoms.begin()]=1;
	molpath.disc.push(*refatoms.begin());
	
	// start a bfs search
	while (!molpath.disc.empty())
	{
		AtomKey front=molpath.disc.front(); molpath.disc.pop();
		AtomList adj(aAdjAtoms(front));
		AtomList newatoms;
		for (AtomList::iterator j=adj.begin(); j!=adj.end(); ++j)
		{
			if (!molpath.visited[*j])
			{
				molpath.visited[*j]=1;
				newatoms.push_back(*j);
			}
		}
		sort(newatoms.begin(), newatoms.end(), gp);
		list<AtomKey> newatomslist(newatoms.begin(), newatoms.end());
		vector<AtomList> minipath(1);
		
		while (!newatomslist.empty())
		{
			list<AtomKey>::iterator j=newatomslist.begin();
			AtomList equalatoms;
			for (list<AtomKey>::iterator k=newatomslist.begin(); k!=newatomslist.end(); ++k)
			{
				if (ep(*j, *k))
				{
					equalatoms.push_back(*k);
				}
			}
			
			for (AtomList::iterator k=equalatoms.begin(); k!=equalatoms.end(); ++k)
			{
				molpath.disc.push(*k);
				newatomslist.remove(*k);
			}
			
			if (equalatoms.size()==1)
			{
				for (vector<AtomList>::iterator k=minipath.begin(); k!=minipath.end(); ++k)
				{
					k->push_back(*equalatoms.begin());
				}
			}
			else
			{
				// need to build up permutation
				sort(equalatoms.begin(), equalatoms.end());
				vector<AtomList> tmpminipath;
				for (vector<AtomList>::iterator k=minipath.begin(); k!=minipath.end(); ++k)
				{
					AtomList seq=*k;
					copy(equalatoms.begin(), equalatoms.end(), back_inserter(seq));
					tmpminipath.push_back(seq);
				}
				minipath=tmpminipath;
			}
		}
		
		if (minipath.size()==1)
		{
			for (AtomList::iterator k=minipath[0].begin(); k!=minipath[0].end(); ++k)
			{
				molpath.path.push_back(*k);
			}
		}
		else
		{
			throw PhoenixError("TERROR");
		}
	}
	
	AtomList molpathsorted=molpath.path;
	
	Float_t rmsd=999999999.0;
	
	for (unsigned int i=0; i!=refpath.size(); ++i)
	{
		vector<Priority> refpriority;
		for (AtomList::iterator j=refpath[i].begin(); j!=refpath[i].end(); ++j)
		{
			refpriority.push_back(refmap[*j]);
		}
		vector<Priority> molpri;
		for (AtomList::iterator j=molpathsorted.begin(); j!=molpathsorted.end(); ++j)
		{
			molpri.push_back(refmap[*j]);
		}
		if (!(refpriority==molpri)) continue;
		
		Float_t tmp_rmsd=0;
		for (unsigned int j=0; j!=refpath[i].size(); ++j)
		{
			tmp_rmsd+=(refpath[i][j]->Coord()-molpathsorted[j]->Coord()).length2();
		}
		tmp_rmsd=sqrt(tmp_rmsd/refpath[i].size());
		if(rmsd>tmp_rmsd) rmsd=tmp_rmsd;
	}
	cout << "RMSD = " << rmsd <<endl;

}

int main(int argc, char *argv[])
{
	bool del_h(true);
	if (string(argv[1])=="-h")
	{
		cout<<"Use hydrogen atoms during RMSD calculation."<<endl;
		del_h=false;
	}
	ChainList lig1, lig2;
	try
	{
		lig1=ReadMOL2(argv[1]);      //read ligand 1
		lig2=ReadMOL2(argv[2]);      //read ligand 2
	}
	catch(...)
	{
		cout<<"Usage: rmsd molecule1 molecule2"<<endl;
		cout<<"Calculates RMSD of heavy atoms only."<<endl;
		cout<<"Usage: rmsd -h molecule1 molecule2"<<endl;
		cout<<"Calculates RMSD of heavy atoms and hydrogen atoms."<<endl;
	}
		
	// delete lonepair
	AtomList allatoms=Atoms();
	for (AtomList::iterator i=allatoms.begin(); i!=allatoms.end(); ++i)
	{
		if ((*i)->Element()==ELE::LP) aDelete(*i);
	}
	
	//annotate with priority
	map<AtomKey, Priority> priorities;
	for (AtomList::iterator i=allatoms.begin(); i!=allatoms.end(); ++i)
	{
		priorities[(*i)]=Priority(*i);
	}
	
	if (del_h)
	{
		for (AtomList::iterator i=allatoms.begin(); i!=allatoms.end(); ++i)
		{
			if ((*i)->Element()==ELE::H) aDelete(*i);
		}
	}
	
	RMSD rmsd(priorities);
	rmsd.initRef(lig1);
	rmsd.calcRMSD(lig2);
}
