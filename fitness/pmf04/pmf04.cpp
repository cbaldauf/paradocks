#include "pmf04.hpp"
#include <fstream>

using namespace mgraph;

//TODO 	
namespace paradocks 
{
	PMF04::PMF04(vector<string>& parameter) : Fitness(parameter), _result(_termscount)
	{
		_SF1=1;	
		_SF2=2;
		_SF3=1;
		_vdw_range=16;
		_vdw_repulsion_cutoff=15;
		_grid_size=0.3;
	}
	
	Float_t PMF04::vdW_type[46]=
	{
		1, // H 0
		2.32556, // C tet 1
		2.1839, // C tri ar 2
		2.01153, // C tri re 3
		1.88987, // C tri 4
		1.8, // C lin 5
		1.83461, // N tet 6
		1.73761, // N tri ar 7
		1.74601, // N tri re 8
		1.72442, // N tri 9
		1.75673, // N lin 10
		1.7129, // O tet 11
		1.57094, // O tri ar 12
		1.57675, // O tri re 13
		1.57801, // O tri 14
		1.55, // O lin 15
		2.11339, // S tet 16
		1.96368, // S tri ar 17
		2, // S tri re 18
		2, // S tri 19
		2, // S lin 20
		1.95237, // P tet 21
		2, // P tri ar 22
		2, // P tri re 23
		2, // P tri 24
		2, // P lin 25
		1.4997, // F 26
		1.71409, // Cl 27
		1.83567, // Br 28
		2.05, // I 29
		2, // Si 30
		1.25, // metall K 31
		1.25, // metall U 32
		1.25, // metall Li 33
		1.25, // metall Na 34
		1.20589, // metall Ca 35
		1.25, // metall Cd 36
		1.25, // metall Co 37
		1.25, // metall Cu 38
		1.2024, // metall Fe 39
		1.25, // metall Hg 40
		1.21929, // metall Mg 41
		1.25, // metall Mn 42
		1.25, // metall Ni 43
		//1.20855, // metall Zn 44
		0.20855, // metall Zn 44
		1.25, // metall Al 45
	};

void PMF04::initProtein(const ChainList& prot, const Vec3_t& as, Float_t r)
{
	CorrectGraphByTopology(prot);
	Fitness::initProtein(prot, as, r);
	
	// atoms for vdw interaction
	AtomList patoms_tmp;
	for (unsigned int i=0; i!=_patoms.size(); ++i)
	{
		if (_patoms[i]->Coord().x()<_center.x()-_radius-9 || _patoms[i]->Coord().y()<_center.y()-_radius-9 || _patoms[i]->Coord().z()<_center.z()-_radius-9 ||
			_patoms[i]->Coord().x()>_center.x()+_radius+9 || _patoms[i]->Coord().y()>_center.y()+_radius+9 || _patoms[i]->Coord().z()>_center.z()+_radius+9)
		{
			continue;
		}
		unsigned int type(vdWType(_patoms[i]));
		if (type)
		{
			_vdw_patom.push_back(i);
			_pvdw.push_back(vdW_type[type]);
			patoms_tmp.push_back(_patoms[i]);
		}
	}
 	// set up neighboorhood list for vdw
 	_vdw_nh.setupList(patoms_tmp, 8, 1);
 	// set up the vdw grid
 	_vdw_grid.setupGrid(as-Vec3_t(r,r,r), as+Vec3_t(r,r,r), _grid_size, 46, 5);	
 	// set up neighboorhood list for sp (with hydrogen atoms)
 	_sp_nh.setupList(_patoms, 9, 1);
 	// set up statistical potential grid
 	_sp_grid.setupGrid(as-Vec3_t(r,r,r), as+Vec3_t(r,r,r), _grid_size, 31, 5);


	set<AtomKey> ua(_patoms.begin(), _patoms.end());
	// search hydrogen HH
	set<AtomKey>::iterator i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::H))
		{
			(*i)->FFType(17);
			ua.erase(i++);
			continue;
		}
		++i;
	}
	//search nonpolar and polar carbon aromatic 	cP and cF
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::C) && (*i)->Geometry() == GEO(GEO::tri))
		{
			BondList bl = aBonds(*i);
			AtomList al = aAdjAtoms(*i);
			bool polar=false, arom=false;
			for(BondList::iterator bl_it=bl.begin();bl_it!=bl.end();++bl_it)
			{
				if((*bl_it)->Type() == BT(BT::ar))
				{
					arom=true;
					break;
				}
			}
			if(arom)
			{
				for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
				{
					if((*al_it)->Element() != ELE(ELE::C) && (*al_it)->Element() != ELE(ELE::H))
					{
						polar=true;
						break;
					}
				}
				if(polar)
				{
					(*i)->FFType(4);			/// polar carbon aromatic cP
					ua.erase(i++);
					continue;
				}
				else
				{
					(*i)->FFType(3);			/// nonpolar carbon aromatic cF
					ua.erase(i++);
					continue;
				}
			}
		}
		++i;
	}
	//search nonpolar aliphatic carbon CF
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::C) && (*i)->Geometry() == GEO(GEO::tet))
		{
			AtomList al = aAdjAtoms(*i);
			bool polar=false;
			for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
			{
				if((*al_it)->Element() != ELE(ELE::C) && (*al_it)->Element() != ELE(ELE::H))
				{
					polar=true;
					break;
				}
			}
			if(!polar)
			{
				(*i)->FFType(1);			/// nonpolar aliphatic carbon CF
				ua.erase(i++);
				continue;
			}
		}
		++i;
	}
	// search negatively charged oxygen (e.g. carboxylate)
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::O) && aBonds(*i).size() == 1)
		{
			BondKey bindung= *(aBonds(*i).begin());
			if(bindung->Type() == BT(BT::s) || bindung->Type() == BT(BT::re))
			{
				bool amid=false;
				AtomKey nachbar = aAdjacentAtom(*i,bindung);
				/// Annahme : negativ gel. O kann nur an C haengen
				AtomList al = aAdjAtoms(nachbar);
				for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
				{
					if((*al_it)->Element() == ELE(ELE::N))
					{
						amid=true;
						break;
					}
				}
				if(amid==false)
				{
					(*i)->FFType(10);
					// search carbon bonded to a negatively charged oxygen
					if(nachbar->Element() == ELE(ELE::C))
					{
						(nachbar)->FFType(5);
						if (ua.find(nachbar)!=ua.end()) ua.erase(nachbar);
					}
					ua.erase(i++);
					continue;
				}
			}
		}
		++i;
	}
	// search positively charged nitrogen (e.g. NH3+, guanidino group)
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::N))
		{
			int resStab = 0, anz_N=0;
			bool guanidino=false;
			BondList bl = aBonds(*i);
			for(BondList::iterator bl_it=bl.begin();bl_it!=bl.end();++bl_it)
			{
				if((*bl_it)->Type() == BT(BT::re))
				{
					++resStab;
				}
			}
			AtomList al = aAdjAtoms(*i);
			for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
			{
				if((*al_it)->Element() == ELE(ELE::C))
				{
					BondList bl = aBonds(*al_it);
					for(BondList::iterator bl_it=bl.begin();bl_it!=bl.end();++bl_it)
					{
						if(aAdjacentAtom(*al_it,*bl_it)->Element() == ELE(ELE::N) && (*bl_it)->Type() == BT(BT::re))
						{
							++anz_N;
						}
					}
					if(anz_N == 3)
					{	
						guanidino=true;
						break;
					}
					else anz_N=0;
				}
			}
			if((aBonds(*i).size() == 4) || (resStab > 1) || guanidino)
			{
				(*i)->FFType(7);
				AtomList al = aAdjAtoms(*i);
				for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
				{
					if((*al_it)->Element() == ELE(ELE::C))
					{
						if (ua.find(*al_it)!=ua.end()) 
						{
							(*al_it)->FFType(6);
							ua.erase(*al_it);
						}
					}
				}
				ua.erase(i++);
				continue;
			}
		}
		++i;
	}
	// search aliphatic sp2 or sp3 carbon bonded to atoms other than carbon or hydrogen (e.g. backbone C or C_alpha) polar!!!
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::C) && ((*i)->Geometry() == GEO(GEO::tri) || (*i)->Geometry() == GEO(GEO::tet)))
		{
			bool hallo=false;
			AtomList bl = aAdjAtoms(*i);
			for(AtomList::iterator al_it=bl.begin();al_it!=bl.end();++al_it)
			{
				if((*al_it)->Element() != ELE(ELE::C) && (*al_it)->Element() != ELE(ELE::H))
				{
					hallo=true;
					break;
				}
			}
			if(hallo)	
			{
				(*i)->FFType(2);
				ua.erase(i++);
				continue;
			}
		}
		++i;
	}
	// search nitrogen	ND and NA 
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::N))
		{
			bool wasserstoff=false;
			AtomList al = aAdjAtoms(*i);
			for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
			{
				if((*al_it)->Element() == ELE(ELE::H))
				{
					wasserstoff=true;
					break;
				}
			}
			// search nitrogen as a hydrogen bond donor
			if(wasserstoff)
			{
				(*i)->FFType(8);
				ua.erase(i++);
				continue;
			}
			// search nitrogen as a hydrogen bond acceptor
			else
			{
				(*i)->FFType(16);
				ua.erase(i++);
				continue;
			}
		}
		++i;
	}
	//search water oxygen	OW
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::O) && aBonds(*i).size()==2)
		{
			int wasserstoff=0;
			AtomList al = aAdjAtoms(*i);
			for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
			{
				if((*al_it)->Element() == ELE(ELE::H)) ++wasserstoff;
			}
			if(wasserstoff==2)
			{
				(*i)->FFType(13);
				ua.erase(i++);
				continue;
			}
		}
		++i;
	}
	// search oxygen as hydrogen bond acceptor and donor	OD and OA
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::O) && aBonds(*i).size() == 1)
		{
			// hb acceptor 
			if((*aBonds(*i).begin())->Type() == BT(BT::d) || (*aBonds(*i).begin())->Type() == BT(BT::re))
			{
				(*i)->FFType(11);
				ua.erase(i++);
				continue;
			}
		}
		// hb donor
		if((*i)->Element() == ELE(ELE::O))
			{
				(*i)->FFType(12);
				ua.erase(i++);
				continue;
			}
		++i;
	}
	// search sulfur as hydrogen bond donor and acceptor	SD and SA
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::S))
		{
			bool wasserstoff=false;
			AtomList al = aAdjAtoms(*i);
			for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
			{
				if((*al_it)->Element() == ELE(ELE::H))
				{	
					wasserstoff=true;
					break;
				}
			}
			// hb donor
			if(wasserstoff)
			{
				(*i)->FFType(15);
				ua.erase(i++);
				continue;
			}	
			// hb acceptor
			else
			{
				(*i)->FFType(14);
				ua.erase(i++);
				continue;
			}
		}
		++i;
	}
	// search metal ion (Zn, Mg, Ca, Fe, Mn, K)	ME
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::Zn) || (*i)->Element() == ELE(ELE::Mg) || (*i)->Element() == ELE(ELE::Ca) || (*i)->Element() == ELE(ELE::Fe) || (*i)->Element() == ELE(ELE::Mn) || (*i)->Element() == ELE(ELE::K))
		{
			(*i)->FFType(9);
			ua.erase(i++);
			continue;
		}
		++i;
	}
}
void PMF04::initLigand(const ChainList& lig)
{
	CorrectGraphByTopology(lig);
	Fitness::initLigand(lig, 2);
	
	set<AtomKey> ua(_latoms.begin(), _latoms.end());
	set<AtomKey>::iterator i=ua.begin();
	// search hydrogen
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::H))
		{
			(*i)->FFType(25);
			ua.erase(i++);
			continue;
		}
		++i;
	}
	// search nonpolar and polar carbon aromatic
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::C) && (*i)->Geometry() == GEO(GEO::tri))
		{
			BondList bl = aBonds(*i);
			AtomList al = aAdjAtoms(*i);
			bool polar=false;
			bool arom=false;
			for(BondList::iterator bl_it=bl.begin();bl_it!=bl.end();++bl_it)
			{
				if((*bl_it)->Type() == BT(BT::ar))
				{
					arom=true;
					break;
				}
			}
			for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
			{
				if((*al_it)->Element() != ELE(ELE::C) && (*al_it)->Element() != ELE(ELE::H))
				{
					polar=true;
					break;
				}
			}
			if(arom && !polar)
			{
				(*i)->FFType(3);
				ua.erase(i++);
				continue;
			}
			if(arom && polar)
			{
				(*i)->FFType(4);
				ua.erase(i++);
				continue;
			}
		}
		++i;
	}
	// search sulfur
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::S))
		{
			// search sulfur in sulfone or sulfoxid
			if(aBonds(*i).size() > 2)
			{
				(*i)->FFType(19);
				ua.erase(i++);
				continue;
			}
			AtomList al = aAdjAtoms(*i);
			bool wasserstoff=false;
			for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
			{
				if((*al_it)->Element() == ELE(ELE::H))
				{
					wasserstoff=true;
					break;
				}
			}
			// search sulfur as hydrogen bond donor
			if(wasserstoff)
			{
				(*i)->FFType(24);
				ua.erase(i++);
				continue;
			}
			// search sulfur as hydrogen bond acceptor
			else
			{
				(*i)->FFType(23);
				ua.erase(i++);
				continue;
			}
		}
		++i;	
	}
	// search oxygen bound to atoms other than carbon or hydrogen (e.g. phosphate)
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::O))
		{
			AtomList al = aAdjAtoms(*i);
			bool OS=false;
			for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
			{
				if((*al_it)->Element() == ELE(ELE::S) || (*al_it)->Element() == ELE(ELE::P))
				{
					OS=true;
					break;
				}
			}
			if(OS)
			{
				(*i)->FFType(20);
				ua.erase(i++);
				continue;
			}
		}
		++i;
	}
	// search negatively charged oxygen (e.g. carboxylate)
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::O) && aBonds(*i).size() == 1)
		{
			//Test auf Amidsauerstoff
			bool amidsauerstoff=false;
			BondKey bindung= *(aBonds(*i).begin());
			AtomKey nachbar = aAdjacentAtom(*i,bindung);
			if((nachbar)->Element() == ELE(ELE::C) && bindung->Type() == BT(BT::re))
			{
				AtomList al = aAdjAtoms(nachbar);
				for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
				{
					if((*al_it)->Element() == ELE(ELE::N) && (*aaFormBond(*al_it,nachbar).begin())->Type() == BT(BT::re))
					{
						amidsauerstoff=true;
						break;
					}
				}
			}

			if((bindung->Type() == BT(BT::s) || bindung->Type() == BT(BT::re)) && amidsauerstoff==false)
			{
				(*i)->FFType(16);
				// search carbon bonded to a negatively charged oxygen
				AtomKey nachbar = aAdjacentAtom(*i,bindung);
				if(nachbar->Element() == ELE(ELE::C))
				{
					(nachbar)->FFType(7);
					if (ua.find(nachbar)!=ua.end()) ua.erase(nachbar);
				}
				ua.erase(i++);
				continue;
			}
		}
		++i;
	}
	// search sp nitrogen
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::N) && (*i)->Geometry() == GEO(GEO::lin) && aBonds(*i).size() == 1)
		{
			(*i)->FFType(14);
			ua.erase(i++);
			continue;
		}
		++i;
	}
	// search sp carbon
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::C) && (*i)->Geometry() == GEO(GEO::lin))
		{
			(*i)->FFType(29);
			ua.erase(i++);
			continue;
		}
		++i;
	}
	// search nitrogen bound to atoms other than H or C and not of type ND
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::N))
		{
			bool notH=true,typeNS=false;
			AtomList nb=aAdjAtoms(*i);
			for(AtomList::iterator ai=nb.begin();ai!=nb.end();ai++)
			{
				if((*ai)->Element() == ELE(ELE::H)) {notH=false;break;}
				if((*ai)->Element() == ELE(ELE::S) || (*ai)->Element() == ELE(ELE::P))	typeNS=true;
			}
			if(notH && typeNS)	
			{
				(*i)->FFType(15);
				ua.erase(i++);
				continue;
			}
		}
		++i;
	}
	// search positively charged nitrogen (e.g. NH3+, guanidino group)
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::N))
		{
			int resStab = 0, anz_N=0, anzBdg=0,sauerstoff=0;
			bool guanidino=false;
			BondList bl = aBonds(*i);
			for(BondList::iterator bl_it=bl.begin();bl_it!=bl.end();++bl_it)
			{
				if((*bl_it)->Type() == BT(BT::re)) ++resStab;
				if((*bl_it)->Type() == BT(BT::s)) anzBdg+=1;
				if((*bl_it)->Type() == BT(BT::d)) anzBdg+=2;
				if((*bl_it)->Type() == BT(BT::t)) anzBdg+=3;
			}
			AtomList al = aAdjAtoms(*i);
			// test guanidino group
			for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
			{
				if((*al_it)->Element() == ELE(ELE::C))
				{
					BondList bl = aBonds(*al_it);
					for(BondList::iterator bl_it=bl.begin();bl_it!=bl.end();++bl_it)
					{
						if(aAdjacentAtom(*al_it,*bl_it)->Element() == ELE(ELE::N) && (*bl_it)->Type() == BT(BT::re) && aBonds(*i).size()==3)
						{
							++anz_N;
						}
					}
					if(anz_N == 3)
					{	
						guanidino=true;
						break;
					}
					else anz_N=0;
				}
				if((*al_it)->Element() == ELE(ELE::O)) sauerstoff+=1;
			}
			if((anzBdg >= 4) || (resStab > 1 && sauerstoff==2) || guanidino)
			{
				(*i)->FFType(9);
				AtomList al = aAdjAtoms(*i);
				for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
				{
					if((*al_it)->Element() == ELE(ELE::C))
					{
						if (ua.find(*al_it)!=ua.end()) 
						{
							(*al_it)->FFType(8);
							ua.erase(*al_it);
						}
					}
				}
				ua.erase(i++);
				continue;
			}
		}
		++i;
	}
	// search polar sp3 carbon bonded to an atom other than carbon or hydrogen
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::C) && (*i)->Geometry() == GEO(GEO::tet))
		{
			AtomList al = aAdjAtoms(*i);
			int counter=0;
			for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
			{
				if((*al_it)->Element() == ELE(ELE::C) || (*al_it)->Element() == ELE(ELE::H))
				{
					++counter;	
				}
			}
			if(counter==4)
			{
				(*i)->FFType(1);
				ua.erase(i++);
				continue;
			}
			// search polar carbon sp3
			(*i)->FFType(2);
			ua.erase(i++);
			continue;
		}
		++i;
	}
	// search carbon sp2
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::C))
		{
			bool polar=false;
			AtomList al = aAdjAtoms(*i);
			for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
			{
				if((*al_it)->Element() != ELE(ELE::C) && (*al_it)->Element() != ELE(ELE::H))
				{
					polar=true;
					break;
				}
			}
			// search polar carbon sp2 not aromatic
			if(polar)
			{
				(*i)->FFType(6);
				ua.erase(i++);
				continue;
			}
			// search nonpolar carbon sp2 not aromatic
			else
			{
				(*i)->FFType(5);
				ua.erase(i++);
				continue;
			}
		}
		++i;
	}
	// search oxygen as hydrogen bond donor (e.g. keto, amid oxygen)
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::O))
		{
			if(aBonds(*i).size()==1)
			{
				(*i)->FFType(17);
				ua.erase(i++);
				continue;
			}
			int wasserstoff=0;
			BondList bl = aBonds(*i);
			for(BondList::iterator bl_it=bl.begin();bl_it!=bl.end();++bl_it)
			{
				if(aAdjacentAtom(*i,*bl_it)->Element() == ELE(ELE::H))
				{
					++wasserstoff;
				}
			}
			if(wasserstoff==1)
			{
				(*i)->FFType(21);
				ua.erase(i++);
				continue;
			}
			else
			{
				(*i)->FFType(18);
				ua.erase(i++);
				continue;
			}
		}
		++i;
	}
	// search nitrogen
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::N))
		{
			int kohlenstoff=0,wasserstoff=0;
			bool arom=false;
			AtomList al = aAdjAtoms(*i);
			for(AtomList::iterator al_it=al.begin();al_it!=al.end();++al_it)
			{
				if((*al_it)->Element() == ELE(ELE::C)) ++kohlenstoff;
				if((*al_it)->Element() == ELE(ELE::H)) ++wasserstoff;
			}
			BondList bl = aBonds(*i);
			for(BondList::iterator bl_it=bl.begin();bl_it!=bl.end();++bl_it)
			{
				if((*bl_it)->Type() == BT(BT::ar))	{arom=true;break;}
			}
			// nitrogen as a hydrogen bond acceptor, not in a ring
			if(wasserstoff==0 && !isRingAtom(*i))
			{
				(*i)->FFType(11);
				ua.erase(i++);
				continue;
			}
			// nitrogen as a hydrogen bond donor, not in a ring
			if(wasserstoff>0 && !isRingAtom(*i) )
			{
				(*i)->FFType(12);
				ua.erase(i++);
				continue;
			}
			// planar nitrogen bonded  to 2 or 3 carbons but not to a hydrogen NP
			if((arom || ((*i)->Geometry() == GEO(GEO::tri))) && (wasserstoff==0) && (kohlenstoff==2 || kohlenstoff==3))
			{
				(*i)->FFType(10);
				ua.erase(i++);
				continue;
			}
			// planar nitrogen in a ring structure
			if(arom || (isRingAtom(*i) && (*i)->Geometry() == GEO(GEO::tri)))
			{
				(*i)->FFType(13);
				ua.erase(i++);
				continue;
			}
			// nitrogen as a hydrogen bond acceptor
			if(wasserstoff==0)
			{
				(*i)->FFType(11);
				ua.erase(i++);
				continue;
			}
			// nitrogen as a hydrogen bond donor
			if(wasserstoff>0)
			{
				(*i)->FFType(12);
				ua.erase(i++);
				continue;
			}
		}
		++i;
	}
	// search fluorine
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::F))
		{
			(*i)->FFType(30);
			ua.erase(i++);
			continue;
		}
		++i;
	}
	// search bromine
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::Br))
		{
			(*i)->FFType(28);
			ua.erase(i++);
			continue;
		}
		++i;
	}
	// search chlorine
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::Cl))
		{
			(*i)->FFType(27);
			ua.erase(i++);
			continue;
		}
		++i;
	}
	// search iron
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::Fe))
		{
			(*i)->FFType(26);
			ua.erase(i++);
			continue;
		}
		++i;
	}
	// search phosphorus
	i=ua.begin();
	while (i!=ua.end())
	{
		if((*i)->Element() == ELE(ELE::P))
		{
			(*i)->FFType(22);
			ua.erase(i++);
			continue;
		}
		++i;
	}
	
	// atoms for vdw interaction
// 	for (AtomList::iterator i=_latoms.begin(); i!=_latoms.end(); ++i)
	_lvdw.clear();
	_vdw_latom.clear();
	for (unsigned int i=0; i!=_latoms.size(); ++i)
	{
		unsigned int type(vdWType(_latoms[i]));
		if (type)
		{
			_lvdw.push_back(vdWType(_latoms[i]));
			_vdw_latom.push_back(i);
		}
	}
	// calculate vdw receptor grid for every ligand type
	for (unsigned int i=0; i!=_vdw_latom.size(); ++i)
	{
		if(_vdw_grid.iscalculated(_lvdw[i])) continue;
		boost::function<Float_t(Vec3_t)> calc = boost::lambda::bind(&PMF04::calcVDW, this, boost::lambda::_1, _lvdw[i]);
		_vdw_grid.calculate(_lvdw[i], calc );
	}
	// calculate statistical potential receptor grid for every ligand type
	for(AtomList::iterator i=_latoms.begin();i!= _latoms.end();i++)
	{
		if(_sp_grid.iscalculated((*i)->FFType())) continue;
		boost::function<Float_t(Vec3_t)> calc = boost::lambda::bind(&PMF04::calcSP, this, boost::lambda::_1, (*i)->FFType());
		_sp_grid.calculate((*i)->FFType(), calc);
	}
	// build neighborhood list for internal ligand vdw clashes penalty
	_internal_vdw_neighborhood_list=vector<vector<unsigned int> >(_vdw_latom.size(), vector<unsigned int>());
	_internal_vdw_neighborhood_list_d0=vector<vector<Float_t> >(_vdw_latom.size(), vector<Float_t>());
	for (unsigned int i=0; i!=_vdw_latom.size(); ++i)
	{
		AtomList tmp=One_N_Neighborhood(_latoms[_vdw_latom[i]], 5);
		for (AtomList::iterator j=tmp.begin(); j!=tmp.end(); ++j)
		{
			if ((*j)->Element()==ELE(ELE::H) || (*j)->Element()==ELE(ELE::LP)) continue;
			// find index of *j in _latoms; the atom number
			unsigned int index_of_j=distance(_latoms.begin(), 
				find(_latoms.begin(), _latoms.end(), *j));
			// find the atom number in _vdw_latom
			unsigned int k=distance(_vdw_latom.begin(), 
				find(_vdw_latom.begin(), _vdw_latom.end(), index_of_j));
			// if the index of the atom number is larger than i, its the first time that we see this pair
			if (k>i)
			{
				_internal_vdw_neighborhood_list[i].push_back(index_of_j);
				_internal_vdw_neighborhood_list_d0[i].push_back(vdW_type[_lvdw[i]]+vdW_type[_lvdw[k]]);
			}
		}
	}
}
Float_t PMF04::calcVDW(const Vec3_t& pos, unsigned int type)
{
	Float_t tmp=0,help;
	const vector<unsigned int> &list(_vdw_nh.getNeighborhood(pos));
	for(vector<unsigned int>::const_iterator j=list.begin(); j!=list.end(); ++j)
	{
		Float_t d=(pos-_patoms[_vdw_patom[*j]]->Coord()).length2();
		if(d<_vdw_range)
		{
			Float_t d0=(vdW_type[type]+_pvdw[*j]);
			d0*=d0;
			d0/=d;
			d0*=d0;
			help=d0*(d0-2);
			if (help < 0) continue;
			if (help >_vdw_repulsion_cutoff) help=_vdw_repulsion_cutoff;
			tmp+=help;
		}
	}
	return tmp;
}
Float_t PMF04::calcSP(const Vec3_t& pos, unsigned int type)
{
	Float_t tmp=0;
	bool lig_c=false,prot_c=false;
	if(( type > 0 && type < 9) || type == 29) lig_c=true;	
	const vector<unsigned int> &list(_sp_nh.getNeighborhood(pos));
	
	for(vector<unsigned int>::const_iterator j=list.begin(); j!=list.end(); ++j)
	{
		if(_patoms[*j]->FFType() > 0 && _patoms[*j]->FFType() < 7) prot_c=true;
		Float_t d=(pos-_patoms[*j]->Coord()).length2();
		if(lig_c && prot_c)
		{
			if(d>36) continue;
			tmp+=_parameter[_patoms[*j]->FFType()][type][_assignDistance[(int)ceil((d*25))]];
			continue;
		}
		if(d>81) continue;
		tmp+=_parameter[_patoms[*j]->FFType()][type][_assignDistance[(int)ceil((d*25))]];
	}
	return tmp;
}

unsigned int PMF04::vdWType(AtomKey a)
{
	if (a->Element()==ELE::H || a->Element()==ELE::LP)
	{
		return 0;
	}
	
	BondList bonds(aBonds(a));
	AtomList abonded(aAdjAtoms(a));
	
	if (a->Element()==ELE::C)
	{
		if (a->Geometry()==GEO::tet)
		{
			return 1;
		}
		if (a->Geometry()==GEO::tri)
		{
			bool ar=false;
			bool re=false;
			for (BondList::iterator j=bonds.begin(); j!=bonds.end(); ++j)
			{
				if ((*j)->Type()==BT::ar)
				{
					ar=true;
				}
				if ((*j)->Type()==BT::re)
				{
					re=true;
				}
			}
			if (ar)
			{
				return 2;
			}
			if (re)
			{
				return 3;
			}
			return 4;
		}
		if (a->Geometry()==GEO::lin)
		{
			return 5;
		}
	}
	if (a->Element()==ELE::N)
	{
		if (a->Geometry()==GEO::tet)
		{
			return 6;
		}
		if (a->Geometry()==GEO::tri)
		{
			bool ar=false;
			bool re=false;
			for (BondList::iterator j=bonds.begin(); j!=bonds.end(); ++j)
			{
				if ((*j)->Type()==BT::ar)
				{
					ar=true;
				}
				if ((*j)->Type()==BT::re)
				{
					re=true;
				}
			}
			if (ar)
			{
				return 7;
			}
			if (re)
			{
				return 8;
			}
			return 9;
		}
		if (a->Geometry()==GEO::lin)
		{
			return 10;
		}
	}
	if (a->Element()==ELE::O)
	{
		if (a->Geometry()==GEO::tet)
		{
			return 11;
		}
		if (a->Geometry()==GEO::tri)
		{
			bool ar=false;
			bool re=false;
			for (BondList::iterator j=bonds.begin(); j!=bonds.end(); ++j)
			{
				if ((*j)->Type()==BT::ar)
				{
					ar=true;
				}
				if ((*j)->Type()==BT::re)
				{
					re=true;
				}
			}
			if (ar)
			{
				return 12;
			}
			if (re)
			{
				return 13;
			}
			return 14;
		}
		if (a->Geometry()==GEO::lin)
		{
			return 15;
		}
	}
	if (a->Element()==ELE::S)
	{
		if (a->Geometry()==GEO::tet)
		{
			return 16;
		}
		if (a->Geometry()==GEO::tri)
		{
			bool ar=false;
			bool re=false;
			for (BondList::iterator j=bonds.begin(); j!=bonds.end(); ++j)
			{
				if ((*j)->Type()==BT::ar)
				{
					ar=true;
				}
				if ((*j)->Type()==BT::re)
				{
					re=true;
				}
			}
			if (ar)
			{
				return 17;
			}
			if (re)
			{
				return 18;
			}
			return 19;
		}
		if (a->Geometry()==GEO::lin)
		{
			return 20;
		}
	}
	if (a->Element()==ELE::P)
	{
		if (a->Geometry()==GEO::tet)
		{
			return 21;
		}
		if (a->Geometry()==GEO::tri)
		{
			bool ar=false;
			bool re=false;
			for (BondList::iterator j=bonds.begin(); j!=bonds.end(); ++j)
			{
				if ((*j)->Type()==BT::ar)
				{
					ar=true;
				}
				if ((*j)->Type()==BT::re)
				{
					re=true;
				}
			}
			if (ar)
			{
				return 22;
			}
			if (re)
			{
				return 23;
			}
			return 24;
		}
		if (a->Geometry()==GEO::lin)
		{
			return 25;
		}
	}
	
	if (a->Element()==ELE::F)
	{
		return 26;
	}
	if (a->Element()==ELE::Cl)
	{
		return 27;
	}
	if (a->Element()==ELE::Br)
	{
		return 28;
	}
	if (a->Element()==ELE::I)
	{
		return 29;
	}
	if (a->Element()==ELE::Si)
	{
		return 30;
	}
	if (a->Element()==ELE::K)
	{
		return 31;
	}
	if (a->Element()==ELE::U)
	{
		return 32;
	}
	if (a->Element()==ELE::Li)
	{
		return 33;
	}
	if (a->Element()==ELE::Na)
	{
		return 34;
	}
	if (a->Element()==ELE::Ca)
	{
		return 35;
	}
	if (a->Element()==ELE::Cd)
	{
		return 36;
	}
	if (a->Element()==ELE::Co)
	{
		return 37;
	}
	if (a->Element()==ELE::Cu)
	{
		return 38;
	}
		if (a->Element()==ELE::Fe)
	{
		return 39;
	}
	if (a->Element()==ELE::Hg)
	{
		return 40;
	}
	if (a->Element()==ELE::Mg)
	{
		return 41;
	}
	if (a->Element()==ELE::Mn)
	{
		return 42;
	}
	if (a->Element()==ELE::Ni)
	{
		return 43;
	}
	if (a->Element()==ELE::Zn)
	{
		return 44;
	}
	if (a->Element()==ELE::Al)
	{
		return 45;
	}
	
	throw FitnessError("no vdw parameter!");
	return 0;
}

AtomList PMF04::One_N_Neighborhood(AtomKey a, unsigned int n)
{
	AtomList result;
	// allready visited atoms
	map<AtomKey, bool> visited;
	// bfs search by using a queue
	queue<pair<AtomKey, unsigned int> > disc;
	disc.push(pair<AtomKey, unsigned int>(a,1));
	
	while (!disc.empty())
	{
		pair<AtomKey, unsigned int> front=disc.front();
		disc.pop();
		AtomList front_adj(aAdjAtoms(front.first));
		for (AtomList::iterator k=front_adj.begin(); k!=front_adj.end(); ++k)
		{
			if (!visited[*k])
			{
				if (front.second>=(n-1)) result.push_back(*k);
				disc.push(pair<AtomKey, unsigned int>(*k,front.second+1));
				visited[*k]=1;
			}
		}
	}
	return result;
}
void PMF04::eval(const Vec3_t& pos, const Quat_t& ori,const vector<Float_t>& angle, vector<Float_t>*& result)
{
 	vector<Vec3_t> ligandcoord(updateLigandCoord(pos, ori, angle));
	Float_t sp=0,tmp_vdw=0,tmp_sp=0;
	for(unsigned int i=0; i!=_vdw_latom.size(); ++i)
	{
		tmp_vdw+=_vdw_grid.interpolate(ligandcoord[_vdw_latom[i]], _lvdw[i]);
	}
	sp+=(_SF1*tmp_vdw);
	for(unsigned int i=0; i!=_latoms.size(); ++i)
	{
		tmp_sp+=_sp_grid.interpolate(ligandcoord[i], _latoms[i]->FFType());
	}
	sp+=(_SF2*tmp_sp);
	// internal ligand vdw clashes penalty
	Float_t ivdw=0;
	for (unsigned int i=0; i!=_internal_vdw_neighborhood_list.size(); ++i)
	{
		for (unsigned int j=0; j!=_internal_vdw_neighborhood_list[i].size(); ++j)
		{
			Float_t d=(ligandcoord[_vdw_latom[i]]-ligandcoord[_internal_vdw_neighborhood_list[i][j]]).length2();
			if(d<64)
			{
				Float_t d0=_internal_vdw_neighborhood_list_d0[i][j];
				d0*=d0;
				d0/=d;
				d0*=d0;
				Float_t clash=d0*(d0-2);
				if (clash>0) ivdw+=clash;
			}
		}
	}
	sp+=(_SF3*ivdw);
// 	cout <<"ivdw : "<<ivdw<<"\ttmp_vdw : "<<tmp_vdw<<"\ttmp_sp : "<<tmp_sp<<"\t"<<"sp :"<<sp<<endl;
	_result[0]=sp; //hier muss das ergebnis rein
	result=&_result;
 }
}
