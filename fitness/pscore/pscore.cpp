#include "pscore.hpp"
#define HBCUTOFF 4.7

using namespace mgraph;

namespace paradocks
{

PScore::PScore(vector<string>& parameter) : Fitness(parameter), _result(_termscount)
{
	if (parameter.size()!=0)
	{
		string emsg("PScore expects no parameter");
		throw FitnessError(emsg);
	}
	_hbfaktor=-14.5847;
	_hbdist=0.910731;
	_hbangle=DegreesToRadians(178.962);
	_inv_hbdist=1/_hbdist;
	_inv_hbangle=1/_hbangle;
}

Float_t PScore::vdW_type[VDWTYPECOUNT]=
{
1,         // H 0
2.36682,   // C tet 1
2.15347,   // C tri ar 2
1.88582,   // C tri re 3
2.23108,   // C tri 4
1.85249,   // C lin 5
1.77486,   // N tet 6
1.54659,   // N tri ar 7
1.83274,   // N tri re 8
1.75951,   // N tri 9
1.54278,   // N lin 10
1.79145,   // O tet 11
1.93169,   // O tri ar 12
1.69397,   // O tri re 13
1.74588,   // O tri 14
2.62296,   // S tet 15
1.91157,   // S tri ar 16
2,         // S tri re 17
2.03912,   // S tri 18
2,         // S lin 19
2.00623,   // P tet 20
2,         // P tri ar 21
2,         // P tri re 22
2,         // P tri 23
2,         // P lin 24
1.70554,   // F 25
1.71643,   // Cl 26
1.76647,   // Br 27
2.05,      // I 28
2,         // Si 29
1.25,      // metall K 30
1.25,      // metall U 31
1.25,      // metall Li 32
1.25,      // metall Na 33
1.18868,   // metall Ca 34
1.25,      // metall Cd 35
1.25,      // metall Co 36
1.25,      // metall Cu 37
0.793814,  // metall Fe 38
1.25,      // metall Hg 39
0.925381,  // metall Mg 40
1.06108,   // metall Mn 41
1.25,      // metall Ni 42
0.9429,    // metall Zn 43
1.25,      // metall Al 44
};

void PScore::initProtein(const ChainList& prot, const Vec3_t& as, Float_t r)
{
	CorrectGraphByTopology(prot);
	// lone pairs are needed
	AddLP(prot);
	
	Fitness::initProtein(prot, as, r);
	
	// atoms for vdw interaction
	AtomList patoms_tmp;
	for (unsigned int i=0; i!=_patoms.size(); ++i)
	{
		if (_patoms[i]->Coord().x()<_center.x()-_radius-8 || _patoms[i]->Coord().y()<_center.y()-_radius-8 || _patoms[i]->Coord().z()<_center.z()-_radius-8 ||
			_patoms[i]->Coord().x()>_center.x()+_radius+8 || _patoms[i]->Coord().y()>_center.y()+_radius+8 || _patoms[i]->Coord().z()>_center.z()+_radius+8)
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
	// set up neighboorhood list
	_vdw_nh.setupList(patoms_tmp, 8, 1);
	// set up the grid
	_vdw_grid.setupGrid(as-Vec3_t(r,r,r), as+Vec3_t(r,r,r), 0.2, 46, 10);
	
	
	//atoms for h-bonding
	patoms_tmp.clear();
	for (unsigned int i=0; i!=_patoms.size(); ++i)
	{
		if (_patoms[i]->Coord().x()<_center.x()-_radius-HBCUTOFF || _patoms[i]->Coord().y()<_center.y()-_radius-HBCUTOFF || _patoms[i]->Coord().z()<_center.z()-_radius-HBCUTOFF ||
			_patoms[i]->Coord().x()>_center.x()+_radius+HBCUTOFF || _patoms[i]->Coord().y()>_center.y()+_radius+HBCUTOFF || _patoms[i]->Coord().z()>_center.z()+_radius+HBCUTOFF)
		{
			continue;
		}
		DA da;
		bool type(daType(_patoms[i], da, _patoms));
		if (type)
		{	
			_da_patom.push_back(i);
			_pda.push_back(da);
			_da_pradius.push_back(vdW_type[vdWType(_patoms[i])]);
			patoms_tmp.push_back(_patoms[i]);
		}
	}
	// set up neighboorhood list
	_da_nh.setupList(patoms_tmp, HBCUTOFF, 1);
}

void PScore::initLigand(const ChainList& lig)
{
	CorrectGraphByTopology(lig);
	// lone pairs are needed
	AddLP(lig);
	
	Fitness::initLigand(lig, 0);
		
	// atoms for vdw interaction
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
	
	// calculate receptor grid for every ligand type
	for (unsigned int i=0; i!=_vdw_latom.size(); ++i)
	{
		if(_vdw_grid.iscalculated(_lvdw[i])) continue;
		boost::function<Float_t(Vec3_t)> calc=
			boost::lambda::bind(&PScore::calcvdw, this, boost::lambda::_1, _lvdw[i]);
		_vdw_grid.calculate(_lvdw[i], calc );
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
			#ifdef DEBUG2
			if (index_of_j>=_latom.size()) throw FitnessError("PScore::initLigand: index larger than vector in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
			#endif
			// find the atom number in _vdw_latom
			unsigned int k=distance(_vdw_latom.begin(), 
				find(_vdw_latom.begin(), _vdw_latom.end(), index_of_j));
			#ifdef DEBUG2
			if (k>=_vdw_latom.size()) throw FitnessError("PScore::initLigand: index larger than vector in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
			#endif
			// if the index of the atom number is larger than i, its the first time that we see this pair
			if (k>i)
			{
				_internal_vdw_neighborhood_list[i].push_back(index_of_j);
				_internal_vdw_neighborhood_list_d0[i].push_back(vdW_type[_lvdw[i]]+vdW_type[_lvdw[k]]);
			}
		}
	}
	
	//atoms for h-bonding
	_da_latom.clear();
	_lda.clear();
	_da_lradius.clear();
	for (unsigned int i=0; i!=_latoms.size(); ++i)
	{
		DA da;
		bool type(daType(_latoms[i], da, _latoms));
		if (type)
		{	
			_da_latom.push_back(i);
			_lda.push_back(da);
			_da_lradius.push_back(vdW_type[vdWType(_latoms[i])]);
		}
	}
	
	log << "initial ligand position: " << _offset.x()<< " "<<_offset.y()<< " "<<_offset.z()<< " "<<endl;
	log << "rotatable bonds "<< _rotcount << endl;
}

Float_t PScore::calcvdw(const Vec3_t& pos, unsigned int type)
{
	Float_t tmp=0;
	const vector<unsigned int> &list(_vdw_nh.getNeighborhood(pos));
	for(vector<unsigned int>::const_iterator j=list.begin(); j!=list.end(); ++j)
	{
		Float_t d=(pos-_patoms[_vdw_patom[*j]]->Coord()).length2();
		if(d<64)
		{
			Float_t d0=vdW_type[type]+_pvdw[*j];
			d0*=d0;
			d0/=d;
			d0*=d0;
			tmp+=d0*(d0-2);
		}
	}
	if (tmp>15) tmp=15;
	return tmp;
}

void PScore::eval(const Vec3_t& pos, const Quat_t& ori, const vector<Float_t>& angle, vector<Float_t>*& result)
{
	vector<Vec3_t> ligandcoord(updateLigandCoord(pos, ori, angle));
	// protein-ligand vdW
	Float_t vdw=0;
	for(unsigned int i=0; i!=_vdw_latom.size(); ++i)
	{
		vdw+=_vdw_grid.interpolate(ligandcoord[_vdw_latom[i]], _lvdw[i]);
	}
	
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
				ivdw+=d0*(d0-2);
			}
		}
	}
	if (ivdw>0) vdw+=ivdw;
	_result[0]=vdw;
	
	// protein ligand h-bonding
	// the range for the decreasing benefit of hbonds
	vector<HBond> allhbond;
	for (unsigned int l=0; l!=_da_latom.size(); ++l)
	{
		const vector<unsigned int> &nhlist(_da_nh.getNeighborhood(ligandcoord[_da_latom[l]]));
		for(unsigned int p=0; p!=nhlist.size(); ++p)
		{
			Float_t d=(ligandcoord[_da_latom[l]]-_patoms[_da_patom[nhlist[p]]]->Coord()).length();
			Float_t d0=_da_lradius[l]+_da_pradius[nhlist[p]];
			// if distance is greater d0, f1 is zero 
			if (d>=d0) continue;
			HBond tmp;
			// between (d0-0.7) and d0 f1 goes from 1 to 0
			tmp.f1=_inv_hbdist*(d0-d);
			if (tmp.f1>1) tmp.f1=1;
			
			//ligand is acceptor, protein is donor
			if (_lda[l].a && _pda[nhlist[p]].d)
			{
				tmp.donor=(_pda[nhlist[p]]);
				tmp.acceptor=(_lda[l]);
				tmp.f2=0;
				tmp.f3=0;
				
				// score for DR-D-A
				if (tmp.donor.state==0)
				{
					// frozen donor D-H-A=180 degree
					// need to check all hydrogen
					for (unsigned int h=0; h!=tmp.donor.h.size(); ++h)
					{
						Vec3_t v1=_patoms[tmp.donor.h[h]]->Coord()-_patoms[tmp.donor.root]->Coord();
						v1.normalize();
						Vec3_t v2=_patoms[tmp.donor.h[h]]->Coord()-ligandcoord[tmp.acceptor.root];
						v2.normalize();
						Float_t angle=acos(v1*v2);
						Float_t f2=_inv_hbangle*(_hbangle-fabs(PI-angle));
						if (f2>tmp.f2) tmp.f2=f2;
					}
					if (tmp.f2<=0) continue;
				}
				if (tmp.donor.state==1)
				{
					// rot donor RD-D-A=109,5 degree
					Vec3_t v1=_patoms[tmp.donor.h[0]]->Coord()-_patoms[tmp.donor.root]->Coord();
					v1.normalize();
					Vec3_t v2=ligandcoord[tmp.acceptor.root]-_patoms[tmp.donor.root]->Coord();
					v2.normalize();
					Float_t angle=acos(v1*v2);
					tmp.f2=_inv_hbangle*(_hbangle-fabs(DegreesToRadians(109.5)-angle));
					if (tmp.f2<=0) continue;	
				}
				if (tmp.donor.state==2)
				{
					tmp.f2=1;
				}
				
				// score for AR-A-D
				if (tmp.acceptor.state==0)
				{
					// frozen acceptor A-LP-D=180 degree
					// need to check all LP
					for (unsigned int lp=0; lp!=tmp.acceptor.lp.size(); ++lp)
					{
						Vec3_t v1=ligandcoord[tmp.acceptor.lp[lp]]-ligandcoord[tmp.acceptor.root];
						v1.normalize();
						Vec3_t v2=ligandcoord[tmp.acceptor.lp[lp]]-_patoms[tmp.donor.root]->Coord();
						v2.normalize();
						Float_t angle=acos(v1*v2);
						Float_t f3=_inv_hbangle*(_hbangle-fabs(PI-angle));
						if (f3>tmp.f3) tmp.f3=f3;
					}
					if (tmp.f3<=0) continue;
				}
				if (tmp.acceptor.state==1)
				{
					// rot donor AR-A-D=109,5 degree
					Vec3_t v1=ligandcoord[tmp.acceptor.lp[0]]-ligandcoord[tmp.acceptor.root];
					v1.normalize();
					Vec3_t v2=_patoms[tmp.donor.root]->Coord()-ligandcoord[tmp.acceptor.root];
					v2.normalize();
					Float_t angle=acos(v1*v2);
					tmp.f3=_inv_hbangle*(_hbangle-fabs(DegreesToRadians(109.5)-angle));
					if (tmp.f3<=0) continue;	
				}
				if (tmp.acceptor.state==2)
				{
					tmp.f3=1;
				}
				tmp.f=tmp.f1*tmp.f2*tmp.f3;
				allhbond.push_back(tmp);
			}
			
			//ligand is donor, protein is acceptor
			if (_lda[l].d && _pda[nhlist[p]].a)
			{
				tmp.donor=(_lda[l]);
				tmp.acceptor=(_pda[nhlist[p]]);
				tmp.f2=0;
				tmp.f3=0;
				
				// score for DR-D-A
				if (tmp.donor.state==0)
				{
					// frozen donor D-H-A=180 degree
					// need to check all hydrogen
					for (unsigned int h=0; h!=tmp.donor.h.size(); ++h)
					{
						Vec3_t v1=ligandcoord[tmp.donor.h[h]]-ligandcoord[tmp.donor.root];
						v1.normalize();
						Vec3_t v2=ligandcoord[tmp.donor.h[h]]-_patoms[tmp.acceptor.root]->Coord();
						v2.normalize();
						Float_t angle=acos(v1*v2);
						Float_t f2=_inv_hbangle*(_hbangle-fabs(PI-angle));
						if (f2>tmp.f2) tmp.f2=f2;
					}
					if (tmp.f2<=0) continue;
				}
				if (tmp.donor.state==1)
				{
					// rot donor RD-D-A=109,5 degree
					Vec3_t v1=ligandcoord[tmp.donor.h[0]]-ligandcoord[tmp.donor.root];
					v1.normalize();
					Vec3_t v2=_patoms[tmp.acceptor.root]->Coord()-ligandcoord[tmp.donor.root];
					v2.normalize();
					Float_t angle=acos(v1*v2);
					tmp.f2=_inv_hbangle*(_hbangle-fabs(DegreesToRadians(109.5)-angle));
					if (tmp.f2<=0) continue;	
				}
				if (tmp.donor.state==2)
				{
					tmp.f2=1;
				}
				
				// score for AR-A-D
				if (tmp.acceptor.state==0)
				{
					// frozen acceptor A-LP-D=180 degree
					// need to check all LP
					for (unsigned int lp=0; lp!=tmp.acceptor.lp.size(); ++lp)
					{
						Vec3_t v1=_patoms[tmp.acceptor.lp[lp]]->Coord()-_patoms[tmp.acceptor.root]->Coord();
						v1.normalize();
						Vec3_t v2=_patoms[tmp.acceptor.lp[lp]]->Coord()-ligandcoord[tmp.donor.root];
						v2.normalize();
						Float_t angle=acos(v1*v2);
						Float_t f3=_inv_hbangle*(_hbangle-fabs(PI-angle));
						if (f3>tmp.f3) tmp.f3=f3;
					}
					if (tmp.f3<=0) continue;
				}
				if (tmp.acceptor.state==1)
				{
					// rot donor AR-A-D=109,5 degree
					Vec3_t v1=_patoms[tmp.acceptor.lp[0]]->Coord()-_patoms[tmp.acceptor.root]->Coord();
					v1.normalize();
					Vec3_t v2=ligandcoord[tmp.donor.root]-_patoms[tmp.acceptor.root]->Coord();
					v2.normalize();
					Float_t angle=acos(v1*v2);
					// diff to optimal
					tmp.f3=_inv_hbangle*(_hbangle-fabs(DegreesToRadians(109.5)-angle));
					if (tmp.f3<=0) continue;	
				}
				if (tmp.acceptor.state==2)
				{
					tmp.f3=1;
				}
				tmp.f=tmp.f1*tmp.f2*tmp.f3;
				allhbond.push_back(tmp);
			}
		}
	}
	Float_t hbscore=0;
	for (vector<HBond>::iterator i=allhbond.begin(); i!=allhbond.end(); ++i)
	{
		hbscore+=i->f;
	}
	hbscore*=_hbfaktor;
	_result[1]=hbscore;
	
	result=&_result;
}

/*
Float_t PScore::vdW_type[45]=
{
1.0, // H 0
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
2.11339, // S tet 15
1.96368, // S tri ar 16
2.0, // S tri re 17
2.0, // S tri 18
2.0, // S lin 19
1.95237, // P tet 20
2.0, // P tri ar 21
2.0, // P tri re 22
2.0, // P tri 23
2.0, // P lin 24
1.4997, // F 25
1.71409, // Cl 26
1.83567, // Br 27
2.05, // I 28
2.0, // Si 29
1.25, // metall K 30
1.25, // metall U 31
1.25, // metall Li 32
1.25, // metall Na 33
1.20589, // metall Ca 34
1.25, // metall Cd 35
1.25, // metall Co 36
1.25, // metall Cu 37
1.2024, // metall Fe 38
1.25, // metall Hg 39
1.21929, // metall Mg 40
1.25, // metall Mn 41
1.25, // metall Ni 42
1.20855, // metall Zn 43
1.25, // metall Al 44
};*/



unsigned int PScore::vdWType(AtomKey a)
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
	}
	if (a->Element()==ELE::S)
	{
		if (a->Geometry()==GEO::tet)
		{
			return 15;
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
				return 16;
			}
			if (re)
			{
				return 17;
			}
			return 18;
		}
		if (a->Geometry()==GEO::lin)
		{
			return 19;
		}
	}
	if (a->Element()==ELE::P)
	{
		if (a->Geometry()==GEO::tet)
		{
			return 20;
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
				return 21;
			}
			if (re)
			{
				return 22;
			}
			return 23;
		}
		if (a->Geometry()==GEO::lin)
		{
			return 24;
		}
	}
	
	if (a->Element()==ELE::F)
	{
		return 25;
	}
	if (a->Element()==ELE::Cl)
	{
		return 26;
	}
	if (a->Element()==ELE::Br)
	{
		return 27;
	}
	if (a->Element()==ELE::I)
	{
		return 28;
	}
	if (a->Element()==ELE::Si)
	{
		return 29;
	}
	if (a->Element()==ELE::K)
	{
		return 30;
	}
	if (a->Element()==ELE::U)
	{
		return 31;
	}
	if (a->Element()==ELE::Li)
	{
		return 32;
	}
	if (a->Element()==ELE::Na)
	{
		return 33;
	}
	if (a->Element()==ELE::Ca)
	{
		return 34;
	}
	if (a->Element()==ELE::Cd)
	{
		return 35;
	}
	if (a->Element()==ELE::Co)
	{
		return 36;
	}
	if (a->Element()==ELE::Cu)
	{
		return 37;
	}
		if (a->Element()==ELE::Fe)
	{
		return 38;
	}
	if (a->Element()==ELE::Hg)
	{
		return 39;
	}
	if (a->Element()==ELE::Mg)
	{
		return 40;
	}
	if (a->Element()==ELE::Mn)
	{
		return 41;
	}
	if (a->Element()==ELE::Ni)
	{
		return 42;
	}
	if (a->Element()==ELE::Zn)
	{
		return 43;
	}
	if (a->Element()==ELE::Al)
	{
		return 44;
	}
	
	throw FitnessError("no vdw parameter!");
	return 0;
}

bool PScore::daType(AtomKey a, DA& type, AtomList& al)
{
	unsigned int aindex=distance(al.begin(), find(al.begin(), al.end(), a));
	if (a->Element()==ELE::O || a->Element()==ELE::N || a->Element()==ELE::S)
	{
		// make frozen
		DA da;
		da.root=aindex;
		da.state=0;
		da.d=0;
		da.a=0;
		AtomList abonded(aAdjAtoms(a));
		vector<unsigned int> root_heavy;
		vector<unsigned int> root_h;
		vector<unsigned int> root_lp;
		for (AtomList::iterator j=abonded.begin(); j!=abonded.end(); ++j)
		{
			unsigned int index=distance(al.begin(), find(al.begin(), al.end(), *j));
			if ((*j)->Element()==ELE::LP)
			{
				root_lp.push_back(index);
			}
			else if ((*j)->Element()==ELE::H)
			{
				root_h.push_back(index);
			}
			else
			{
				root_heavy.push_back(index);
			}
		}
		if (!root_heavy.size())
		{
			// water
			if (al[da.root]->Element()==ELE::O && root_h.size()==2)
			{
				da.d=2;
				da.a=2;
				da.state=2;
				type=da;
				return 1;
			}
			throw FitnessError("missing D/A root in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));	
		}
		if (root_heavy.size()==1)
		{
			BondList rootbondlist=aaFormBond(al[root_heavy[0]], al[da.root]);
			if (rootbondlist.size()!=1) throw FitnessError("rootbondlist.size()!=1 in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
			if (rootbondlist[0]->Type()==BT::s) da.state=1;
			else da.state=0;
		}
		else
		{
			da.state=0;
		}
		da.d=root_h.size();
		if (da.d)
		{
			// find donor with respect to atoms
			if (da.state==0)
			{
				// thats easy, vector goes from root to hydrogen
				da.h=root_h;
			}
			else
			{
				// we do not save the hydrogen, but the heavy atom neighbor of the rootatom
				da.h.push_back(root_heavy[0]);
			}
		}
		
		da.a=root_lp.size();
		if (da.a)
		{
			// find acceptor with respect to atoms
			if (da.state==0)
			{
				//vector goes from root to lp
				da.lp=root_lp;
			}
			else
			{
				// we do not save the lp, but the heavy atom neighbor of the rootatom
				da.lp.push_back(root_heavy[0]);
			}
		}
		
		if (da.a || da.d)
		{
			type=da;
			return 1;
		}
	}
	if (a->Element()==ELE::Zn)
	{
		DA da;
		da.root=aindex;
		da.state=2;
		da.d=4;
		da.a=0;
		type=da;
		return 1;
	}
	return 0;
}

AtomList PScore::One_N_Neighborhood(AtomKey a, unsigned int n)
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

void PScore::AddLP(const ChainList& cl)
{
	AtomList atoms=clAtoms(cl);
	for (unsigned int i=0; i!=atoms.size(); ++i)
	{
		if (atoms[i]->Element()==ELE::O || atoms[i]->Element()==ELE::N || atoms[i]->Element()==ELE::S)
		{
			AtomKey root=atoms[i];
			AtomList abonded(aAdjAtoms(root));
			AtomList root_heavy_neighbour;
			AtomList root_hydrogen;
			for (AtomList::iterator j=abonded.begin(); j!=abonded.end(); ++j)
			{
				if ((*j)->Element()!=ELE::H && (*j)->Element()!=ELE::LP)
				{
					root_heavy_neighbour.push_back(*j);
				}
				if ((*j)->Element()==ELE::H)
				{
					root_hydrogen.push_back(*j);
				}
			}
			
			// construct lps for hbonds
			int n=root_heavy_neighbour.size()+root_hydrogen.size();
			if (root->Geometry()==GEO::tet)
			{
				if (n==3)
				{
					//construct one lp
					Vec3_t lp;
					for (AtomList::iterator j=root_hydrogen.begin(); j!=root_hydrogen.end(); ++j)
					{
						Vec3_t tmp=(*j)->Coord()-root->Coord();
						tmp.normalize();
						lp+=tmp;
					}
					for (AtomList::iterator j=root_heavy_neighbour.begin(); j!=root_heavy_neighbour.end(); ++j)
					{
						Vec3_t tmp=(*j)->Coord()-root->Coord();
						tmp.normalize();
						lp+=tmp;
					}
					lp*=-1;
					lp.normalize();
					AtomKey lpa=aCreate(root->Parent(), ELE(ELE::LP), root->Coord()-lp, GEO(GEO::none));
					BondKey lpb=bCreate(root, lpa, BT(BT::s));
					lpa->Name("tet3");
					lpa->MOL2Type("LP");
					lpb->MOL2Type("1");
				}
				if (n==2)
				{
					//construct two lp
					vector<Vec3_t> root_neighbour_vector;
					for (AtomList::iterator j=root_hydrogen.begin(); j!=root_hydrogen.end(); ++j)
					{
						Vec3_t tmp=(*j)->Coord()-root->Coord();
						tmp.normalize();
						root_neighbour_vector.push_back(tmp);
					}
					for (AtomList::iterator j=root_heavy_neighbour.begin(); j!=root_heavy_neighbour.end(); ++j)
					{
						Vec3_t tmp=(*j)->Coord()-root->Coord();
						tmp.normalize();
						root_neighbour_vector.push_back(tmp);
					}
					if (root_neighbour_vector.size()!=2) throw FitnessError("can not compute lp with root_neighbour_vector.size()!=2 in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
					// vector between the lps
					Vec3_t tmp1=root_neighbour_vector[0]+root_neighbour_vector[1];
					tmp1*=-1;
					tmp1.normalize();
					// vector orthogonal two root_neighbour_vector
					Vec3_t tmp2=root_neighbour_vector[0]^root_neighbour_vector[1];
					tmp2.normalize();
					// rotation
					Quat_t rot;
					rot.makeRotate (tmp1, tmp2);
					Quat_t slerp1;
					slerp1.slerp((54.75/90.0), Quat_t(), rot);
					Matrix_t mat1(slerp1);
					slerp1.slerp((-54.75/90.0), Quat_t(), rot);
					Matrix_t mat2(slerp1);
					
					Vec3_t tmp3=tmp1*mat1;
					Vec3_t tmp4=tmp1*mat2;
					
					AtomKey lpa1=aCreate(root->Parent(), ELE(ELE::LP), root->Coord()+tmp3, GEO(GEO::none));
					BondKey lpb1=bCreate(root, lpa1, BT(BT::s));
					lpa1->Name("tet2");
					lpa1->MOL2Type("LP");
					lpb1->MOL2Type("1");
					
					AtomKey lpa2=aCreate(root->Parent(), ELE(ELE::LP), root->Coord()+tmp4, GEO(GEO::none));
					BondKey lpb2=bCreate(root, lpa2, BT(BT::s));
					lpa2->Name("tet2");
					lpa2->MOL2Type("LP");
					lpb2->MOL2Type("1");
				}
				if (n==1)
				{
					//construct three lp
					vector<Vec3_t> root_neighbour_vector;
					for (AtomList::iterator j=root_hydrogen.begin(); j!=root_hydrogen.end(); ++j)
					{
						Vec3_t tmp=(*j)->Coord()-root->Coord();
						tmp.normalize();
						root_neighbour_vector.push_back(tmp);
					}
					for (AtomList::iterator j=root_heavy_neighbour.begin(); j!=root_heavy_neighbour.end(); ++j)
					{
						Vec3_t tmp=(*j)->Coord()-root->Coord();
						tmp.normalize();
						root_neighbour_vector.push_back(tmp);
					}
					if (root_neighbour_vector.size()!=1) throw FitnessError("can not compute lp with root_neighbour_vector.size()!=1 in FILE:"+string(__FILE__)+" LINE:"+boost::lexical_cast<string>(__LINE__));
					// TODO make this right
					AtomKey lpa1=aCreate(root->Parent(), ELE(ELE::LP), root->Coord()-root_neighbour_vector[0], GEO(GEO::none));
					BondKey lpb1=bCreate(root, lpa1, BT(BT::s));
					lpa1->Name("tet1");
					lpa1->MOL2Type("LP");
					lpb1->MOL2Type("1");
					
					AtomKey lpa2=aCreate(root->Parent(), ELE(ELE::LP), root->Coord()-root_neighbour_vector[0], GEO(GEO::none));
					BondKey lpb2=bCreate(root, lpa2, BT(BT::s));
					lpa2->Name("tet1");
					lpa2->MOL2Type("LP");
					lpb2->MOL2Type("1");
					
					AtomKey lpa3=aCreate(root->Parent(), ELE(ELE::LP), root->Coord()-root_neighbour_vector[0], GEO(GEO::none));
					BondKey lpb3=bCreate(root, lpa3, BT(BT::s));
					lpa3->Name("tet1");
					lpa3->MOL2Type("LP");
					lpb3->MOL2Type("1");
				}
			}
			if (root->Geometry()==GEO::tri)
			{
				if (n==2)
				{
					//construct one lp
					Vec3_t lp;
					for (AtomList::iterator j=root_hydrogen.begin(); j!=root_hydrogen.end(); ++j)
					{
						Vec3_t tmp=(*j)->Coord()-root->Coord();
						tmp.normalize();
						lp+=tmp;
					}
					for (AtomList::iterator j=root_heavy_neighbour.begin(); j!=root_heavy_neighbour.end(); ++j)
					{
						Vec3_t tmp=(*j)->Coord()-root->Coord();
						tmp.normalize();
						lp+=tmp;
					}
					lp*=-1;
					lp.normalize();
					AtomKey lpa1=aCreate(root->Parent(), ELE(ELE::LP), root->Coord()+lp, GEO(GEO::none));
					BondKey lpb1=bCreate(root, lpa1, BT(BT::s));
					lpa1->Name("tri2");
					lpa1->MOL2Type("LP");
					lpb1->MOL2Type("1");
				}
				if (n==1)
				{
					//construct two lp
					Vec3_t tmp1=root_heavy_neighbour[0]->Coord()-root->Coord();
					AtomList planeatoms=aAdjAtoms(root_heavy_neighbour[0]);
					Vec3_t tmp2(1,1,1);
					for (AtomList::iterator j=planeatoms.begin(); j!=planeatoms.end(); ++j)
					{
						if ((*j)!=root)
						{
							tmp2=(*j)->Coord()-root->Coord();
						}
					}
					Vec3_t tmp3=tmp1^tmp2;
					Vec3_t tmp4=tmp1^tmp3;
					Quat_t rot;
					rot.makeRotate (tmp1, tmp4);
					Quat_t slerp1;
					slerp1.slerp((120.0/90.0), Quat_t(), rot);
					Matrix_t mat1(slerp1);
					Vec3_t tmp5=tmp1*mat1;
					tmp5.normalize();
					Vec3_t tmp6=(tmp1*mat1)*mat1;
					tmp6.normalize();
					
					AtomKey lpa1=aCreate(root->Parent(), ELE(ELE::LP), root->Coord()+tmp5, GEO(GEO::none));
					BondKey lpb1=bCreate(root, lpa1, BT(BT::s));
					lpa1->Name("tri1");
					lpa1->MOL2Type("LP");
					lpb1->MOL2Type("1");
					
					AtomKey lpa2=aCreate(root->Parent(), ELE(ELE::LP), root->Coord()+tmp6, GEO(GEO::none));
					BondKey lpb2=bCreate(root, lpa2, BT(BT::s));
					lpa2->Name("tri1");
					lpa2->MOL2Type("LP");
					lpb2->MOL2Type("1");
				}
			}
		}
	}
}

}
