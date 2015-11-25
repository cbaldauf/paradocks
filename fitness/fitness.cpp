#include "fitness.hpp"
#include <queue>
#include <cmath>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <iostream>

using namespace boost::lambda;

namespace paradocks
{

unsigned int BFSDepth(AtomKey a)
{
	unsigned int result=0;
	// allready visited atoms
	map<AtomKey, bool> visited;
	// bfs search by using a queue
	queue<pair<AtomKey, unsigned int> > disc;
	disc.push(pair<AtomKey, unsigned int>(a,0));
	
	while (!disc.empty())
	{
		pair<AtomKey, unsigned int> front=disc.front();
		disc.pop();
		AtomList front_adj(aAdjAtoms(front.first));
		for (AtomList::iterator k=front_adj.begin(); k!=front_adj.end(); ++k)
		{
			if (!visited[*k])
			{
				disc.push(pair<AtomKey, unsigned int>(*k,front.second+1));
				visited[*k]=true;
				if (front.second>result) result=front.second;
			}
		}
	}
	return result;
}

void Fitness::initProtein(const ChainList& prot, const Vec3_t& as, Float_t r)
{
	// active site
	_center=as;
	_radius=r;
	
	// assign protein atoms	
	_patoms=clAtoms(prot);
	if (_patoms.empty()) throw FitnessError("Protein is empty!");
}

void Fitness::initLigand(const ChainList& lig, unsigned int h_handling)
{
	_latoms=clAtoms(lig);
	if (!_latoms.size()) throw FitnessError("Ligand is empty!");
	AtomList connected(aExtendConnected(*(_latoms.begin())));
	if (_latoms.size()!=connected.size()) throw FitnessError("Ligand is not connected!");
	
	// find the graph center fragment
	AtomKey center=*(min_element(_latoms.begin(), _latoms.end(), (bind(BFSDepth, _1) < bind(BFSDepth, _2))));
	AtomList centerfragment(RigidFragment(center, h_handling));
	
	// calculate offset
	_offset=Vec3_t();
	for (AtomList::iterator i=centerfragment.begin(); i!=centerfragment.end(); ++i)
	{
		_offset=_offset+(*i)->Coord();
	}
	_offset=_offset/centerfragment.size();
	
	// move the ligan d to the origin
	_lcoord.clear();
	for (AtomList::iterator i=_latoms.begin(); i!=_latoms.end(); ++i)
	{
		(*i)->Coord((*i)->Coord()-_offset);
		// save initial coordinates
		_lcoord.push_back((*i)->Coord());
	}
	
	BondList bonds=alBonds(_latoms);
	BondList rot;
	// find rotatable bonds
	for (BondList::iterator i=bonds.begin(); i!=bonds.end(); ++i)
	{
		if ((*i)->Type()==BT::s && !(BranchEnd((*i)->Atom1(), h_handling)) &&
		!(BranchEnd((*i)->Atom2(), h_handling)) && !isRingBond(*i))
		{
			rot.push_back(*i);
		}
	}
		
	// build up the fragments for bond rotation
	// axis of rotation is from -> to
	_allrotfragments.clear();
	for (BondList::iterator i=rot.begin(); i!=rot.end(); ++i)
	{
		rotfragment rf;
		int d1, d2;
		d1=GraphDistance(center,(*i)->Atom1());
		d2=GraphDistance(center,(*i)->Atom2());
		if (d1<d2)
		{
			vector<AtomKey>::iterator pos=find(_latoms.begin(), _latoms.end(), (*i)->Atom1());
			if (pos==_latoms.end())
			{
				throw FitnessError("index error in initLigand!");
			}
			rf.from=distance(_latoms.begin(), pos);
			rf.dist=d1;
			pos=find(_latoms.begin(), _latoms.end(), (*i)->Atom2());
			if (pos==_latoms.end())
			{
				throw FitnessError("index error in initLigand!");
			}
			rf.to=distance(_latoms.begin(), pos);
		}
		else
		{
			vector<AtomKey>::iterator pos=find(_latoms.begin(), _latoms.end(), (*i)->Atom2());
			if (pos==_latoms.end())
			{
				throw FitnessError("index error in initLigand!");
			}
			rf.from=distance(_latoms.begin(), pos);
			rf.dist=d2;
			pos=find(_latoms.begin(), _latoms.end(), (*i)->Atom1());
			if (pos==_latoms.end())
			{
				throw FitnessError("index error in initLigand!");
			}
			rf.to=distance(_latoms.begin(), pos);
		}
		vector<AtomKey> atoms=CollectBranch(_latoms[rf.from], _latoms[rf.to]);
		for (vector<AtomKey>::iterator j=atoms.begin(); j!=atoms.end(); ++j)
		{
			vector<AtomKey>::iterator pos=find(_latoms.begin(), _latoms.end(), *j);
			if (pos==_latoms.end())
			{
				throw FitnessError("index error in initLigand!");
			}
			unsigned int index= distance(_latoms.begin(), pos);
			rf.atomindex.push_back(index);
		}
		
		_allrotfragments.push_back(rf);
	}
	
	sort(_allrotfragments.begin(), _allrotfragments.end(), less_rotfragment());
	
	_rotcount=rot.size();
}

struct less_bond : public binary_function<boost::tuple<int, int, string>, boost::tuple<int, int, string> ,bool> {
bool operator()(boost::tuple<int, int, string> x, boost::tuple<int, int, string> y) 
{
	if (x.get<0>() == y.get<0>())
	{
		return x.get<1>() < y.get<1>();
	}
	else
	{
		return x.get<0>() < y.get<0>();
	}
}};

void Fitness::writeOptimum(const Vec3_t& pos, const Quat_t& ori, 
	const vector<Float_t>& angle, ostream& outfile)
{
	vector<Vec3_t> ligandcoord(updateLigandCoord(pos, ori, angle));
	
	vector<Float_t>* dummy;
	eval(pos, ori, angle, dummy);
	
	BondList lb=alBonds(_latoms);
	
	unsigned int latomssize=0, lbsize=0;
	for (unsigned int i=0; i!=_latoms.size(); ++i)
	{
		if (_latoms[i]->Element()!=ELE::LP)
		{
			++latomssize;
		}
	}
	for (BondList::iterator i=lb.begin(); i!=lb.end(); ++i)
	{
		if ((*i)->Atom1()->Element()!=ELE::LP && (*i)->Atom2()->Element()!=ELE::LP)
		{
			++lbsize;
		}
	}

	outfile<<"@<TRIPOS>MOLECULE" << endl;
	outfile<<(*_latoms.begin())->Parent()->Parent()->Name()<<endl;
	outfile<<latomssize<<" "<<lbsize<<endl;
	outfile<< endl << endl << endl <<endl;
	outfile<<"@<TRIPOS>ATOM"<< endl;
	
	map<AtomKey, unsigned int> atomtoid;
	unsigned int bid=1;
	for (unsigned int i=0; i!=_latoms.size(); ++i)
	{
		if (_latoms[i]->Element()!=ELE::LP)
		{
			outfile<<i+1<<" "<< _latoms[i]->Name()<<" "<<
			ligandcoord[i].x()<<" "<<
			ligandcoord[i].y()<<" "<<
			ligandcoord[i].z()<<" "<<
			_latoms[i]->MOL2Type() << endl;
			atomtoid[_latoms[i]]=i+1;
		}
	}
	outfile<<"@<TRIPOS>BOND"<< endl;
	vector <boost::tuple<int, int, string> > bondvector;
	for (BondList::iterator i=lb.begin(); i!=lb.end(); ++i)
	{
		if ((*i)->Atom1()->Element()!=ELE::LP && (*i)->Atom2()->Element()!=ELE::LP)
		{
			if (atomtoid[(*i)->Atom1()]<atomtoid[(*i)->Atom2()])
			{
				bondvector.push_back(boost::make_tuple(atomtoid[(*i)->Atom1()], atomtoid[(*i)->Atom2()], (*i)->MOL2Type()));
			}
			else
			{
				bondvector.push_back(boost::make_tuple(atomtoid[(*i)->Atom2()], atomtoid[(*i)->Atom1()], (*i)->MOL2Type()));
			}
		}
	}
	
	sort(bondvector.begin(), bondvector.end(), less_bond());
	
	for (vector <boost::tuple<int, int, string> >::iterator i=bondvector.begin(); i!=bondvector.end(); ++i)
	{
		outfile<<bid<<" "<<i->get<0>()<<" "<<i->get<1>()<<" "<<i->get<2>()<<endl;
		++bid;
	}
}

void Fitness::restoreLigand()
{
	for (unsigned int i=0; i!=_latoms.size(); ++i)
	{
		_latoms[i]->Coord(_lcoord[i]);
	}
}

AtomList Fitness::RigidFragment(AtomKey a, unsigned int h_handling)
{
	AtomList result;
	result.push_back(a);
	// allready visited atoms
	map<AtomKey, bool> visited;
	// bfs search by using a queue
	queue<AtomKey> disc;
	disc.push(a);
	
	while (!disc.empty())
	{
		AtomKey front=disc.front();
		disc.pop();
		AtomList front_adj(aAdjAtoms(front));
		for (AtomList::iterator k=front_adj.begin(); k!=front_adj.end(); ++k)
		{
			if (!visited[*k])
			{
				visited[*k]=true;
				BondKey b=*(aaFormBond(front, *k).begin());
				// the atom belongs to the rigid fragment if the new atom is connected by
				// a not rotatable bond and the atom is not a branch end
				if (b->Type()!=BT::s || isRingBond(b) || 
					BranchEnd(*k, h_handling) || BranchEnd(front, h_handling))
				{
					disc.push(*k);
					result.push_back(*k);
				}
			}
		}
	}
	return result;
}

bool Fitness::BranchEnd(AtomKey a, unsigned int h_handling)
{
	switch (h_handling)
	{
		case 0:
		{
			if (aHeavyCount(a)>1) return false;
			else return true;
			break;
		}
		case 1:
		{
			int heavy=aHeavyCount(a);
			if (heavy>1) return false;
			else if (((a)->Element()>ELE::H) && ((a)->Element()!=ELE::C))
			{
				if ((heavy+aHCount(a))>1) return false;
			}
			else return true;
			break;
		}
		case 2:
		{
			if ((aHeavyCount(a)+aHCount(a))>1) return false;
			else return true;
			break;
		}
		default:
		{
			#ifdef DEBUG
			throw FitnessError("Fitness::BranchEnd: wrong h_handling value.");
			#endif
		}
	}
	return false;
}

AtomList Fitness::CollectBranch(AtomKey from, AtomKey to)
{
	AtomList result;
	// allready visited atoms
	map<AtomKey, bool> visited;
	visited[from]=true;
	visited[to]=true;
	// bfs search by using a queue
	queue<AtomKey> disc;
	disc.push(to);
	
	while (!disc.empty())
	{
		AtomKey front=disc.front();
		disc.pop();
		AtomList front_adj(aAdjAtoms(front));
		for (AtomList::iterator k=front_adj.begin(); k!=front_adj.end(); ++k)
		{
			if (!visited[*k])
			{
				visited[*k]=true;
				disc.push(*k);
				result.push_back(*k);
			}
		}
	}
	return result;
}

vector<Vec3_t> Fitness::updateLigandCoord(const Vec3_t& pos, const Quat_t& ori, const vector<Float_t>& angle)
{
	vector<Vec3_t> result(_lcoord);
	// restore original ligand
	if (_allrotfragments.size()!=angle.size()) 
		throw FitnessError("Fitness::updateLigandCoord: _allrotfragments.size()!=angle.size()");
	for (unsigned int i=0 ; i!=_allrotfragments.size(); ++i)
	{
		Matrix_t t1, t2, final;
		t1.makeTranslate(-(result[_allrotfragments[i].from]));
		t2.makeTranslate(result[_allrotfragments[i].from]);
		Matrix_t rotm(Quat_t(angle[i], 
			result[_allrotfragments[i].to]-result[_allrotfragments[i].from]));
		final=(t1*rotm)*t2;
		for (vector<unsigned int>::iterator j=_allrotfragments[i].atomindex.begin(); j!=_allrotfragments[i].atomindex.end(); ++j)
		{
			result[*j]=result[*j]*final;
		}
	}
	Matrix_t transm;
	transm.makeTranslate(pos);
	Matrix_t centermatrix=Matrix_t(ori)*transm;
	for (unsigned int i=0; i!=_latoms.size(); ++i)
	{
		result[i]=result[i]*centermatrix;
	}
	return result;
}

void Grid::setupGrid(const Vec3_t& p1, const Vec3_t& p2, Float_t spacing, unsigned int typecount, Float_t outside)
{
	_spacing=spacing;
	_outside=outside;
	if (_spacing<0.1) throw FitnessError("Grid::Grid: Grid lattice spacing needs to be at least 0.1 A");
	
	Vec3_t dimension;
	if (p1.x()<p2.x())
	{
		_offset.x()=p1.x();
		dimension.x()=(p2.x()-p1.x());
	}
	else
	{
		_offset.x()=p2.x();
		dimension.x()=(p1.x()-p2.x());
	}
	
	if (p1.y()<p2.y())
	{
		_offset.y()=p1.y();
		dimension.y()=(p2.y()-p1.y());
	}
	else
	{
		_offset.y()=p2.y();
		dimension.y()=(p1.y()-p2.y());
	}
	
	if (p1.z()<p2.z())
	{
		_offset.z()=p1.z();
		dimension.z()=(p2.z()-p1.z());
	}
	else
	{
		_offset.z()=p2.z();
		dimension.z()=(p1.z()-p2.z());
	}
	_xsize=(unsigned int)ceil(dimension.x()/_spacing);
	_ysize=(unsigned int)ceil(dimension.y()/_spacing);
	_zsize=(unsigned int)ceil(dimension.z()/_spacing);
	if (_xsize==0) _xsize=1;
	if (_ysize==0) _ysize=1;
	if (_zsize==0) _zsize=1;
	_end.x()=_offset.x()+_xsize*_spacing;
	_end.y()=_offset.y()+_ysize*_spacing;
	_end.z()=_offset.z()+_zsize*_spacing;
		
	_grid=vector<vector<vector<vector<Float_t> > > >(typecount);
}

void Grid::calculate(unsigned int type, boost::function<Float_t(Vec3_t)>& calc)
{
	if (type>_grid.size()) 
		throw FitnessError("Grid::Calculate: max count of grid types exceeded");
	_grid[type]=vector<vector<vector<Float_t> > >(_xsize+1, vector<vector<Float_t> >
	(_ysize+1, vector<Float_t>(_zsize+1)));
	for(unsigned int x=0; x<_xsize+1; ++x)
	{
		for(unsigned int y=0; y<_ysize+1; ++y)
		{
			for(unsigned int z=0; z<_zsize+1; ++z)
			{
				_grid[type][x][y][z]=calc(_offset+Vec3_t(x*_spacing,y*_spacing,z*_spacing));
			}
		}
	}	
}

bool Grid::iscalculated(unsigned int type)
{
	return !(_grid[type].empty());
}

Float_t Grid::interpolate(const Vec3_t& p, unsigned int type)
{
	int x, y, z;
	Vec3_t index((p-_offset)/_spacing);
	// calculate index and delta
	Float_t tmp;
	Float_t dx(modf(index.x(), &tmp));
	x=(int)tmp;
	if (x<0 || x>=(int)_xsize) return _outside;
	Float_t dy(modf(index.y(), &tmp));
	y=(int)tmp;
	if (y<0 || y>=(int)_ysize) return _outside;
	Float_t dz(modf(index.z(), &tmp));
	z=(int)tmp;
	if (z<0 || z>=(int)_ysize) return _outside;
	Float_t x11=_grid[type][x][y][z];
	Float_t x12=_grid[type][x+1][y][z]-x11;
	Float_t x21=_grid[type][x][y][z+1];
	Float_t x22=_grid[type][x+1][y][z+1]-x21;
	Float_t x31=_grid[type][x][y+1][z];	
	Float_t x32=_grid[type][x+1][y+1][z]-x31;
	Float_t x41=_grid[type][x][y+1][z+1];
	Float_t x42=_grid[type][x+1][y+1][z+1];
	Float_t x32x12=x32-x12;
	Float_t x31x11=x31-x11;
	// put all together
	return (x11+x12*dx+(x31x11+(x32x12)*dx)*dy)+
	(x21-x11+(x22-x12)*dx+(x41-x21-x31x11+(x42-x41-x22-x32x12)*dx)*dy)*dz;
}

void StaticNeighborhood::setupList(const AtomList& al, Float_t cutoff, Float_t spacing)
{
	_spacing=spacing;
	
	if (al.empty()) throw FitnessError("StaticNeighborhood::setupList: AtomList empty!");
	if (cutoff<1) throw FitnessError("StaticNeighborhood::setupList: cutoff < 1!");
	if (spacing<0.1) throw FitnessError("StaticNeighborhood::setupList: spacing < 0.1!");
	
	// calculate the dimensions
	Float_t xmin, xmax, ymin, ymax, zmin, zmax;
	xmin=xmax=al[0]->Coord().x();
	ymin=ymax=al[0]->Coord().y();
	zmin=zmax=al[0]->Coord().z();
	for (unsigned int i=1; i!=al.size(); ++i)
	{
		if (al[i]->Coord().x()<xmin) xmin=al[i]->Coord().x();
		if (al[i]->Coord().y()<ymin) ymin=al[i]->Coord().y();
		if (al[i]->Coord().z()<zmin) zmin=al[i]->Coord().z();
		if (al[i]->Coord().x()>xmax) xmax=al[i]->Coord().x();
		if (al[i]->Coord().y()>ymax) ymax=al[i]->Coord().y();
		if (al[i]->Coord().z()>zmax) zmax=al[i]->Coord().z();
	}
	
	_offset=Vec3_t(xmin, ymin, zmin)-Vec3_t(cutoff, cutoff, cutoff)-Vec3_t(_spacing, _spacing, _spacing);
	_end=Vec3_t(xmax, ymax, zmax)+Vec3_t(cutoff, cutoff, cutoff)+Vec3_t(_spacing, _spacing, _spacing);
	Vec3_t dimension=_end-_offset;
	
	_xsize=(unsigned int)ceil(dimension.x()/_spacing);
	_ysize=(unsigned int)ceil(dimension.y()/_spacing);
	_zsize=(unsigned int)ceil(dimension.z()/_spacing);
	if (_xsize==0) _xsize=1;
	if (_ysize==0) _ysize=1;
	if (_zsize==0) _zsize=1;
	_end.x()=_offset.x()+_xsize*_spacing;
	_end.y()=_offset.y()+_ysize*_spacing;
	_end.z()=_offset.z()+_zsize*_spacing;
	
	// allocate neighborhood list
	_neighborhoodlist=vector<vector<vector<vector<unsigned int> > > >
	( _xsize, vector<vector<vector<unsigned int> > >
	( _ysize, vector<vector<unsigned int> >
	( _zsize ,vector<unsigned int>())));
	
	
	// all "boxes" in _boxes_below_cutoff are below the cutoff
	// and the proteinatoms should be put in the neighborhood list
	// of all boxes
	vector<boost::tuple<int, int, int> > below_cutoff;
	int boxsize=(int)ceil(cutoff/_spacing)+1;
	Float_t cutoff_square=cutoff*cutoff;
	for (int i=-boxsize; i<=boxsize; ++i)
	{
		for (int j=-boxsize; j<=boxsize; ++j)
		{
			for (int k=-boxsize; k<=boxsize; ++k)
			{
				Float_t xdiff=0, ydiff=0, zdiff=0;
				if (i) xdiff=(abs(i)-1)*_spacing;
				if (j) ydiff=(abs(j)-1)*_spacing;
				if (k) zdiff=(abs(k)-1)*_spacing;
				
				if ((xdiff*xdiff+ydiff*ydiff+zdiff*zdiff)<cutoff_square)
				{
					below_cutoff.push_back(boost::tuple<int, int, int>(i, j, k));
				}
			}
		}
	}

	// put every atom in the in the boxes of the neighborhood list
	for (unsigned int i=0; i!=al.size(); ++i)
	{
		Vec3_t index((al[i]->Coord()-_offset)/_spacing);
		unsigned int x=(unsigned int)floor(index.x());
		unsigned int y=(unsigned int)floor(index.y());
		unsigned int z=(unsigned int)floor(index.z());
		for (vector<boost::tuple<int, int, int> >::iterator j=below_cutoff.begin(); j!=below_cutoff.end(); ++j)
		{
			_neighborhoodlist[x+j->get<0>()][y+j->get<1>()][z+j->get<2>()].push_back(i);
		}
	}
}

const vector<unsigned int> StaticNeighborhood::getNeighborhood(const Vec3_t& p)
{
	Vec3_t index((p-_offset)/_spacing);
	int x=(int)floor(index.x());
	int y=(int)floor(index.y());
	int z=(int)floor(index.z());
	if (x<0 || y<0 || z<0 || x>=(int)_xsize || y>=(int)_ysize || z>=(int)_zsize)
	{
		return _empty;
	}
	return _neighborhoodlist[x][y][z];
}

}

