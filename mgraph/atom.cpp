#include "atom.hpp"

using namespace boost::assign;

namespace mgraph 
{

const string ELE::_convertele2str[107] = {"LP", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
"Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb",
"Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba",
"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
"Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
"Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
"No", "Lr", "Rf", "DU", "*"};

const map<string, ELE::_ELE> ELE::_convertstr2ele = map_list_of
("LP", ELE::LP) ("H",  ELE::H)  ("He", ELE::He) ("Li", ELE::Li) ("Be", ELE::Be)
("B",  ELE::B)  ("C",  ELE::C)  ("N",  ELE::N)  ("O",  ELE::O)  ("F",  ELE::F) 
("Ne", ELE::Ne) ("Na", ELE::Na) ("Mg", ELE::Mg) ("Al", ELE::Al) ("Si", ELE::Si)
("P",  ELE::P)  ("S",  ELE::S)  ("Cl", ELE::Cl) ("Ar", ELE::Ar) ("K",  ELE::K) 
("Ca", ELE::Ca) ("Sc", ELE::Sc) ("Ti", ELE::Ti) ("V",  ELE::V)  ("Cr", ELE::Cr)
("Mn", ELE::Mn) ("Fe", ELE::Fe) ("Co", ELE::Co) ("Ni", ELE::Ni) ("Cu", ELE::Cu)
("Zn", ELE::Zn) ("Ga", ELE::Ga) ("Ge", ELE::Ge) ("As", ELE::As) ("Se", ELE::Se)
("Br", ELE::Br) ("Kr", ELE::Kr) ("Rb", ELE::Rb) ("Sr", ELE::Sr) ("Y",  ELE::Y)
("Zr", ELE::Zr) ("Nb", ELE::Nb) ("Mo", ELE::Mo) ("Tc", ELE::Tc) ("Ru", ELE::Ru)
("Rh", ELE::Rh) ("Pd", ELE::Pd) ("Ag", ELE::Ag) ("Cd", ELE::Cd) ("In", ELE::In)
("Sn", ELE::Sn) ("Sb", ELE::Sb) ("Te", ELE::Te) ("I",  ELE::I)  ("Xe", ELE::Xe)
("Cs", ELE::Cs) ("Ba", ELE::Ba) ("La", ELE::La) ("Ce", ELE::Ce) ("Pr", ELE::Pr)
("Nd", ELE::Nd) ("Pm", ELE::Pm) ("Sm", ELE::Sm) ("Eu", ELE::Eu) ("Gd", ELE::Gd)
("Tb", ELE::Tb) ("Dy", ELE::Dy) ("Ho", ELE::Ho) ("Er", ELE::Er) ("Tm", ELE::Tm)
("Yb", ELE::Yb) ("Lu", ELE::Lu) ("Hf", ELE::Hf) ("Ta", ELE::Ta) ("W",  ELE::W)
("Re", ELE::Re) ("Os", ELE::Os) ("Ir", ELE::Ir) ("Pt", ELE::Pt) ("Au", ELE::Au)
("Hg", ELE::Hg) ("Tl", ELE::Tl) ("Pb", ELE::Pb) ("Bi", ELE::Bi) ("Po", ELE::Po)
("At", ELE::At) ("Rn", ELE::Rn) ("Fr", ELE::Fr) ("Ra", ELE::Ra) ("Ac", ELE::Ac)
("Th", ELE::Th) ("Pa", ELE::Pa) ("U",  ELE::U)  ("Np", ELE::Np) ("Pu", ELE::Pu)
("Am", ELE::Am) ("Cm", ELE::Cm) ("Bk", ELE::Bk) ("Cf", ELE::Cf) ("Es", ELE::Es)
("Fm", ELE::Fm) ("Md", ELE::Md) ("No", ELE::No) ("Lr", ELE::Lr) ("Rf", ELE::Rf)
("DU", ELE::DU) ("*", ELE::Any);

const string GEO::_convertg2str[7] = {"none", "lin", "tri", "tet", "bip", "oct", "UNK"};

const map<string, GEO::_GEO> GEO::_convertstr2g = map_list_of
("none", GEO::none) ("lin", GEO::lin)
("tri", GEO::tri) ("tet", GEO::tet)
("bip", GEO::bip) ("oct", GEO::oct)
("UNK", GEO::UNK);

}

