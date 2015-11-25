#include "residue.hpp"

using namespace boost::assign;

namespace mgraph {

const string RT::_convertrt2str[34] = {"UNK", "ALA", "ARG", "ASN", "ASP", "CYS",
"GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
"THR", "TRP", "TYR", "VAL", "ASZ", "GLZ", "HIE", "HIP", "LYZ", "CYM", "CYX",
"TYM", "AMN", "AMI", "CXL", "CXC"};

const map<string, RT::_RT> RT::_convertstr2rt = map_list_of
("UNK", RT::UNK)
("A", RT::ALA) ("ALA", RT::ALA) ("R", RT::ARG) ("ARG", RT::ARG)
("N", RT::ASN) ("ASN", RT::ASN) ("D", RT::ASP) ("ASP", RT::ASP)
("C", RT::CYS) ("CYS", RT::CYS) ("Q", RT::GLN) ("GLN", RT::GLN)
("E", RT::GLU) ("GLU", RT::GLU) ("G", RT::GLY) ("GLY", RT::GLY)
("H", RT::HIS) ("HIS", RT::HIS) ("I", RT::ILE) ("ILE", RT::ILE)
("L", RT::LEU) ("LEU", RT::LEU) ("K", RT::LYS) ("LYS", RT::LYS)
("M", RT::MET) ("MET", RT::MET) ("F", RT::PHE) ("PHE", RT::PHE)
("P", RT::PRO) ("PRO", RT::PRO) ("S", RT::SER) ("SER", RT::SER)
("T", RT::THR) ("THR", RT::THR) ("W", RT::TRP) ("TRP", RT::TRP)
("Y", RT::TYR) ("TYR", RT::TYR) ("V", RT::VAL) ("VAL", RT::VAL)
("ASZ", RT::ASZ) ("ASH", RT::ASZ)
("GLZ", RT::GLZ) ("GLH", RT::GLZ)
("HID", RT::HIS)
("HIE", RT::HIE)
("HIP", RT::HIP)
("LYZ", RT::LYZ) ("LYN", RT::LYZ)
("CYM", RT::CYM)
("CYX", RT::CYX)
("TYM", RT::TYM)
("AMN", RT::AMN)
("AMI", RT::AMI)
("CXL", RT::CXL)
("CXC", RT::CXC)
;

}
