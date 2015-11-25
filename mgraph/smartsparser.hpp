#ifndef RESIDUETYPE_H
#define RESIDUETYPE_H

#include <string>

/**
\defgroup smarts SMARTS parser
	This SMARTS parser is implemented using Spirit an object-oriented recursive-descent parser generator framework. The grammar is derived from the Daylight Theory Manual (Version 4.9). The following grammar explanations are written in Extended Backus–Naur Form (EBNF).\n \n
	The following features differ from the original Daylight specification.\n
	<ol>
	 <li>isotopic specifications such as carbon-12, carbon-13 are omitted </li>
	 <li>disconnected compounds are omitted</li>
	 <li>chirality is omitted</li>
	 <li>all atoms, bonds, charges and attached hydrogens have to be explicitely specified</li>
	 <li>an additional <b>return_sign<b> "§" was introduced, which specifies atoms whose atom keys are returned after a subgraph match
	</ol>
	
**/

/**
\brief mySMARTS parser
\n
<b>any_atom ::=</b> "*" \n
<b>atom_symbol ::=</b> any_atom | "H" | "He" | "Li" | "Be" | "B" | "C" | "N" | "O" | "F" | "Ne" | "Na" | "Mg" | "Al" | "Si" | "P" | "S" | "Cl" | "Ar" | "K" | "Ca" | "Sc" | "Ti" | "V" | "Cr" | "Mn" | "Fe" | "Co" | "Ni" | "Cu" | "Zn" | "Ga" | "Ge" | "As" | "Se" | "Br" | "Kr" | "Rb" | "Sr" | "Y" | "Zr" | "Nb" | "Mo" | "Tc" | "Ru" | "Rh" | "Pd" | "Ag" | "Cd" | "In" | "Sn" | "Sb" | "Te" | "I" | "Xe" | "Cs" | "Ba" | "La" | "Ce" | "Pr" | "Nd" | "Pm" | "Sm" | "Eu" | "Gd" | "Tb" | "Dy" | "Ho" | "Er" | "Tm" | "Yb" | "Lu" | "Hf" | "Ta" | "W" | "Re" | "Os" | "Ir" | "Pt" | "Au" | "Hg" | "Tl" | "Pb" | "Bi" | "Po" | "At" | "Rn" | "Fr" | "Ra" | "Ac" | "Th" | "Pa" | "U" | "Np" | "Pu" | "Am" | "Cm" | "Bk" | "Cf" | "Es" | "Fm" | "Md" | "No" | "Lr" | "Rf" | "Db" | "Sg" | "Bh" | "Hs" | "Mt" | "Ds" | "Rg" \n
<b>atom_number ::=</b> "#1" | "#2" | ... | "#110" | "#111" \n
<b>atom ::=</b> open_square_bracket ( ([negation | return_sign] atom_symbol) | ([negation | return_sign] atom_number)) [atom_property] close_square_bracket \n
<b>open_square_bracket ::=</b> "[" \n
<b>close_square_bracket ::=</b> "]" \n
<b>negation ::=</b> "!" \n
<b>open_branch ::=</b> "(" \n 
<b>close_branch ::=</b> ")" \n
<b>bond ::=</b> [negation] "-" | [negation] "=" | [negation] "#" | [negation] ":" | "~" \n
<b>digit ::=</b> "0" | "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9" \n
<b>atom_property ::=</b> charge | attached_hydrogens | (charge attached_hydrogens) | (attached_hydrogens charge)\n 
<b>charge ::=</b> "+"[digit] | "-"[digit] \n
<b>attached_hydrogens ::=</b> "H"[digit] \n
<b>simple_branch ::=</b> open_branch bond atom[digit] (bond atom[digit])* close_branch \n
<b>multiple_branch ::=</b> open_branch bond atom[digit] (simple_branch | multiple_branch | bond atom[digit])* close_branch \n
<b>molecule ::=</b> atom[digit] ( multiple_branch | simple_branch | (bond atom [digit]) )* \n
<b>operator ::=</b> "or" \n
<b>query ::=</b> molecule (operator molecule)* \n
<b>return_sign ::=</b> "§"
\ingroup smarts	
**/



bool parse_smarts(char const* str);
// ChainList mySmartsQuery(string s);
void mySmartsQuery(char const* str);

#endif