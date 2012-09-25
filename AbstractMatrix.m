(* ::Package:: *)

(* AbstractMatrix Mathematica package 
   This is a Mathematica for doing abstract matrix calculations. It
   treats any non-numeric symbols as matrices.
    
   Copyright 2009 Barry Wardell
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *)

BeginPackage["AbstractMatrix`"]

(* Exported functions *)
AbstractMatrix::usage = "AbstractMatrix is a package which provides abstract matrix operations such as dot product and trace."
AbstractDot::usage = "AbstractDot[X,Y] is the abstract matrix dot product of X and Y."
AbstractTrace::usage = "AbstractTrace[X] is the abstract matrix trace of the matrix X."
SimplifyTrace::usage = "SimplifyTrace[X] puts any occurance of AbstractTrace in X into normal form. This allows for the simplification of expressions with AbstractTrace[] in them."
CyclicPermutations::usage = "CyclicPermutations[x] calculates the cyclic permutations of X."
InsertCyclic::usage = "InsertCyclic[list, elem] produces a list of lists with each element corresponding to elem inserted into list at a different position."
ListInsert::usage = "ListInsert[list, elem, pos] produces a list of lists with element inserted into list at positions pos."
Contraction::usage = "Contraction[x, {i1, i2}] symbolizes a contraction of indices i1 and i2 together."

(* Error Messages *)
SimplifyTrace::notCanonical = "Warning, not necessarily in canonical form ``."

(* Is the symbol a bitensor or not? *)
BitensorQ::usage = "Tests if a symbol is a bitensor."
BiscalarQ::usage = "Tests if a symbol is a bitensor."
NotBitensorQ::usage = "Tests if a symbol is not a bitensor."

Begin["`Private`"]
(* Implementation of the package *)

(******************************* AbstractDot *********************************)
(* Dot product with a number is just normal multiplication*)
AbstractDot[a_?NumericQ, x_] := a x;
AbstractDot[x_, a_?NumericQ] := a x;

(* Dot products are associative but not distributive (Maeder pp. 178) 
 * FIXME: 
 * 1. Why does setting the OneIdentity Attribute make this an order of
 *    magnitude slower?
 * 2. Why does putting this before the first two definitions make this 
 *    several orders of magnitude slower? It doesn't seem to fully
 *    expand things out without it set. *)
SetAttributes[AbstractDot, {Flat, OneIdentity}];

(* Numbers pull through the dot product *)
AbstractDot[a_?NumericQ x_, y_] := a AbstractDot[x,y];
AbstractDot[x_, a_?NumericQ y_] := a AbstractDot[x,y];

AbstractDot[m_?NotBitensorQ^(a_),x_] := m^(a) AbstractDot[x];
AbstractDot[m_?NotBitensorQ^(a_) x_, y_] := m^(a) AbstractDot[x,y];
AbstractDot[x_, m_?NotBitensorQ^(a_) y_] := m^(a) AbstractDot[x,y];

(* Biscalars pull through dot product *)
AbstractDot[a_?BiscalarQ x_, y_] := a AbstractDot[x,y];
AbstractDot[x_, a_?BiscalarQ y_] := a AbstractDot[x,y];
AbstractDot[x_ a_?BiscalarQ, y_] := a AbstractDot[x,y];
AbstractDot[x_, y_ a_?BiscalarQ] := a AbstractDot[x,y];

AbstractDot[m_?BiscalarQ^(a_),x_] := m^(a) AbstractDot[x];
AbstractDot[m_?BiscalarQ^(a_) x_, y_] := m^(a) AbstractDot[x,y];
AbstractDot[x_, m_?BiscalarQ^(a_) y_] := m^(a) AbstractDot[x,y];
AbstractDot[x_ m_?BiscalarQ^(a_), y_] := m^(a) AbstractDot[x,y];
AbstractDot[x_, y_ m_?BiscalarQ^(a_)] := m^(a) AbstractDot[x,y];

(* The dot product is distributive *)
e : AbstractDot[_, _Plus] := Distribute[Unevaluated[e]]
e : AbstractDot[_Plus, _] := Distribute[Unevaluated[e]]

(* Print AbstractDot as a cross inside a circle *)
Format[AbstractDot[x_, y_]] := Infix[AbstractDot[x, y],"\[CircleTimes]"];

(******************************* Contraction *********************************)
(*Contraction[x_,{_,_}] := x;*)
(* Numbers pull through the contraction *)
Contraction[a_?NumericQ x_, l_List] := a Contraction[x, l];

(* The contraction is distributive *)
e : Contraction[_Plus, {_,_}] := Distribute[Unevaluated[e]];

Format[Contraction[x_, {i1_, i2_}]] :=



\!\(\*SubscriptBox["\"\<C\>\"", 
RowBox[{"{", 
RowBox[{"i1", ",", "i2"}], "}"}]]\)[x];

(****************************** AbstractTrace ********************************)
(* AbstractTrace (aB)=a*AbstractTrace (B), a is a number, B is an AbstractMatrix *)
AbstractTrace[a_?NumericQ x_] := a*AbstractTrace[x];

AbstractTrace[m_?NotBitensorQ^(a_)] := m^(a);
AbstractTrace[m_?NotBitensorQ^(a_) x_] := m^(a) AbstractTrace[x];

AbstractTrace[a_?BiscalarQ] := a;
AbstractTrace[a_?BiscalarQ x_] := a*AbstractTrace[x];
AbstractTrace[x_ a_?BiscalarQ] := a*AbstractTrace[x];
AbstractTrace[m_?BiscalarQ^(a_)] := m^(a);
AbstractTrace[m_?BiscalarQ^(a_) x_] := m^(a) AbstractTrace[x];
AbstractTrace[x_ m_?BiscalarQ^(a_)] := m^(a) AbstractTrace[x];

(* AbstractTrace (a)=a, where a is a number *)
AbstractTrace[a_?NumericQ] := a;

(* Trace is distributive *)
e : AbstractTrace[_Plus] := Distribute[Unevaluated[e]];

(* Print AbstractTrace[x] as a tr[x] *)
Format[AbstractTrace[x_]] := "tr"[x];

(****************************** SimplifyTrace ********************************)
(* Convert symbols representing matrices to a sequence of numbers in the range
 * [0,x-1] where x is the number of unique matrices present, eg.
 * AbstractDot[x[2],x[5],x[1],x[2]] -> {1,2,0,1}  }*)
SetAttributes[ExtractNumbers, Listable]
ExtractNumbers[x_] := Module[{seq = {Sequence @@ x}, elems},
	(* Get an ordered list of the unique elements *)
	elems = Union[seq];
	
	(* Replace our sequence of matrices with a sequence of numbers *)
	seq /. Thread[Sort[seq] -> Range[0, Length[seq] - 1]]
]

(* Calculate all cyclic permutations of x *)
CyclicPermutations[x_] := With[{len = Length[x]},
  Table[RotateLeft[x, n], {n, 0, len - 1}]
]

InsertCyclic[l_, s_] := ListInsert[l, s, Range[Length[l]+1]];

ListInsert[l_, s_, pos_] := Map[Insert[l, s, #]&, pos];

(* Calculate the "weight" of the permutation x *)
PermWeight[x_List] := With[{n = Length[x]},
  Sum[ n^(n-i) x[[i]], {i, 1, n}]
]

(* Perform simplification on the trace given that it is invariant
 * under cyclic permutations of the matrices *)
SimplifyTrace[x_AbstractTrace] := With[{matrix = First[x]},
	Module[{perms, permweights, min, canonicalPermPos},
		(* Calculate cyclic permutations *)
		perms = CyclicPermutations[matrix];
		
		(* Get the weight of all the permutations *)
		permweights = PermWeight /@ ExtractNumbers[perms];
		
		(* Find the minimum permutation weight *)		
		min = Min[permweights];

		(* Find the position in the list of permutations of the first
		 * occurance of the permutation with the minimum weight *)
		canonicalPermPos = First[Position[permweights, min]];
		
		(* Print an error if we're not getting a canonical form		
		With[{numequivperms = Count[permweights, Min[permweights]]},
			If[(numequivperms>1) && (numequivperms<Length[permweights]),
				Message[SimplifyTrace::notCanonical, x]
			];
		]*)

		(* Return the permutation with the minimum weight *)
		AbstractTrace[Extract[perms,canonicalPermPos]]
	]
]

SimplifyTrace[x_Plus] := Map[SimplifyTrace, x]

SimplifyTrace[x_Times] := Map[SimplifyTrace, x]

SimplifyTrace[a_] := a /; NumberQ[a]

SimplifyTrace[x_Power] := Power[SimplifyTrace[x[[1]]], x[[2]]]

End[]

EndPackage[]




