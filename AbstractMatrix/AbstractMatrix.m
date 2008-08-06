(* Mathematica package *)

(* Created by the Wolfram Workbench 14-Jul-2008 *)

BeginPackage["AbstractMatrix`"]
(* Exported symbols added here with SymbolName::usage *) 

AbstractMatrix::usage = "AbstractMatrix is a package which provides abstract matrix operations such as dot product and trace."
AbstractDot::usage = "AbstractDot[X,Y] is the abstract matrix dot product of X and Y."
AbstractTrace::usage = "AbstractTrace[X] is the abstract matrix trace of the matrix X."
SimplifyTrace::usage = "SimplifyTrace[X] puts any occurance of AbstractTrace in X into normal form."

(* Error Messages *)
SimplifyTrace::notCanonical = "Warning, not necessarily in canonical form ``."

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

(* The dot product is distributive *)
e : AbstractDot[_, _Plus] := Distribute[Unevaluated[e]]
e : AbstractDot[_Plus, _] := Distribute[Unevaluated[e]]

(* Print AbstractDot as a space *)
Format[AbstractDot[x_, y_]] := Infix[AbstractDot[x, y],"\[CircleTimes]"];

(****************************** AbstractTrace ********************************)
(* AbstractTrace(aB)=a*AbstractTrace(B), a is a number, B is an AbstractMatrix *)
AbstractTrace[a_?NumericQ x_] := a*AbstractTrace[x];

(* AbstractTrace(a)=a, where a is a number *)
AbstractTrace[a_?NumericQ] := a;

(* Trace is distributive *)
e : AbstractTrace[_Plus] := Distribute[Unevaluated[e]];

(* Print AbstractTrace[x] as a tr[x] *)
Format[AbstractTrace[x_]] := "tr"x;

(*AbstractTrace[x_]:=AbstractTrace[Sort[x]]/;!OrderedQ[x]&&(Length[x]==2)*)(* \
AbstractTrace(AB)=AbstractTrace(BA), A and B maAbstractTraceices *)
(*AbstractTrace[x_]:=AbstractTrace[RotateLeft[x]]/;!OrderedQ[x]*)
(*AbstractTrace[x_]:=AbstractTrace[RotateLeft[x,Ordering[x,1]-1]]/;x[[Ordering[x,1]]][[1,1]]\
!=x[[{1}]][[1,1]]*)

(****************************** SimplifyTrace ********************************)
(* Extract just the numeric quantities x *)
SetAttributes[ExtractNumbers, Listable]
ExtractNumbers[x_] := Extract[x, Position[x, _Integer]]

(* Calculate all cyclic permutations of x *)
CyclicPermutations[x_] := With[{len = Length[x]},
  Table[RotateLeft[x, n], {n, 0, len - 1}]
]

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

(*SimplifyTrace[Times[a_,x_AbstractTrace]]:=Times[a,SimplifyTrace[x]]/;NumberQ[a]\
*)
End[]

EndPackage[]
