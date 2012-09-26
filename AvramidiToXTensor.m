(* ::Package:: *)

(* AvramidiToXTensor Mathematica package 
   This is a Mathematica package for converting expressions from Avramidi notation
   to a the notation used by xTensor.
    
   Copyright 2009-2011 Barry Wardell
   
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

BeginPackage["CovariantSeries`AvramidiToXTensor`",
 {
  "xAct`xTensor`",
  "xAct`xPerm`",
  "CovariantSeries`",
  "CovariantSeries`AbstractMatrix`"
 }
]

AvramidiToXTensor::usage = 
  "AvramidiToXTensor is a package to convert expressions from the \
Avramidi notation (with \!\(\*SubscriptBox[\"\[ScriptCapitalK]\", \"n\
\"]\)) to Riemann tensors and their derivatives contracted with \
derivatives of \[Sigma] in xTensor form.";

RiemannPart::usage = 
  "RiemannPart[\!\(\*SubscriptBox[\"\[ScriptCapitalK]\", \"n\"]\)] \
gives the Riemann tensor part of \!\(\*SubscriptBox[\"s\", \"n\"]\) \
in xTensor form.";

$PrePrint=ScreenDollarIndices;
$CovDFormat="Postfix";
$DefInfoQ=False;
FreeIndices=IndexRange[\[Alpha],\[Nu]];
SigmaIndices=RotateLeft[IndexRange[a,l]]; (* RotateLeft so that 'a' is used for generating new indices *)
DefManifold[M, 4, Join[FreeIndices,SigmaIndices]];
DefMetric[-1, metric[-a, -b], CD, {";", "\[Del]"}, CurvatureRelations -> False];
(*DefTensor[\[Sigma][a], M];*)
DefTensor[PotentialP[], M];
DefTensor[HadamardW[], M];
DefTensor[Symmetric\[ScriptCapitalB][], M];
Symmetric\[ScriptCapitalB][ args__] := Symmetric\[ScriptCapitalB]@@Sort[{args}] /; ! OrderedQ[{args}];
IndexVBundle=VBundleOfIndex[a];

PrintAs[metric] ^= "g";
PrintAs[epsilonmetric] ^= "\[Epsilon]";
PrintAs[RiemannCD] ^= "R";
PrintAs[RicciCD] ^= "R";
PrintAs[RicciScalarCD] ^= "R";
PrintAs[WeylCD] ^= "W";
PrintAs[TFRicciCD] ^= "S";
PrintAs[PotentialP] ^= "\[ScriptCapitalP]";
PrintAs[HadamardW] ^= "W";
PrintAs[Symmetric\[ScriptCapitalB]] ^= "\[ScriptCapitalB]";

NumFreeIndices::usage="";


Begin["`Private`"] (* Begin Private Context *) 
(* Calculate the number of indices (contracted with sigma^a) used by an expression *)
NumSigmaIndices[x_] := Module[{positions, numIndices},
  (* Find where the K[n]'s are*)
  positions = Position[x, \[ScriptCapitalK][_]];
  
  (* Pull out the values for n and add them together *)
  numIndices = Plus @@ (Part[x, Sequence @@ #, 1] & /@ positions);
  
  (* Look for AddFreeIndices *)
  numIndices = numIndices - NumAddFreeIndices[x];
  
  (* Look for R[n] *)
  positions = Position[x, \[ScriptCapitalR][_]];
  
  (* Pull out the values for n and add them together *)
  numIndices = numIndices + Plus @@ (Part[x, Sequence @@ #, 1] & /@ positions);
 
  (* Look for \[ScriptCapitalP][n] *)
  positions = Position[x, \[ScriptCapitalP][_]];
  
  (* Pull out the values for n and add them together *)
  numIndices = numIndices + Plus @@ (Part[x, Sequence @@ #, 1] & /@ positions);

  (* Look for W[n] *)
  positions = Position[x, W[_]];

  (* Look for \[ScriptCapitalB][n] *)
  positions = Position[x, \[ScriptCapitalB][_]];

  (* Pull out the values for n and add them together *)
  numIndices = numIndices + Plus @@ (Part[x, Sequence @@ #, 1] & /@ positions);

  (* Look for g *)
  positions = Position[x, g];

  (* 2 extra indices per g *)
  numIndices = numIndices + 2 Length[positions];

  (* Look for Power *)
  positions = Position[x, Power[_,_]];

  numIndices = numIndices + (Plus @@ ((Part[x, Sequence @@ #, 2]-1)NumSigmaIndices[Part[x, Sequence @@ #, 1]] & /@ positions));

  numIndices
]

NumAddFreeIndices[x_] := Module[{positions, numIndices},
  (* Look for AddFreeIndices *)
  positions = Position[x, AddFreeIndex[_,_]];

  (* Pull out the values for n and add them together *)
  numIndices = (Plus @@ (Part[x, Sequence @@ #, 2] & /@ positions));

  numIndices
]

NumFreeIndices[\[ScriptCapitalR][_]] := 3;
NumFreeIndices[\[ScriptCapitalK][_]] := 2;
NumFreeIndices[\[ScriptCapitalP][_]] := 0;
NumFreeIndices[W[_]] := 0;
NumFreeIndices[\[ScriptCapitalB][_]] := 0;
NumFreeIndices[g] := 0;
NumFreeIndices[x_AbstractDot] := Plus@@(NumFreeIndices/@(List@@x)) - 2 (Length[x]-1);
NumFreeIndices[x_AbstractTrace]:= NumFreeIndices[x[[1]]]-2;
NumFreeIndices[x_Contraction]:= NumFreeIndices[x[[1]]]-2;
NumFreeIndices[x_Times] := Plus @@(NumFreeIndices/@(List@@x));
NumFreeIndices[x_Power] := x[[2]]NumFreeIndices[x[[1]]];
NumFreeIndices[x_AddFreeIndex] := x[[2]]+NumFreeIndices[x[[1]]];
NumFreeIndices[n_?NumericQ] := 0;

PartitionIndices[indicesIn_, partition_]:= Module[{iter, expr, indices=indicesIn},
  expr = {};
    For[iter = 1, iter <= Length[partition], iter++,
      AppendTo[expr, Take[indices, partition[[iter]]]];
      indices = Drop[indices, partition[[iter]]];
    ];
  expr
]


(* Convert K_n to n derivatives of Riemann in xTensor form *)
RiemannPart[\[ScriptCapitalK][num_], a_?AIndexQ, b_?AIndexQ, indices_IndexList] := Module[{vbundle, CD, expr, iter},
  (* Get the vbundle corresponding to the index a *)
  vbundle = IndexVBundle;

  (* Get the covariant derivative 
     FIXME: Should we worry about getting a different CD for each index? 
     Probably not, but maybe there is a case where it would be important. *)
  CD = CovDOfMetric[First[MetricsOfVBundle[vbundle]]];
  
  (* First, create the Riemann tensor *)
  expr = Riemann[CD][a, -indices[[1]], b, -indices[[2]]];
  
  (* Add covariant derivatives *)
  For[iter = 3, iter <= num, iter++,
   expr = CD[-indices[[iter]]]@expr;
  ];
  
  expr
]

(* Convert R_n to n derivatives of Ricci in xTensor form *)
RiemannPart[\[ScriptCapitalR][num_], a_?AIndexQ, b_?AIndexQ, c_?AIndexQ, indices_IndexList] := Module[{vbundle, CD, expr, iter},
  (* Get the vbundle corresponding to the index a *)
  vbundle = IndexVBundle;

  (* Get the covariant derivative
     FIXME: Should we worry about getting a different CD for each index?
     Probably not, but maybe there is a case where it would be important. *)
  CD = CovDOfMetric[First[MetricsOfVBundle[vbundle]]];

  (* First, create the Riemann tensor *)
  expr = Riemann[CD][a, b, -indices[[1]], c];

  (* Add covariant derivatives *)
  For[iter = 2, iter <= num, iter++,
   expr = CD[-indices[[iter]]]@expr;
  ];

  expr
]

(* In a sum, we treat each term independently *)
expr : AvramidiToXTensor[_Plus, _, _] := Distribute[Unevaluated[expr]]
expr : AvramidiToXTensor[_Plus, _] := Distribute[Unevaluated[expr]]


AvramidiToXTensor[a_?NumericQ, _] := a
AvramidiToXTensor[a_?NumericQ x_, y_] := a AvramidiToXTensor[x,y]

AvramidiToXTensor[CovariantSeries`m^(a_), _] := CovariantSeries`m^(a)
AvramidiToXTensor[CovariantSeries`m^(a_) x_, y_] := CovariantSeries`m^(a) AvramidiToXTensor[x,y]

AvramidiToXTensor[x_, vbundle_?VBundleQ] :=
  Module[{expr, sigmaIndices, freeIndices, nsi, nfi},
  (* Find how many indices we need *)
  nsi = NumSigmaIndices[x];

  (* Find out how many free indices we need *)
  nfi = NumFreeIndices[x];(* nfi = Length[FreeIndices]; (*FIXME*)*)

  (* Get some free indices *)
  (*freeIndices = FreeIndices[[1;;nfi]]; *)(*GetIndicesOfVBundle[vbundle, 2];*) (* FIXME: this could be different from 2 *)
  (*freeIndices[[2]] = -freeIndices[[2]];*)
  freeIndices = GetIndicesOfVBundle[vbundle, nfi, SigmaIndices];

  (* Get enough indices *)
  sigmaIndices = GetIndicesOfVBundle[vbundle, nsi, FreeIndices];

  expr = AvramidiToXTensor[x, IndexList@@freeIndices, IndexList@@sigmaIndices];
  
  expr
]

(* Multiplication is a bit tricky since we want all terms to have unique indices *)
AvramidiToXTensor[x_Times, freeIndices_IndexList, sigmaIndices_IndexList] :=
  Module[{parts, indicesPerTerm, termIndices, freeIndicesPerTerm, termFreeIndices},
  (* Separate multiplication into a list of each term *)
  parts = List @@ x;
  
  (* Figure out how many indices each term uses *)
  indicesPerTerm = Map[NumSigmaIndices, parts];
  freeIndicesPerTerm = Map[NumFreeIndices, parts];

  (* And divide the indices up between each term *)
  termIndices = PartitionIndices[sigmaIndices, indicesPerTerm];
  termFreeIndices = PartitionIndices[freeIndices, freeIndicesPerTerm];
  
  Times@@MapThread[AvramidiToXTensor[#1, #2, #3] &, {parts, List@@termFreeIndices, List @@ termIndices}]
]

(* Power *)
AvramidiToXTensor[Power[x_AbstractTrace, pow_?Positive], freeIndices_IndexList, sigmaIndices_IndexList] :=
  Module[{expr, partitionedIndices, partitionedFreeIndices, numFreeIndices, iter},
  partitionedIndices = Partition[sigmaIndices, NumSigmaIndices[x]];

  numFreeIndices = NumFreeIndices[x];
  partitionedFreeIndices = 
    If[numFreeIndices>0,
      Partition[freeIndices, numFreeIndices],
      ConstantArray[IndexList[], pow]];

  expr = 1;

  For[iter = 1, iter <= pow, iter++,
    expr = expr AvramidiToXTensor[x, partitionedFreeIndices[[iter]], partitionedIndices[[iter]]];
  ];

  expr
]

(* AbstractDot *)
AvramidiToXTensor[x_AbstractDot, freeIndices_IndexList, sigmaIndices_IndexList] :=
  Module[{vbundle, numTerms, numFreeIndices, numContractedIndices, contractedIndices, iter, expr, sigmaIndicesPerTerm, indicesUsed, parts, freeIndicesPerTerm, termFreeIndices},
  (* Get the vbundle corresponding to the index a *)
  vbundle = IndexVBundle;
  numTerms = Length[x];
  numContractedIndices = numTerms - 1;
  numFreeIndices = NumFreeIndices[x];

  (* Figure out how many free indices each term uses *)
  parts = List @@ x;
  freeIndicesPerTerm = Map[NumFreeIndices, parts]-{1,Sequence@@ConstantArray[2,Length[x]-2],1};

  (* And divide the indices up between each term *)
  termFreeIndices = PartitionIndices[freeIndices, freeIndicesPerTerm];

  contractedIndices = Table[DummyIn[vbundle],{numContractedIndices}];
  
  sigmaIndicesPerTerm = Map[NumSigmaIndices, List @@ x];

  expr = AvramidiToXTensor[x[[1]], IndexList[ Sequence@@termFreeIndices[[1]], contractedIndices[[1]] ], sigmaIndices[[ 1 ;; sigmaIndicesPerTerm[[1]] ]] ];
  indicesUsed = sigmaIndicesPerTerm[[1]];
  
  For[iter=2 , iter<numTerms, iter++,
  	expr = expr AvramidiToXTensor[x[[iter]], IndexList[ contractedIndices[[iter-1]], Sequence@@termFreeIndices[[iter]], contractedIndices[[iter]] ], sigmaIndices[[ indicesUsed+1 ;; indicesUsed+1+sigmaIndicesPerTerm[[iter]] ]] ];
  	indicesUsed += sigmaIndicesPerTerm[[iter]];
  ];
  
  expr = expr AvramidiToXTensor[x[[-1]], IndexList[contractedIndices[[-1]], Sequence@@termFreeIndices[[-1]]], sigmaIndices[[ indicesUsed + 1;; -1 ]] ]
]

(* AbstractTrace *)
AvramidiToXTensor[AbstractTrace[x_], freeIndices_IndexList, sigmaIndices_IndexList] := Module[{vbundle, a},
  vbundle = IndexVBundle;
  a = DummyIn[vbundle];
  AvramidiToXTensor[x, IndexList[a, a], sigmaIndices]
]

(* Contraction *)
AvramidiToXTensor[Contraction[x_, pos_List], freeIndices_IndexList, sigmaIndices_IndexList] :=
  Module[{vbundle, a, spos = Sort[pos], freeIndicesList},
  vbundle = IndexVBundle;
  a = DummyIn[vbundle];
  
  freeIndicesList = IndexList@@Riffle[List@@freeIndices, a, {spos[[1]], spos[[2]], spos[[2]]-spos[[1]]}];
  freeIndicesList[[pos[[1]]]] = - freeIndicesList[[pos[[1]]]];
  AvramidiToXTensor[x, freeIndicesList, sigmaIndices]
]

(* \[ScriptCapitalK] and \[ScriptCapitalR] and \[ScriptCapitalP] *)
AvramidiToXTensor[\[ScriptCapitalK][n_], inds_IndexList, sigmaIndices_IndexList] := RiemannPart[\[ScriptCapitalK][n], inds[[1]], -inds[[2]], sigmaIndices]
AvramidiToXTensor[\[ScriptCapitalR][n_], inds_IndexList, sigmaIndices_IndexList] := RiemannPart[\[ScriptCapitalR][n], -inds[[1]], -inds[[2]], -inds[[3]], sigmaIndices]
AvramidiToXTensor[\[ScriptCapitalP][0], inds_IndexList, sigmaIndices_IndexList] := PotentialP[]
AvramidiToXTensor[W[0], inds_IndexList, sigmaIndices_IndexList] := HadamardW[]
AvramidiToXTensor[\[ScriptCapitalB][0], inds_IndexList, sigmaIndices_IndexList] := Symmetric\[ScriptCapitalB][]
AvramidiToXTensor[\[ScriptCapitalP][0]^n_, inds_IndexList, sigmaIndices_IndexList] := PotentialP[]^n
AvramidiToXTensor[\[ScriptCapitalP][n_], inds_IndexList, sigmaIndices_IndexList] := PotentialP[Sequence@@(Times[-1,#]&/@sigmaIndices)]
AvramidiToXTensor[W[n_], inds_IndexList, sigmaIndices_IndexList] := HadamardW[Sequence@@(Times[-1,#]&/@sigmaIndices)]
AvramidiToXTensor[\[ScriptCapitalB][n_], inds_IndexList, sigmaIndices_IndexList] := Symmetric\[ScriptCapitalB][Sequence@@(Times[-1,#]&/@sigmaIndices)]
AvramidiToXTensor[g, inds_IndexList, sigmaIndices_IndexList] := metric[Sequence@@(Times[-1,#]&/@sigmaIndices)]

AvramidiToXTensor[\[ScriptCapitalP][k_]^pow_, IndexList[], sigmaIndices_IndexList] :=
  Module[{expr, partitionedIndices, iter},
  partitionedIndices = Partition[sigmaIndices, NumSigmaIndices[\[ScriptCapitalP][k]]];

  expr = 1;
  For[iter = 1, iter <= pow, iter++,
    expr = expr AvramidiToXTensor[\[ScriptCapitalP][k], IndexList[], partitionedIndices[[iter]]];
  ];

  expr
]


(*AddFreeIndices*)
AvramidiToXTensor[x:(AddFreeIndex[_,2]), freeIndices_IndexList, sigmaIndices_IndexList] :=
  Module[{nsi, vbundle, a, b},
  (* Get the vbundle corresponding to the index a *)
  vbundle = IndexVBundle;

  nsi = NumSigmaIndices[x];

  a = freeIndices[[-2]];
  b = freeIndices[[-1]];

  Apply[Plus, AvramidiToXTensor[x[[1]], freeIndices[[1;;-3]], #]& /@ 
    Flatten[CyclicPermutations[Join[IndexList[-a],#]]&/@
      CyclicPermutations[Join[IndexList[b],sigmaIndices[[1;;nsi]]]],1] ]/((nsi+2)(nsi+1))
]

AvramidiToXTensor[AbstractTrace[x:(AddFreeIndex[_,2])], freeIndices_IndexList, sigmaIndices_IndexList] :=
  Module[{nsi, newSigmaIndices, range, vbundle, contractedIndex},
  (* Get the vbundle corresponding to the index a *)
  vbundle = IndexVBundle;

  nsi = NumSigmaIndices[x];

  contractedIndex = DummyIn[vbundle];

  (* Insert the contracted index once cyclically *)
  newSigmaIndices = InsertCyclic[sigmaIndices, contractedIndex];

  (* Get the list of places to insert the second contracted index - from the first contracted index to the end *)
  range = -Reverse[Range[Range[nsi+1]]];

  (* Insert the second contracted index *)
  newSigmaIndices = Flatten[MapThread[ListInsert[#1, -contractedIndex, #2]&, {newSigmaIndices, range}],1];

  (* Produce the xTensor expression. The factor of 2 / ((nsi+2)(nsi+1)) accounts for the number of terms coming from symmetrization *)
  2 / ((nsi+2)(nsi+1)) Apply[Plus, AvramidiToXTensor[x[[1]], freeIndices, #]& /@ newSigmaIndices]
]

(* FIXME: This works for the term in V_ 1. Check it also works for all other cases. *)
AvramidiToXTensor[x:(AddFreeIndex[_,1]), freeIndices_IndexList, sigmaIndices_IndexList] :=
  Module[{nsi, vbundle, a},
  (* Get the vbundle corresponding to the index a *)
  vbundle = IndexVBundle;

  nsi = NumSigmaIndices[x];

  a = freeIndices[[-1]];

  Apply[Plus,AvramidiToXTensor[x[[1]], freeIndices[[1;;-2]],
    #]& /@ CyclicPermutations[Join[sigmaIndices[[1;;nsi]], IndexList[-a]]]] / (nsi+1)
]

AvramidiToXTensor[AbstractDot[AddFreeIndex[x_,1],y_], freeIndices_IndexList, sigmaIndices_IndexList] := Module[
{n1, n2, vbundle, contractedIndex, xnfi, ynfi, termFreeIndices},
  (*Get the vbundle corresponding to the index a*)
  vbundle = IndexVBundle;

  contractedIndex = DummyIn[vbundle];
 
  (* Figure out how many free indices each term uses *)
  xnfi = NumFreeIndices[AddFreeIndex[x,1]]-1;
  ynfi = NumFreeIndices[y]-1;

  (* And divide the indices up between each term *)
  termFreeIndices = PartitionIndices[freeIndices, {xnfi,ynfi}];

  n1 = NumSigmaIndices[x]-1;
  n2 = NumSigmaIndices[y];

  Apply[Plus, 
AvramidiToXTensor[x,termFreeIndices[[1]],#]& /@ CyclicPermutations[Join[sigmaIndices[[1;;n1]], IndexList[contractedIndex]]]]*
  AvramidiToXTensor[y,IndexList[contractedIndex,Sequence@@termFreeIndices[[2]]],sigmaIndices[[n1+1;;n1+n2]]] / (n1+1)
]

AvramidiToXTensor[AbstractDot[AddFreeIndex[x_,2],y_], freeIndices_IndexList, sigmaIndices_IndexList] := Module[
{n1, n2, vbundle, contractedIndex},
  (*Get the vbundle corresponding to the index a*)
  vbundle = IndexVBundle;

  contractedIndex = DummyIn[vbundle];

  n1 = NumSigmaIndices[x]-2;
  n2 = NumSigmaIndices[y];

  Apply[Plus, AvramidiToXTensor[x,IndexList[],#]& /@ Flatten[CyclicPermutations[Join[IndexList[-freeIndices[[1]]],#]]&/@CyclicPermutations[Join[IndexList[contractedIndex],sigmaIndices[[1;;n1]]]],1] ]*
  AvramidiToXTensor[y,IndexList[contractedIndex,freeIndices[[2]]],sigmaIndices[[n1+1;;n1+n2]]] / ((n1+2)(n1+1))
]

AvramidiToXTensor[AbstractDot[y_,AddFreeIndex[x_,1]], freeIndices_IndexList, sigmaIndices_IndexList] := Module[
{n1, n2, vbundle, contractedIndex, xnfi, ynfi, termFreeIndices},
  (*Get the vbundle corresponding to the index a*)
  vbundle = IndexVBundle;

  contractedIndex = DummyIn[vbundle];
 
  (* Figure out how many free indices each term uses *)
  xnfi = NumFreeIndices[AddFreeIndex[x,1]]-1;
  ynfi = NumFreeIndices[y]-1;

  (* And divide the indices up between each term *)
  termFreeIndices = PartitionIndices[freeIndices, {xnfi,ynfi}];

  n1 = NumSigmaIndices[x]-1;
  n2 = NumSigmaIndices[y];

  Apply[Plus, AvramidiToXTensor[x,termFreeIndices[[1]],#]& /@ CyclicPermutations[Join[sigmaIndices[[1;;n1]], IndexList[-contractedIndex]]]]*
  AvramidiToXTensor[y,IndexList[contractedIndex,termFreeIndices[[2]]],sigmaIndices[[n1+1;;n1+n2]]]/ (n1+1)
]

AvramidiToXTensor[AbstractDot[y_,AddFreeIndex[x_,2]], freeIndices_IndexList, sigmaIndices_IndexList] := Module[
{n1, n2, vbundle, contractedIndex},
  (*Get the vbundle corresponding to the index a*)
  vbundle = IndexVBundle;

  contractedIndex = DummyIn[vbundle];

  n1 = NumSigmaIndices[x]-2;
  n2 = NumSigmaIndices[y];

  Apply[Plus, AvramidiToXTensor[x,IndexList[],#]& /@ Flatten[CyclicPermutations[Join[IndexList[-freeIndices[[2]]],#]]&/@CyclicPermutations[Join[IndexList[contractedIndex],sigmaIndices[[1;;n1]]]],1] ]*
  AvramidiToXTensor[y,IndexList[freeIndices[[1]],contractedIndex],sigmaIndices[[n1+1;;n1+n2]]]/ ((n1+2)(n1+1))
]

(*AvramidiToXTensor[AbstractDot[\[ScriptCapitalR][1], \[ScriptCapitalK][2], AddFreeIndex[AbstractTrace[\[ScriptCapitalK][2]], 1]], IndexList[a_?AIndexQ, b_?AIndexQ], sigmaIndices_IndexList] :=
  Module[{vbundle, q, r, s, t},
  vbundle = IndexVBundle;
  q = NewIndexIn[vbundle]; r = NewIndexIn[vbundle];
  s = NewIndexIn[vbundle]; t = NewIndexIn[vbundle];

  
  (* First, create the Riemann tensor *)
  expr = Riemann[CD][-r, -q, -s, -sigmaIndices[[1]]]Riemann[CD][s, -sigmaIndices[[2]], r, -sigmaIndices[[3]]]Riemann[CD][t, q, -t, -sigmaIndices[[4]]]
  
]*)

(* FIXME: We need to treat this as a special case because otherwise we get {-beta, -beta} as free indices rather than {alpha, -beta}.
   Ideally we shouldn't have to worry about this special case *)
AvramidiToXTensor[AbstractTrace[x:(AbstractDot[AddFreeIndex[_, 1], _])], IndexList[], sigmaIndices_IndexList] := Module[
{vbundle, a},
  (*Get the vbundle corresponding to the index a*)
  vbundle = IndexVBundle;

  a = DummyIn[vbundle];

  AvramidiToXTensor[x, IndexList[a, -a], sigmaIndices]
]

AvramidiToXTensor[AbstractDot[Contraction[\[ScriptCapitalR][k_],pos_List],AddFreeIndex[x_, 1], y_], IndexList[], sigmaIndices_IndexList] := Module[
{n1, n2, vbundle, contractedIndex1, contractedIndex2},
  (*Get the vbundle corresponding to the index a*)
  vbundle = IndexVBundle;

  contractedIndex1 = DummyIn[vbundle];
  contractedIndex2 = DummyIn[vbundle];
  n1 = NumSigmaIndices[x]-1;
  n2 = NumSigmaIndices[y];

  AvramidiToXTensor[Contraction[\[ScriptCapitalR][k],pos],IndexList[-contractedIndex2],sigmaIndices[[1;;k]]]*
    Apply[Plus, 
      AvramidiToXTensor[x,IndexList[],#]& /@ CyclicPermutations[Join[sigmaIndices[[k+1;;k+n1]], IndexList[contractedIndex1]]]]*
      AvramidiToXTensor[y,IndexList[contractedIndex1,contractedIndex2],sigmaIndices[[k+n1+1;;k+n1+n2]]] / (n1+1)
]

AvramidiToXTensor[AbstractDot[Contraction[AbstractDot[\[ScriptCapitalR][k_],x_],pos_List],AddFreeIndex[y_, 1], z_], IndexList[], sigmaIndices_IndexList] := Module[
{n1, n2, n3, vbundle, contractedIndex1, contractedIndex2},
  (*Get the vbundle corresponding to the index a*)
  vbundle = IndexVBundle;

  contractedIndex1 = DummyIn[vbundle];
  contractedIndex2 = DummyIn[vbundle];

  n1 = NumSigmaIndices[x];
  n2 = NumSigmaIndices[y]-1;
  n3 = NumSigmaIndices[z];

  AvramidiToXTensor[Contraction[AbstractDot[\[ScriptCapitalR][k],x],pos],IndexList[-contractedIndex2],sigmaIndices[[1;;k+n1]]]*
    Apply[Plus,
      AvramidiToXTensor[y,IndexList[],#]& /@ CyclicPermutations[Join[sigmaIndices[[k+n1+1;;k+n1+n2]], IndexList[contractedIndex1]]]]*
      AvramidiToXTensor[z,IndexList[contractedIndex1,contractedIndex2],sigmaIndices[[k+n1+n2+1;;k+n1+n2+n3]]] / (n2+1)
]

End[] (* End Private Context *)

EndPackage[]



