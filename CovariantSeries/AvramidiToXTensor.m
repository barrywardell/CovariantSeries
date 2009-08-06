(* AvramidiToXTensor Mathematica package 
   This is a Mathematica package for converting expressions from Avramidi notation
   to a the notation used by xTensor.
    
   Copyright 2009 Barry Wardell
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *)

BeginPackage["CovariantSeries`AvramidiToXTensor`", { "xAct`xTensor`", "CovariantSeries`", "AbstractMatrix`" }]
(* Exported symbols added here with SymbolName::usage *)  

AvramidiToXTensor::usage = 
  "AvramidiToXTensor is a package to convert expressions from the \
Avramidi notation (with \!\(\*SubscriptBox[\"\[ScriptCapitalK]\", \"n\
\"]\)) to Riemann tensors and their derivatives contracted with \
derivatives of \[Sigma] in xTensor form.";

RiemannPart::usage = 
  "RiemannPart[\!\(\*SubscriptBox[\"\[ScriptCapitalK]\", \"n\"]\)] \
gives the Riemann tensor part of \!\(\*SubscriptBox[\"s\", \"n\"]\) \
in xTensor form.";

\[Sigma]::usage = "\[Sigma] is Synge's world function";

$PrePrint=ScreenDollarIndices;
$CovDFormat="Postfix";

DefManifold[M, 4, {\[Alpha], \[Beta], a, b, c, d, e, f, h, i, j, k, l}];
DefMetric[-1, g[-a, -b], CD, {";", "\[Del]"}, CurvatureRelations -> False];
DefTensor[\[Sigma][a], M];

PrintAs[g] ^= "g";
PrintAs[epsilong] ^= "\[Epsilon]";
PrintAs[RiemannCD] ^= "R";
PrintAs[RicciCD] ^= "R";
PrintAs[RicciScalarCD] ^= "R";
PrintAs[WeylCD] ^= "W";
PrintAs[TFRicciCD] ^= "S";


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

PartitionIndices[indicesIn_, partition_]:= Module[{iter, expr, indices=indicesIn},
  expr = {};
    For[iter = 1, iter <= Length[partition], iter++,
      AppendTo[expr, Take[indices, partition[[iter]]]];
      indices = Drop[indices, partition[[iter]]];
    ];
  expr
]


(* Convert K_n to n derivatives of Riemann in xTensor form *)
RiemannPart[\[ScriptCapitalK][n_], a_?AIndexQ, b_?AIndexQ, indices_IndexList] := Module[{vbundle, CD, expr, i},
  (* Get the vbundle corresponding to the index a *)
  vbundle = VBundleOfIndex[a];

  (* Get the covariant derivative 
     FIXME: Should we worry about getting a different CD for each index? 
     Probably not, but maybe there is a case where it would be important. *)
  CD = CovDOfMetric[First[MetricsOfVBundle[vbundle]]];
  
  (* First, create the Riemann tensor *)
  expr = Riemann[CD][a, -indices[[1]], b, -indices[[2]]];
  
  (* Add covariant derivatives *)
  For[i = 3, i <= n, i++,
   expr = CD[-indices[[i]]]@expr;
  ];
  
  expr
]

(* In a sum, we treat each term independently *)
e : AvramidiToXTensor[_Plus, _, _] := Distribute[Unevaluated[e]]
e : AvramidiToXTensor[_Plus, _] := Distribute[Unevaluated[e]]


AvramidiToXTensor[a_?NumericQ, _] := a
AvramidiToXTensor[a_?NumericQ x_, y_] := a AvramidiToXTensor[x,y]

AvramidiToXTensor[CovariantSeries`m^(a_), _] := CovariantSeries`m^(a)
AvramidiToXTensor[CovariantSeries`m^(a_) x_, y_] := CovariantSeries`m^(a) AvramidiToXTensor[x,y]

AvramidiToXTensor[x_, vbundle_?VBundleQ] := Module[{expr, sigmaIndices, freeIndices, addFreeIndices, n, nf},
  (* Find how many indices we need *)
  n = NumSigmaIndices[x];
  nf = NumAddFreeIndices[x];
    
  (* Get some free indices *)
  freeIndices = {\[Alpha], -\[Beta]}; (*GetIndicesOfVBundle[vbundle, 2];*) (* FIXME: this could be different from 2 *)
  (*freeIndices[[2]] = -freeIndices[[2]];*)
  
  (* Get enough indices *)
  sigmaIndices = GetIndicesOfVBundle[vbundle, n, freeIndices];
  addFreeIndices = GetIndicesOfVBundle[vbundle, nf, Join[freeIndices, sigmaIndices]];

  expr = AvramidiToXTensor[x, IndexList@@freeIndices, IndexList @@ sigmaIndices, IndexList@@addFreeIndices];
  
  expr
]

(* Multiplication is a bit tricky since we want all terms to have unique indices *)
AvramidiToXTensor[x_Times, freeIndices_IndexList, sigmaIndices_IndexList, addFreeIndices_IndexList] := Module[{parts, indicesPerTerm, termIndices},
  (* Separate multiplication into a list of each term *)
  parts = List @@ x;
  
  (* Figure out how many indices each term uses *)
  indicesPerTerm = Map[NumSigmaIndices, parts];

  (* And divide the indices up between each term *)
  termIndices = PartitionIndices[sigmaIndices, indicesPerTerm];
  
  (* FIXME: we should have something other than freeIndices here *)
  Times@@MapThread[AvramidiToXTensor[#1, freeIndices, #2, addFreeIndices] &, {parts, List @@ termIndices}]
]

(* Powers - FIXME: assumes no free indices *)
AvramidiToXTensor[Power[x_AbstractTrace, pow_?Positive], freeIndices_IndexList, sigmaIndices_IndexList, addFreeIndices_IndexList] := Module[{expr, partitionedIndices, iter},
  partitionedIndices = Partition[sigmaIndices, NumSigmaIndices[x]];

  expr = 1;

  For[iter = 1, iter <= pow, iter++,
    expr = expr AvramidiToXTensor[x, freeIndices, partitionedIndices[[iter]], addFreeIndices];
  ];

  expr
]

(* AbstractDot *)
AvramidiToXTensor[x_AbstractDot, freeIndices_IndexList, sigmaIndices_IndexList, addFreeIndices_IndexList] := Module[{vbundle, numTerms, numContractedIndices, contractedIndices, iter, expr, sigmaIndicesPerTerm, indicesUsed},
  (* Get the vbundle corresponding to the index a *)
  vbundle = VBundleOfIndex[freeIndices[[1]]];
  numTerms = Length[x];
  numContractedIndices = numTerms - 1;

  contractedIndices = {};
  For[iter=1, iter<=numContractedIndices, iter++,
    contractedIndices = Append[contractedIndices, NewIndexIn[vbundle]];
  ];
  
  sigmaIndicesPerTerm = Map[NumSigmaIndices, List @@ x];

  expr = AvramidiToXTensor[x[[1]], IndexList[ freeIndices[[1]], -contractedIndices[[1]] ], sigmaIndices[[ 1 ;; sigmaIndicesPerTerm[[1]] ]], addFreeIndices ];
  indicesUsed = sigmaIndicesPerTerm[[1]];
  
  For[iter=2 , iter<numTerms, iter++,
  	expr = expr AvramidiToXTensor[x[[iter]], IndexList[ contractedIndices[[iter-1]], -contractedIndices[[iter]]], sigmaIndices[[ indicesUsed+1 ;; indicesUsed+1+sigmaIndicesPerTerm[[iter]] ]], addFreeIndices ];
  	indicesUsed += sigmaIndicesPerTerm[[iter]];
  ];
  
  expr = expr AvramidiToXTensor[x[[-1]], IndexList[contractedIndices[[-1]], freeIndices[[-1]]], sigmaIndices[[ indicesUsed + 1;; -1 ]], addFreeIndices ]
]

(* AddFreeIndices *)
AvramidiToXTensor[x_AddFreeIndex, freeIndices_IndexList, sigmaIndices_IndexList, addFreeIndices_IndexList] := Module[{res, n, nf, extraSigmaIndices, vbundle, metric, allSigmaIndices, originalExpr, expr, iter, iter2},
  (* Get the vbundle corresponding to the index a *)
  vbundle = VBundleOfIndex[freeIndices[[1]]];
  metric = First[MetricsOfVBundle[vbundle]];
  
  n = NumSigmaIndices[x];
  nf = NumAddFreeIndices[x];
  
  res = Insert[sigmaIndices, addFreeIndices[[1]], #] & /@ Range[1, n+1];
  extraSigmaIndices = IndexList@@GetIndicesOfVBundle[vbundle, nf, Join[freeIndices, sigmaIndices, addFreeIndices]];
  
  allSigmaIndices = Join[sigmaIndices, extraSigmaIndices];
  
  originalExpr = 1/(n+nf) AvramidiToXTensor[x[[1]], freeIndices, allSigmaIndices, addFreeIndices[[1;;nf]]];
  expr = 0;
  iter2=1;
  For[iter2=1, iter2<=nf, iter2++,
    For[iter=1, iter<=(n+nf), iter++,
  	  expr = expr + originalExpr metric[ allSigmaIndices[[iter]], -addFreeIndices[[iter2]] ];
    ];
  ];

  {allSigmaIndices,addFreeIndices,expr,res}
]

AvramidiToXTensor[\[ScriptCapitalK][n_], IndexList[a_?AIndexQ, b_?AIndexQ], sigmaIndices_IndexList, addFreeIndices_IndexList] := RiemannPart[\[ScriptCapitalK][n], a, b, sigmaIndices]

AvramidiToXTensor[AbstractTrace[x_], IndexList[a_?AIndexQ, b_?AIndexQ], sigmaIndices_IndexList, addFreeIndices_IndexList] := ReplaceDummies[ AvramidiToXTensor[x, IndexList[a, b], sigmaIndices, addFreeIndices] g[-a, -b] ]

End[] (* End Private Context *)

EndPackage[]



