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

Begin["`Private`"] (* Begin Private Context *) 

(* Convert K_n to n derivatives of Riemann in xTensor form *)
RiemannPart[\[ScriptCapitalK][n_], a_?AIndexQ, b_?AIndexQ] := Module[{vbundle, indices, CD, expr, i},
  (* Get the vbundle corresponding to the index a *)
  vbundle = VBundleOfIndex[a];

  (* Get some indices *)
  indices = GetIndicesOfVBundle[vbundle, n, {a, b, -a, -b}];
  
  (* Get the covariant derivative *)
  CD = CovDOfMetric[First[MetricsOfVBundle[vbundle]]];
  
  (* First, create the Riemann tensor *)
  expr = Riemann[CD][a, -indices[[1]], b, -indices[[2]]];
  
  (* Add covariant derivatives *)
  For[i = 3, i <= n, i++,
   expr = CD[-indices[[i]]]@expr;
  ];
  
  expr
]

AvramidiToXTensor[\[ScriptCapitalK][n_], a_?AIndexQ, b_?AIndexQ] := Module[{expr, indices, i},
  (* Riemann tensor part *)
  expr = RiemannPart[\[ScriptCapitalK][n],a,b];
  
  (* Extract the indices we need for the \[Sigma]'s *)
  indices = Cases[IndicesOf[Free][expr], Except[a|b]];

  (* Add sigmas *)
  For[i = 1, i <= Length[indices], i++,
    expr = expr \[Sigma][-indices[[i]]];
  ];
  
  ReplaceDummies[expr]
]

(* Convert AbstractDot[...] to Riemann in xTensor Form *)
AvramidiToXTensor[x_AbstractDot, a_?AIndexQ, b_?AIndexQ] := 
 Module[{vbundle, indices, i, expr, n},
  n = Length[x] - 1;
  
  (* Get the vbundle corresponding to the index a *)
  vbundle = VBundleOfIndex[a];

  (* Get some indices *)
  indices = Join[{a},GetIndicesOfVBundle[vbundle, n, {a, b, -a, -b}],{-b}];
  
  expr = 1;
  
  For[i = 1, i <= Length[x], i++,
  	expr = expr AvramidiToXTensor[x[[i]],indices[[i]],-indices[[i+1]]];
   ];
   
  expr
]

(* Convert AbstractTrace[...] to Riemann in xTensor form *)
AvramidiToXTensor[x_AbstractTrace, vbundle_?VBundleQ] := Module[{indices},
  (* Get some indices *)
  indices = GetIndicesOfVBundle[vbundle, 1];
  
  AvramidiToXTensor[#, indices[[1]], -indices[[1]]]& @@ x
]

(* In a sum, we treat each term independently *)
e : AvramidiToXTensor[_Plus, _, _] := Distribute[Unevaluated[e]]
e : AvramidiToXTensor[_Plus, _] := Distribute[Unevaluated[e]]

AvramidiToXTensor[x_Times, a_?AIndexQ, b_AIndexQ] := Map[AvramidiToXTensor[#,a,b]&, x]
AvramidiToXTensor[x_Times, v_?VBundleQ ] := Map[AvramidiToXTensor[#, v]&, x]

AvramidiToXTensor[a_?NumericQ, _?AIndexQ, _?AIndexQ] := a
AvramidiToXTensor[a_?NumericQ, v_?VBundleQ] := a

End[] (* End Private Context *)

EndPackage[]



