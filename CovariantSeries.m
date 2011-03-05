(* ::Package:: *)

(* CovariantSeries Mathematica package 
   This is a Mathematica for calculating covariant series expansions
   of bitensors.
    
   Copyright 2009-2010 Barry Wardell
   
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
 
BeginPackage["CovariantSeries`", {"AbstractMatrix`"}]

(* Kappa and R are the bitensors in terms of which all other bitensors are expanded *)
\[ScriptCapitalK]::usage = "\[ScriptCapitalK][n] is the matrix \!\(\*SubscriptBox[\"K\", \"(n)\"]\) of Avramidi."
\[ScriptCapitalR]::usage = "\[ScriptCapitalR][n] is the tensor \!\(\*SubscriptBox[\[ScriptCapitalR], \"(n)\"]\) of Avramidi."
W::usage = "W[n] is the n-th term in the covariant expansion of the general bitensor, W."
\[ScriptCapitalB]::usage = "\[ScriptCapitalB][n] is the n-th term in the covariant expansion of the bitensor, \[ScriptCapitalB]."
\[ScriptCapitalP]::usage = "\[ScriptCapitalP][n] is the n-th term in the covariant expansion of the scalar potential, \[ScriptCapitalP]."
m::usage = "m is the field mass."
g::usage = "g is the metric tensor."

(* We can calculate covariant series expansions for the following bitensors *)
Bitensors = { 
	{GammaBitensor, "\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"a'\"], b]\)"},
	{EtaBitensor,"\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"a\"], \"b'\"]\)"},
	{LambdaBitensor,"\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"a\"], b]\)"},
	{XiBitensor,"\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"a'\"], \"b'\"]\)"},
	{ABitensor,"\!\(\*SubscriptBox[SuperscriptBox[g,\"a'\"], \"b ; c\"]\)"},
	{BBitensor,"\!\(\*SubscriptBox[SuperscriptBox[g,\"a'\"], \"b ; c'\"]\)"},
	{ZetaBitensor,"\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"\[Alpha]'\"], \"\[Alpha]'\"]\)"},
	{SqrtDeltaBitensor,"\!\(\*SuperscriptBox[\[CapitalDelta],\"1/2\"]\)"},
	{SqrtDeltaInvBitensor,"\!\(\*SuperscriptBox[\[CapitalDelta],\"-1/2\"]\)"},
	{CDSqrtDeltaBitensor,"\!\(\*SuperscriptBox[CD\[CapitalDelta],\"1/2\"]\)"},
	{BoxSqrtDeltaBitensor,"\!\(\*SuperscriptBox[Box\[CapitalDelta],\"1/2\"]\)"},
	{WBitensor,"W(x,x')"},
	{\[ScriptCapitalB]Bitensor,"\[ScriptCapitalB](x,x')"},
	{SqrtDeltaInvDalWBitensor,"\!\(\*SuperscriptBox[\[CapitalDelta],\"-1/2\"]BoxW(x,x')\)"},
	{SqrtDeltaInvDalSqrtDeltaBitensor,"\!\(\*SuperscriptBox[\[CapitalDelta],\"-1/2\"] \*SuperscriptBox[Box\[CapitalDelta],\"1/2\"]\)"},
	{VTildeBitensor,"Vtilde(x,x')"},
	{VBitensor,"V(x,x')"},
	{TauBitensor,"\[Tau](x,x')"},
	{TauPBitensor,"\[Tau]'(x,x')"},
	{DeWittABitensor,"a(x,x')"}
}

Bitensors::usage = "Bitensors is a list of bitensors which can be expanded using CovariantSeries. The list is currently:\n"<>
	ToString[#[[1]] &/@ Bitensors]

(* Usage messages for all the bitensors we support *)
(Evaluate[#[[1]]]::"usage" = ToString[#[[1]]] <> " is the bitensor " <> #[[2]] <> ".") & /@ Bitensors

AbstractCovD::usage = "Covariant derivative of bitensor."
AbstractDal::usage = "Covariant d'Alembertian of bitensor."
SigmaCD::usage = "SigmaCD[x] is the derivative operator D on x."
SigmaCDPlus::usage = "SigmaCDPlus[x] is the part of the derivative operator D on x that increases the order in \!\(\*SuperscriptBox[\[Sigma], \"\[Mu]\"]\)."
SigmaCDSame::usage = "SigmaCDSame[x] is the part of the derivative operator D on x that preserves the order in \!\(\*SuperscriptBox[\[Sigma], \"\[Mu]\"]\)."

(* The main functions we provide *)
CovariantSeries::usage = "CovariantSeries[X, n] calculates the covariant series of X up to "<>
	"order n, where X is one of the bitensors: \n" <> ToString[#[[1]] &/@ Bitensors];

CovariantSeriesCoefficient::usage = "CovariantSeriesCoefficient[X, n] calculates the order n "<>
	"coefficient of the covariant series of X,  where X is one of the bitensors: \n" <> ToString[#[[1]] &/@ Bitensors];

xTensorNotation::usage = "Specifies whether CovariantSeriesCoefficient should output in xTensor notation. Default: FALSE."
Options[CovariantSeriesCoefficient] = {xTensorNotation -> False};

AddFreeIndex::usage = "AddFreeIndex[] replaces one of the \[Sigma]-contracted indices with a free index."

SetRicciFlat::usage = "SetRicciFlat[] tries to simplify and speed up calculations by ignoring terms which are 0 in Ricci-flat spacetimes";

Begin["`Private`"]
(* Implementation of the package *)

(Evaluate[#[[1]]]/: BitensorQ[Evaluate[#[[1]]]] = True) &/@ Bitensors
VBitensor/: BitensorQ[VBitensor[_]] = True
BitensorQ[\[ScriptCapitalK]] = True
BitensorQ[\[ScriptCapitalR]] = True
BitensorQ[\[ScriptCapitalP]] = True
BiscalarQ[\[ScriptCapitalP][_]] = True
BitensorQ[_] = False
NotBitensorQ[x_] = If[BitensorQ[x],False,True]

(* Format Kappa_n and R_n nicely *)
Format[\[ScriptCapitalK][n_]] := Subscript[\[ScriptCapitalK], n];
Format[\[ScriptCapitalR][n_]] := Subscript[\[ScriptCapitalR], n];
Format[W[n_]] := Subscript[W, n];
Format[\[ScriptCapitalP][n_]] := Subscript[\[ScriptCapitalP], n];
Format[AddFreeIndex[x_,a_]] := Subscript[x, -a];

(* For Ricci-flat spacetimes, we can set tr (K[_]) = 0 *)
SetRicciFlat[] := Module[{},
  AbstractTrace[\[ScriptCapitalK][_]] = 0;
  Contraction[\[ScriptCapitalR][_],{2,3}] = 0;
  ];

(************************************** D *************************************)
(* Part that keeps the order the same *)
SigmaCDSame[AbstractTrace[x_]] := AbstractTrace[SigmaCDSame[x]]
SigmaCDSame[AbstractDot[x_, y_]] :=  AbstractDot[SigmaCDSame[x], y] + AbstractDot[x, SigmaCDSame[y]]
SigmaCDSame[\[ScriptCapitalK][n_]] := n \[ScriptCapitalK][n]
SigmaCDSame[\[ScriptCapitalP][n_]] := n \[ScriptCapitalP][n]
SigmaCDSame[a_?NumericQ x_] := a SigmaCDSame[x]
e: SigmaCDSame[_Plus] :=  Distribute[Unevaluated[e]]
SigmaCDSame[a_?NumericQ] := 0

(* Part that increases the order *)
SigmaCDPlus[AbstractTrace[x_]] := AbstractTrace[SigmaCDPlus[x]]
SigmaCDPlus[AbstractDot[x_, y_]] :=  AbstractDot[SigmaCDPlus[x], y] + AbstractDot[x, SigmaCDPlus[y]]
SigmaCDPlus[\[ScriptCapitalK][n_]] := \[ScriptCapitalK][n + 1]
SigmaCDPlus[\[ScriptCapitalP][n_]] := \[ScriptCapitalP][n + 1]
SigmaCDPlus[a_?NumericQ x_] := a SigmaCDPlus[x]
e: SigmaCDPlus[_Plus] := Distribute[Unevaluated[e]]
SigmaCDPlus[a_?NumericQ] := 0

(* Total *)
SigmaCD[x_] := SigmaCDSame[x] + SigmaCDPlus[x]

(******************************** Add Free Index ******************************)
AddFreeIndex[x_] := AddFreeIndex[x,1]
AddFreeIndex[AddFreeIndex[x_,a_],b_] := AddFreeIndex[x, a+b]
AddFreeIndex[a_?NumericQ,_] := a
AddFreeIndex[a_?NumericQ x_,b_] := a AddFreeIndex[x,b]
AddFreeIndex[m^(a_),_] := m^(a)
AddFreeIndex[m^(a_) x_,b_] := m^(a) AddFreeIndex[x,b]
AddFreeIndex[\[ScriptCapitalP][0],_] := \[ScriptCapitalP][0]
AddFreeIndex[\[ScriptCapitalP][0] x_,b_] := \[ScriptCapitalP][0] AddFreeIndex[x,b]
e: AddFreeIndex[_Plus,a_] := Distribute[Unevaluated[e]]

(*AddFreeIndex[x_] := x /. {(*\[ScriptCapitalL]->\[ScriptCapitalM],\[ScriptCapitalK]->\[ScriptCapitalL],*)
	 \[ScriptCapitalK][n_,m_]:>\[ScriptCapitalK][n,m+1],\[ScriptCapitalK][n_]:>\[ScriptCapitalK][n,1],
	 W[n_,m_]:>W[n,m+1], W[n_]:>W[n,1]} (* X->Y, W->X}*)*)

CovariantSeriesCoefficient[bitensor_, n_, OptionsPattern[]]:=
  If[OptionValue[xTensorNotation] == True,
    If[FreeQ[$ContextPath, "CovariantSeries`AvramidiToXTensor`"],
      Print["You must load the AvramidiToXTensor package first. Try '<<CovariantSeries`AvramidiToXTensor`'."];,
      CovariantSeries`AvramidiToXTensor`AvramidiToXTensor[CovariantSeriesCoefficient[bitensor,n], CovariantSeries`AvramidiToXTensor`TangentM]
    ],
    CovariantSeriesCoefficient[bitensor,n]
  ]

(************************************ gamma ***********************************)
GammaBitensor /: CovariantSeries[GammaBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[GammaBitensor, i],{i,0,n}]

GammaBitensor /: CovariantSeriesCoefficient[GammaBitensor, 0] = -1;
GammaBitensor /: CovariantSeriesCoefficient[GammaBitensor, 1] = 0;

GammaBitensor /: CovariantSeriesCoefficient[GammaBitensor, n_]:= 
	GammaBitensor /: CovariantSeriesCoefficient[GammaBitensor, n] = 
	Expand[(-((n - 1)/(n + 1)))*Sum[Binomial[n - 2, k]*AbstractDot[\[ScriptCapitalK][n - k], CovariantSeriesCoefficient[GammaBitensor,k]], {k, 0,n  - 2}]];

(************************************* eta ************************************)
EtaBitensor /: CovariantSeries[EtaBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[EtaBitensor, i],{i,0,n}]

EtaBitensor /: CovariantSeriesCoefficient[EtaBitensor, 0] = -1;
EtaBitensor /: CovariantSeriesCoefficient[EtaBitensor, 1] = 0;

EtaBitensor /: CovariantSeriesCoefficient[EtaBitensor, n_]:= 
	EtaBitensor /: CovariantSeriesCoefficient[EtaBitensor, n] = 
	Expand[Sum[Binomial[n, k]*AbstractDot[CovariantSeriesCoefficient[GammaBitensor,k], CovariantSeriesCoefficient[EtaBitensor,n - k]], {k, 2, n}]];

(*********************************** lambda ***********************************)
LambdaBitensor /: CovariantSeries[LambdaBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[LambdaBitensor, i],{i,0,n}]

LambdaBitensor /: CovariantSeriesCoefficient[LambdaBitensor, 0] = 1;
LambdaBitensor /: CovariantSeriesCoefficient[LambdaBitensor, 1] = 0;

LambdaBitensor /: CovariantSeriesCoefficient[LambdaBitensor, n_]:= 
	LambdaBitensor /: CovariantSeriesCoefficient[LambdaBitensor, n] = 
	Expand[Sum[Binomial[n, k]*AbstractDot[ (n-k)*SigmaCDPlus[CovariantSeriesCoefficient[EtaBitensor,n-k-1]]
		 - SigmaCDSame[CovariantSeriesCoefficient[EtaBitensor,n-k]],
		  CovariantSeriesCoefficient[GammaBitensor,k]], {k, 0, n-2}]];

(************************************* Xi *************************************)
XiBitensor /: CovariantSeries[XiBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[XiBitensor, i],{i,0,n}]

XiBitensor /: CovariantSeriesCoefficient[XiBitensor, 0] = 1;
XiBitensor /: CovariantSeriesCoefficient[XiBitensor, 1] = 0;

XiBitensor /: CovariantSeriesCoefficient[XiBitensor, n_]:= 
	XiBitensor /: CovariantSeriesCoefficient[XiBitensor, n] = 
		Expand[n*CovariantSeriesCoefficient[EtaBitensor, n] -
			 Sum[Binomial[n, k]*k*AbstractDot[CovariantSeriesCoefficient[GammaBitensor, n-k],
			 	 CovariantSeriesCoefficient[EtaBitensor, k]], {k, 2, n - 2}]];

(************************************** A *************************************)
ABitensor /: CovariantSeries[ABitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[ABitensor, i],{i,0,n}]

ABitensor /: CovariantSeriesCoefficient[ABitensor, 0] = 0;

ABitensor /: CovariantSeriesCoefficient[ABitensor, n_]:= 
	ABitensor /: CovariantSeriesCoefficient[ABitensor, n] = 
		Expand[(1/(n + 1))*(n*\[ScriptCapitalR][n] 
			- Sum[Binomial[n, k]*AbstractDot[CovariantSeriesCoefficient[ABitensor,k],
				CovariantSeriesCoefficient[XiBitensor, n-k]], {k, 0, n - 2}])];

(************************************** B *************************************)
BBitensor /: CovariantSeries[BBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[BBitensor, i],{i,0,n}]

BBitensor /: CovariantSeriesCoefficient[BBitensor, 0] = 0;

BBitensor /: CovariantSeriesCoefficient[BBitensor, n_]:= 
	BBitensor /: CovariantSeriesCoefficient[BBitensor, n] = 
		Expand[(1/n)*(CovariantSeriesCoefficient[ABitensor,n]
			- Sum[Binomial[n, k]*AbstractDot[CovariantSeriesCoefficient[ABitensor, k],
				CovariantSeriesCoefficient[EtaBitensor, n-k]], {k, 0, n - 2}])];

(************************************ Zeta ************************************)
ZetaBitensor /: CovariantSeries[ZetaBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[ZetaBitensor, i],{i,0,n}]

ZetaBitensor /: CovariantSeriesCoefficient[ZetaBitensor, 0] = 4;

ZetaBitensor /: CovariantSeriesCoefficient[ZetaBitensor, n_]:= 
	ZetaBitensor /: CovariantSeriesCoefficient[ZetaBitensor, n] = 
		SimplifyTrace[
			Expand[-(1/(2 n)) AbstractTrace[CovariantSeriesCoefficient[XiBitensor,n]]]
		];

(********************************** Delta^1/2 *********************************)
SqrtDeltaBitensor /: CovariantSeries[SqrtDeltaBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[SqrtDeltaBitensor, i],{i,0,n}]

SqrtDeltaBitensor /: CovariantSeriesCoefficient[SqrtDeltaBitensor, 0] = 1;
SqrtDeltaBitensor /: CovariantSeriesCoefficient[SqrtDeltaBitensor, 1] = 0;

SqrtDeltaBitensor /: CovariantSeriesCoefficient[SqrtDeltaBitensor, n_]:= 
	SqrtDeltaBitensor /: CovariantSeriesCoefficient[SqrtDeltaBitensor, n] = 
		SimplifyTrace[
			Expand[(1/n)*Sum[Binomial[n, k] k CovariantSeriesCoefficient[ZetaBitensor, k]*
				CovariantSeriesCoefficient[SqrtDeltaBitensor, n-k], {k, 1, n}]]
		];

(********************************** Delta^-1/2 *********************************)
SqrtDeltaInvBitensor /: CovariantSeries[SqrtDeltaInvBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[SqrtDeltaInvBitensor, i],{i,0,n}]

SqrtDeltaInvBitensor /: CovariantSeriesCoefficient[SqrtDeltaInvBitensor, 0] = 1;
SqrtDeltaInvBitensor /: CovariantSeriesCoefficient[SqrtDeltaInvBitensor, 1] = 0;

SqrtDeltaInvBitensor /: CovariantSeriesCoefficient[SqrtDeltaInvBitensor, n_] :=
	SqrtDeltaInvBitensor /: CovariantSeriesCoefficient[SqrtDeltaInvBitensor, n] =
        SimplifyTrace[
			Expand[-(1/n)*Sum[Binomial[n, k] k CovariantSeriesCoefficient[ZetaBitensor, k]*
				CovariantSeriesCoefficient[SqrtDeltaInvBitensor, n-k], {k, 1, n}]]
        ];

(**************************** Cov Deriv of Delta^1/2 **************************)
CDSqrtDeltaBitensor /: CovariantSeries[CDSqrtDeltaBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[CDSqrtDeltaBitensor, i],{i,0,n}]

CDSqrtDeltaBitensor /: CovariantSeriesCoefficient[CDSqrtDeltaBitensor, 0] = 0;

CDSqrtDeltaBitensor /: CovariantSeriesCoefficient[CDSqrtDeltaBitensor, n_]:= 
	CDSqrtDeltaBitensor /: CovariantSeriesCoefficient[CDSqrtDeltaBitensor, n] = 
		-Expand[Sum[Binomial[n, k-1]*AbstractDot[ AddFreeIndex[CovariantSeriesCoefficient[SqrtDeltaBitensor,k]],
			CovariantSeriesCoefficient[EtaBitensor,n+1-k] ], {k, 1, n+1}]];

(************************* D'alembertian of Delta^1/2 *************************)
BoxSqrtDeltaBitensor /: CovariantSeries[BoxSqrtDeltaBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[BoxSqrtDeltaBitensor, i],{i,0,n}]

BoxSqrtDeltaBitensor /: CovariantSeriesCoefficient[BoxSqrtDeltaBitensor, 0] :=  AddFreeIndex[CovariantSeriesCoefficient[CDSqrtDeltaBitensor,1]];

BoxSqrtDeltaBitensor /: CovariantSeriesCoefficient[BoxSqrtDeltaBitensor, n_]:= 
	BoxSqrtDeltaBitensor /: CovariantSeriesCoefficient[BoxSqrtDeltaBitensor, n] = 
		Expand[-Sum[Binomial[n, k-1]*AbstractDot[ AddFreeIndex[CovariantSeriesCoefficient[CDSqrtDeltaBitensor,k]],
			CovariantSeriesCoefficient[EtaBitensor,n+1-k] ], {k, 1, n+1}]
			-Sum[Binomial[n, k]*AbstractDot[Contraction[CovariantSeriesCoefficient[ABitensor,k],{2,3}],
				CovariantSeriesCoefficient[CDSqrtDeltaBitensor,n-k] ], {k, 1, n+1}]
		];

(**************** Covariant Derivative of a general bitensor ******************)
AbstractCovD /: CovariantSeries[AbstractCovD[x_?BitensorQ], n_] := Sum[(-1)^i / i! CovariantSeriesCoefficient[AbstractCovD[x], i],{i,0,n}]

AbstractCovD /: CovariantSeriesCoefficient[AbstractCovD[x_?BitensorQ], 0] :=
	AbstractCovD /: CovariantSeriesCoefficient[AbstractCovD[x], 0] = AddFreeIndex[CovariantSeriesCoefficient[x,1]];

AbstractCovD /: CovariantSeriesCoefficient[AbstractCovD[x_?BitensorQ], n_]:= 
	AbstractCovD /: CovariantSeriesCoefficient[AbstractCovD[x], n] = 
		Expand[-Sum[Binomial[n, k]*AbstractDot[ AddFreeIndex[CovariantSeriesCoefficient[x,k+1]],
			CovariantSeriesCoefficient[EtaBitensor,n-k] ], {k, 0, n}]];

(********************** D'alembertian of a general bitensor *******************)
AbstractDal /: CovariantSeries[AbstractDal[x_?BitensorQ], n_] := Sum[(-1)^i / i! CovariantSeriesCoefficient[AbstractDal[x], i],{i,0,n}]

AbstractDal /: CovariantSeriesCoefficient[AbstractDal[x_?BitensorQ], 0] := 
	AbstractDal /: CovariantSeriesCoefficient[AbstractDal[x], 0] = AbstractTrace[AddFreeIndex[CovariantSeriesCoefficient[AbstractCovD[x],1]]];

AbstractDal /: CovariantSeriesCoefficient[AbstractDal[x_?BitensorQ], n_]:= 
	AbstractDal /: CovariantSeriesCoefficient[AbstractDal[x], n] = 
		Expand[-Sum[Binomial[n, k]*AbstractTrace[AbstractDot[ AddFreeIndex[CovariantSeriesCoefficient[AbstractCovD[x],k+1]],
			CovariantSeriesCoefficient[EtaBitensor,n-k] ]], {k, 0, n}]
			-Sum[Binomial[n, k]*AbstractDot[Contraction[CovariantSeriesCoefficient[ABitensor,k],{2,3}],
				CovariantSeriesCoefficient[AbstractCovD[x],n-k] ], {k, 1, n}]
		];

(***************************** Genera Bitensor, W *****************************)
WBitensor /: CovariantSeries[WBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[WBitensor, i],{i,0,n}]

WBitensor /: CovariantSeriesCoefficient[WBitensor, n_] := W[n];

(***************************** Genera Bitensor, \[ScriptCapitalB] *****************************)
\[ScriptCapitalB]Bitensor /: CovariantSeries[\[ScriptCapitalB]Bitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[\[ScriptCapitalB]Bitensor, i],{i,0,n}]

\[ScriptCapitalB]Bitensor /: CovariantSeriesCoefficient[\[ScriptCapitalB]Bitensor, n_] /; n<0 = 0;
\[ScriptCapitalB]Bitensor /: CovariantSeriesCoefficient[\[ScriptCapitalB]Bitensor, n_] := \[ScriptCapitalB][n];

(***************************** Delta^-1/2 AbstractDal[W] ******************************)
SqrtDeltaInvDalWBitensor /: CovariantSeries[SqrtDeltaInvDalWBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[SqrtDeltaInvDalWBitensor, i],{i,0,n}]

SqrtDeltaInvDalWBitensor /: CovariantSeriesCoefficient[SqrtDeltaInvDalWBitensor, n_]:= 
	SqrtDeltaInvDalWBitensor /: CovariantSeriesCoefficient[SqrtDeltaInvDalWBitensor, n] = 
		Expand[Sum[Binomial[n, k] CovariantSeriesCoefficient[SqrtDeltaInvBitensor, k]*
				(CovariantSeriesCoefficient[AbstractDal[WBitensor], n-k]
				- m^2 CovariantSeriesCoefficient[WBitensor,n-k]), {k, 0, n}]
		];
		
(************************* Delta^-1/2 Dal[Delta^1/2] **************************)
SqrtDeltaInvDalSqrtDeltaBitensor /: CovariantSeries[SqrtDeltaInvDalSqrtDeltaBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[SqrtDeltaInvDalSqrtDeltaBitensor, i],{i,0,n}]

SqrtDeltaInvDalSqrtDeltaBitensor /: CovariantSeriesCoefficient[SqrtDeltaInvDalSqrtDeltaBitensor, n_]:= 
	SqrtDeltaInvDalSqrtDeltaBitensor /: CovariantSeriesCoefficient[SqrtDeltaInvDalSqrtDeltaBitensor, n] = 
		Expand[Sum[Binomial[n, k] CovariantSeriesCoefficient[SqrtDeltaInvBitensor, k]*
				(CovariantSeriesCoefficient[AbstractDal[SqrtDeltaBitensor], n-k]
				- m^2 CovariantSeriesCoefficient[SqrtDeltaBitensor,n-k]), {k, 0, n}]
		];

(********************************* VTilde *************************************)
VTildeBitensor /: CovariantSeries[VTildeBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[VTildeBitensor, i],{i,0,n}]

VTildeBitensor /: CovariantSeriesCoefficient[VTildeBitensor, n_]:= 
	VTildeBitensor /: CovariantSeriesCoefficient[VTildeBitensor, n] = 
		Expand[-1/(n+1) (1/2 CovariantSeriesCoefficient[SqrtDeltaInvDalSqrtDeltaBitensor,n](* -
			Binomial[n, 2] g CovariantSeriesCoefficient[SqrtDeltaInvDalWBitensor, n - 2]*))
		];

(************************************ V ****************************************)
VBitensor /: CovariantSeries[VBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[VBitensor, i],{i,0,n}]

VBitensor /: CovariantSeriesCoefficient[VBitensor, n_]:= 
	VBitensor /: CovariantSeriesCoefficient[VBitensor, n] = 
		Expand[Sum[Binomial[n, k] CovariantSeriesCoefficient[SqrtDeltaBitensor, k]*
				CovariantSeriesCoefficient[VTildeBitensor,n-k], {k, 0, n}] 
				+ g CovariantSeriesCoefficient[\[ScriptCapitalB]Bitensor,n-2]
		];

(*********************************** tau ***************************************)
TauBitensor /: CovariantSeries[TauBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[TauBitensor, i],{i,0,n}]

TauBitensor /: CovariantSeriesCoefficient[TauBitensor, 0] = 0;
TauBitensor /: CovariantSeriesCoefficient[TauBitensor, 1] = 0;

TauBitensor /: CovariantSeriesCoefficient[TauBitensor, n_]:= 
	TauBitensor /: CovariantSeriesCoefficient[TauBitensor, n] = 
	Expand[SimplifyTrace[-n*SigmaCDPlus[CovariantSeriesCoefficient[ZetaBitensor,n-1]]
		 + SigmaCDSame[CovariantSeriesCoefficient[ZetaBitensor,n]]]
	];

(*********************************** tau ***************************************)
TauPBitensor /: CovariantSeries[TauPBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[TauPBitensor, i],{i,0,n}]

TauPBitensor /: CovariantSeriesCoefficient[TauPBitensor, 0] = 0;
TauPBitensor /: CovariantSeriesCoefficient[TauPBitensor, 1] = 0;

TauPBitensor /: CovariantSeriesCoefficient[TauPBitensor, n_]:= 
	TauPBitensor /: CovariantSeriesCoefficient[TauPBitensor, n] = 
	Expand[n CovariantSeriesCoefficient[ZetaBitensor,n]];
	
(************************************ V0 ***************************************)
VBitensor /: CovariantSeries[VBitensor[0], n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[VBitensor[0], i],{i,0,n}]

VBitensor /: CovariantSeriesCoefficient[VBitensor[0], n_]:= 
	VBitensor /: CovariantSeriesCoefficient[VBitensor[0], n] = 
		Expand[1/(n+1)( Sum[Binomial[n, k] CovariantSeriesCoefficient[VBitensor[0], k]*
				CovariantSeriesCoefficient[TauPBitensor,n-k], {k, 0, n-2}]
				- 1/2 (CovariantSeriesCoefficient[AbstractDal[SqrtDeltaBitensor], n] 
				- m^2 CovariantSeriesCoefficient[SqrtDeltaBitensor, n]-
                Sum[Binomial[n, k] CovariantSeriesCoefficient[SqrtDeltaBitensor, k]*
				\[ScriptCapitalP][n-k], {k, 0, n}]))
		];

(************************************ Vl ***************************************)
VBitensor /: CovariantSeries[VBitensor[l_Integer?Positive], n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[VBitensor[l], i],{i,0,n}]

VBitensor /: CovariantSeriesCoefficient[VBitensor[l_Integer?Positive], n_]:= 
	VBitensor /: CovariantSeriesCoefficient[VBitensor[l], n] = 
		Expand[1/(n+l+1)( Sum[Binomial[n, k] CovariantSeriesCoefficient[VBitensor[l], k]*
				CovariantSeriesCoefficient[TauPBitensor,n-k], {k, 0, n-2}]
				- 1/(2 l) (CovariantSeriesCoefficient[AbstractDal[VBitensor[l-1]], n] 
				- m^2 CovariantSeriesCoefficient[VBitensor[l-1], n]-
                Sum[Binomial[n, k] CovariantSeriesCoefficient[VBitensor[l-1], k]*
				\[ScriptCapitalP][n-k], {k, 0, n}]))
		];

(************************************ a_k **************************************)
DeWittABitensor /: CovariantSeries[DeWittABitensor[l_Integer?NonNegative], n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[DeWittABitensor[l], i],{i,0,n}]

DeWittABitensor /: CovariantSeriesCoefficient[DeWittABitensor[0], 0] = 1;
DeWittABitensor /: CovariantSeriesCoefficient[DeWittABitensor[0], _] = 0;

DeWittABitensor /: CovariantSeriesCoefficient[DeWittABitensor[l_Integer?Positive], n_]:=
	DeWittABitensor /: CovariantSeriesCoefficient[DeWittABitensor[l], n]=
		Expand[2^l (l-1)! (-1)^l Sum[ Binomial[n,k] CovariantSeriesCoefficient[SqrtDeltaInvBitensor, k] ReplaceAll[CovariantSeriesCoefficient[VBitensor[l-1],n-k], m->0], {k,0,n}]];

End[]

EndPackage[]
