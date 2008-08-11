(* Mathematica Package *)

(* Created by the Wolfram Workbench 24-Jul-2008 *)

BeginPackage["CovariantSeries`", "AbstractMatrix`"]
(* Exported symbols added here with SymbolName::usage *) 

(* The main functions we provide *)
CovariantSeries::usage = "CovariantSeries[X, n] calculates the covariant series of X up to "<>
	"order n, where X is one of the bitensors: \n" <> ToString[#[[1]] &/@ Bitensors];

CovariantSeriesCoefficient::usage = "CovariantSeriesCoefficient[X, n] calculates the order n "<>
	"coefficient of the covariant series of X,  where X is one of the bitensors: \n" <> ToString[#[[1]] &/@ Bitensors];

(* Kappa and R are the bitensors in terms of which all other bitensors are expanded *)
\[ScriptCapitalK]::usage = "\[ScriptCapitalK][n] is the \!\(\*SubscriptBox[K, \"(n)\"]\) of Avramidi."
\[ScriptCapitalR]::usage = "\[ScriptCapitalR][n] is the \!\(\*SubscriptBox[\[ScriptCapitalR], \"(n)\"]\) of Avramidi."
W::usage = "W[n] is the n-th term in the covariant expansion of the general bitensor, W."
m::usage = "m is the field mass."
g::usage = "g is the metric tensor."

(* Format Kappa_n and R_n nicely *)
Format[\[ScriptCapitalK][n_]] := Subscript[\[ScriptCapitalK], n];
Format[\[ScriptCapitalR][n_]] := Subscript[\[ScriptCapitalR], n];
Format[W[n_]] := Subscript[W, n];
Format[AddFreeIndex[x_,a_]] := Subscript[x, -a];

(* We can calculate covariant series expansions for the following bitensors *)
Bitensors = { 
	{GammaBitensor, "\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"\[Mu]'\"], \[Nu]]\)"},
	{EtaBitensor,"\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"\[Mu]'\"], \[Nu]]\)"},
	{LambdaBitensor,"\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"\[Mu]\"], \[Nu]]\)"},
	{ZBitensor,"\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"\[Mu]'\"], \"\[Nu]'\"]\)"},
	{ABitensor,"\!\(\*SubscriptBox[SuperscriptBox[g,\"\[Mu]'\"], \"\[Nu] ; \[Tau]\"]\)"},
	{BBitensor,"\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"\[Mu]'\"], \"\[Nu] ; \[Tau]'\"]\)"},
	{ZetaBitensor,"\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"\[Mu]'\"], \"\[Mu]'\"]\)"},
	{SqrtDeltaBitensor,"\!\(\*SuperscriptBox[\[CapitalDelta],\"1/2\"]\)"},
	{SqrtDeltaInvBitensor,"\!\(\*SuperscriptBox[\[CapitalDelta],\"-1/2\"]\)"},
	{CDSqrtDeltaBitensor,"\!\(\*SuperscriptBox[CD\[CapitalDelta],\"1/2\"]\)"},
	{BoxSqrtDeltaBitensor,"\!\(\*SuperscriptBox[Box\[CapitalDelta],\"1/2\"]\)"},
	{WBitensor,"\!\(\*SuperscriptBox[Box\[CapitalDelta],\"1/2\"]\)"},
	{SqrtDeltaInvDalWBitensor,"\!\(\*SuperscriptBox[Box\[CapitalDelta],\"1/2\"]\)"},
	{SqrtDeltaInvDalSqrtDeltaBitensor,"\!\(\*SuperscriptBox[Box\[CapitalDelta],\"1/2\"]\)"},
	{VTildeBitensor,"\!\(\*SuperscriptBox[Box\[CapitalDelta],\"1/2\"]\)"},
	{VBitensor,"\!\(\*SuperscriptBox[Box\[CapitalDelta],\"1/2\"]\)"},
	{TauBitensor,"\!\(\*SuperscriptBox[Box\[CapitalDelta],\"1/2\"]\)"},
	{TauPBitensor,"\!\(\*SuperscriptBox[Box\[CapitalDelta],\"1/2\"]\)"}
}

Bitensors::usage = "Bitensors is a list of bitensors which can be expanded using CovariantSeries. The list is currently:\n"<>
	ToString[#[[1]] &/@ Bitensors]

(Evaluate[#[[1]]]/: BitensorQ[Evaluate[#[[1]]]] = True) &/@ Bitensors
VBitensor/: BitensorQ[VBitensor[n_]] = True

(* Usage messages for all the bitensors we support *)
(Evaluate[#[[1]]]::"usage" = ToString[#[[1]]] <> " is the bitensor " <> #[[2]] <> ".") & /@ Bitensors

CovD::usage = "Covariant derivative of bitensor."
Dal::usage = "Covariant d'Alembertian of bitensor."
SigmaCD::usage = "SigmaCD[x] is the derivative operator D on x."
SigmaCDPlus::usage = "SigmaCDPlus[x] is the part of the derivative operator D on x that increases the order in \!\(\*SuperscriptBox[\[Sigma], \"\[Mu]\"]\)."
SigmaCDSame::usage = "SigmaCDSame[x] is the part of the derivative operator D on x that preserves the order in \!\(\*SuperscriptBox[\[Sigma], \"\[Mu]\"]\)."

Begin["`Private`"]
(* Implementation of the package *)

(************************************** D *************************************)
(* Part that keeps the order the same *)
SigmaCDSame[AbstractTrace[x_]] := AbstractTrace[SigmaCDSame[x]]
SigmaCDSame[AbstractDot[x_, y_]] :=  AbstractDot[SigmaCDSame[x], y] + AbstractDot[x, SigmaCDSame[y]]
SigmaCDSame[\[ScriptCapitalK][n_]] := n \[ScriptCapitalK][n]
SigmaCDSame[a_?NumericQ x_] := a SigmaCDSame[x]
e: SigmaCDSame[_Plus] :=  Distribute[Unevaluated[e]]
SigmaCDSame[a_?NumericQ] := a

(* Part that increases the order *)
SigmaCDPlus[AbstractTrace[x_]] := AbstractTrace[SigmaCDPlus[x]]
SigmaCDPlus[AbstractDot[x_, y_]] :=  AbstractDot[SigmaCDPlus[x], y] + AbstractDot[x, SigmaCDPlus[y]]
SigmaCDPlus[\[ScriptCapitalK][n_]] := \[ScriptCapitalK][n + 1]
SigmaCDPlus[a_?NumericQ x_] := a SigmaCDPlus[x]
e: SigmaCDPlus[_Plus] := Distribute[Unevaluated[e]]
SigmaCDPlus[a_?NumericQ] := a

(* Total *)
SigmaCD[x_] := SigmaCDSame[x] + SigmaCDPlus[x]

(******************************** Add Free Index ******************************)
AddFreeIndex[x_] := AddFreeIndex[x,1]
AddFreeIndex[AddFreeIndex[x_,a_],b_] := AddFreeIndex[x, a+b]
AddFreeIndex[a_?NumericQ,_] := a
AddFreeIndex[a_?NumericQ x_,b_] := a AddFreeIndex[x,b]
e: AddFreeIndex[_Plus,a_] := Distribute[Unevaluated[e]]

(*AddFreeIndex[x_] := x /. {(*\[ScriptCapitalL]->\[ScriptCapitalM],\[ScriptCapitalK]->\[ScriptCapitalL],*)
	 \[ScriptCapitalK][n_,m_]:>\[ScriptCapitalK][n,m+1],\[ScriptCapitalK][n_]:>\[ScriptCapitalK][n,1],
	 W[n_,m_]:>W[n,m+1], W[n_]:>W[n,1]} (* X->Y, W->X}*)*)

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

LambdaBitensor /: CovariantSeriesCoefficient[LambdaBitensor, 0] = -1;
LambdaBitensor /: CovariantSeriesCoefficient[LambdaBitensor, 1] = 0;

LambdaBitensor /: CovariantSeriesCoefficient[LambdaBitensor, n_]:= 
	LambdaBitensor /: CovariantSeriesCoefficient[LambdaBitensor, n] = 
	Expand[Sum[Binomial[n, k]*AbstractDot[ (n-k)*SigmaCDPlus[CovariantSeriesCoefficient[EtaBitensor,n-k-1]]
		 - SigmaCDSame[CovariantSeriesCoefficient[EtaBitensor,n-k]],
		  CovariantSeriesCoefficient[GammaBitensor,k]], {k, 0, n-2}]];

(************************************** Z *************************************)
ZBitensor /: CovariantSeries[ZBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[ZBitensor, i],{i,0,n}]

ZBitensor /: CovariantSeriesCoefficient[ZBitensor, 0] = -1;
ZBitensor /: CovariantSeriesCoefficient[ZBitensor, 1] = 0;

ZBitensor /: CovariantSeriesCoefficient[ZBitensor, n_]:= 
	ZBitensor /: CovariantSeriesCoefficient[ZBitensor, n] = 
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
				CovariantSeriesCoefficient[ZBitensor, n-k]], {k, 0, n - 2}])];

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

ZetaBitensor /: CovariantSeriesCoefficient[ZetaBitensor, 0] = 0;

ZetaBitensor /: CovariantSeriesCoefficient[ZetaBitensor, n_]:= 
	ZetaBitensor /: CovariantSeriesCoefficient[ZetaBitensor, n] = 
		SimplifyTrace[
			Expand[-(1/(2 n)) AbstractTrace[CovariantSeriesCoefficient[ZBitensor,n]]]
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

SqrtDeltaInvBitensor /: CovariantSeriesCoefficient[SqrtDeltaInvBitensor, n_?EvenQ] := 
	CovariantSeriesCoefficient[SqrtDeltaInvBitensor, n] = CovariantSeriesCoefficient[SqrtDeltaBitensor, n];
SqrtDeltaInvBitensor /: CovariantSeriesCoefficient[SqrtDeltaInvBitensor, n_?OddQ] := 
	CovariantSeriesCoefficient[SqrtDeltaInvBitensor, n] = -CovariantSeriesCoefficient[SqrtDeltaBitensor, n];

(**************************** Cov Deriv of Delta^1/2 **************************)
CDSqrtDeltaBitensor /: CovariantSeries[CDSqrtDeltaBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[CDSqrtDeltaBitensor, i],{i,0,n}]

CDSqrtDeltaBitensor /: CovariantSeriesCoefficient[CDSqrtDeltaBitensor, 0] = 0;

CDSqrtDeltaBitensor /: CovariantSeriesCoefficient[CDSqrtDeltaBitensor, n_]:= 
	CDSqrtDeltaBitensor /: CovariantSeriesCoefficient[CDSqrtDeltaBitensor, n] = 
		-Expand[Sum[Binomial[n, k-1]*AbstractDot[ AddFreeIndex[CovariantSeriesCoefficient[SqrtDeltaBitensor,k]],
			CovariantSeriesCoefficient[EtaBitensor,n+1-k] ], {k, 2, n+1}]];

(************************* D'alembertian of Delta^1/2 *************************)
BoxSqrtDeltaBitensor /: CovariantSeries[BoxSqrtDeltaBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[BoxSqrtDeltaBitensor, i],{i,0,n}]

BoxSqrtDeltaBitensor /: CovariantSeriesCoefficient[BoxSqrtDeltaBitensor, 0] :=  -AddFreeIndex[CovariantSeriesCoefficient[CDSqrtDeltaBitensor,1]];

BoxSqrtDeltaBitensor /: CovariantSeriesCoefficient[BoxSqrtDeltaBitensor, n_]:= 
	BoxSqrtDeltaBitensor /: CovariantSeriesCoefficient[BoxSqrtDeltaBitensor, n] = 
		Expand[-Sum[Binomial[n, k-1]*AbstractDot[ AddFreeIndex[CovariantSeriesCoefficient[CDSqrtDeltaBitensor,k]],
			CovariantSeriesCoefficient[EtaBitensor,n+1-k] ], {k, 2, n+1}]
			-Sum[Binomial[n, k]*AbstractDot[CovariantSeriesCoefficient[ABitensor,k],
				CovariantSeriesCoefficient[CDSqrtDeltaBitensor,n-k] ], {k, 2, n+1}]
		];

(**************** Covariant Derivative of a general bitensor ******************)
CovD /: CovariantSeries[CovD[x_?BitensorQ], n_] := Sum[(-1)^i / i! CovariantSeriesCoefficient[CovD[x], i],{i,0,n}]

CovD /: CovariantSeriesCoefficient[CovD[x_?BitensorQ], 0] :=
	CovD /: CovariantSeriesCoefficient[CovD[x], 0] = AddFreeIndex[CovariantSeriesCoefficient[x,1]];

CovD /: CovariantSeriesCoefficient[CovD[x_?BitensorQ], n_]:= 
	CovD /: CovariantSeriesCoefficient[CovD[x], n] = 
		Expand[-Sum[Binomial[n, k-1]*AbstractDot[ AddFreeIndex[CovariantSeriesCoefficient[x,k]],
			CovariantSeriesCoefficient[EtaBitensor,n+1-k] ], {k, 1, n+1}]];

(********************** D'alembertian of a general bitensor *******************)
Dal /: CovariantSeries[Dal[x_?BitensorQ], n_] := Sum[(-1)^i / i! CovariantSeriesCoefficient[Dal[x], i],{i,0,n}]

Dal /: CovariantSeriesCoefficient[Dal[x_?BitensorQ], 0] := 
	Dal /: CovariantSeriesCoefficient[Dal[x], 0] = AddFreeIndex[CovariantSeriesCoefficient[CovD[x],1]];

Dal /: CovariantSeriesCoefficient[Dal[x_?BitensorQ], n_]:= 
	Dal /: CovariantSeriesCoefficient[Dal[x], n] = 
		Expand[-Sum[Binomial[n, k-1]*AbstractDot[ AddFreeIndex[CovariantSeriesCoefficient[CovD[x],k]],
			CovariantSeriesCoefficient[EtaBitensor,n+1-k] ], {k, 1, n+1}]
			-Sum[Binomial[n, k]*AbstractDot[CovariantSeriesCoefficient[ABitensor,k],
				CovariantSeriesCoefficient[CovD[x],n-k] ], {k, 1, n}]
		];

(***************************** Genera Bitensor, W *****************************)
WBitensor /: CovariantSeries[WBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[WBitensor, i],{i,0,n}]

WBitensor /: CovariantSeriesCoefficient[WBitensor, n_] := W[n];

(***************************** Delta^-1/2 Dal[W] ******************************)
SqrtDeltaInvDalWBitensor /: CovariantSeries[SqrtDeltaInvDalWBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[SqrtDeltaInvDalWBitensor, i],{i,0,n}]

SqrtDeltaInvDalWBitensor /: CovariantSeriesCoefficient[SqrtDeltaInvDalWBitensor, n_]:= 
	SqrtDeltaInvDalWBitensor /: CovariantSeriesCoefficient[SqrtDeltaInvDalWBitensor, n] = 
		Expand[Sum[Binomial[n, k] CovariantSeriesCoefficient[SqrtDeltaInvBitensor, k]*
				(CovariantSeriesCoefficient[Dal[WBitensor], n-k]
				- m^2 CovariantSeriesCoefficient[WBitensor,n-k]), {k, 0, n}]
		];
		
(************************* Delta^-1/2 Dal[Delta^1/2] **************************)
SqrtDeltaInvDalSqrtDeltaBitensor /: CovariantSeries[SqrtDeltaInvDalSqrtDeltaBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[SqrtDeltaInvDalSqrtDeltaBitensor, i],{i,0,n}]

SqrtDeltaInvDalSqrtDeltaBitensor /: CovariantSeriesCoefficient[SqrtDeltaInvDalSqrtDeltaBitensor, n_]:= 
	SqrtDeltaInvDalSqrtDeltaBitensor /: CovariantSeriesCoefficient[SqrtDeltaInvDalSqrtDeltaBitensor, n] = 
		Expand[Sum[Binomial[n, k] CovariantSeriesCoefficient[SqrtDeltaInvBitensor, k]*
				(CovariantSeriesCoefficient[Dal[SqrtDeltaBitensor], n-k]
				- m^2 CovariantSeriesCoefficient[SqrtDeltaBitensor,n-k]), {k, 0, n}]
		];

(********************************* VTilde *************************************)
VTildeBitensor /: CovariantSeries[VTildeBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[VTildeBitensor, i],{i,0,n}]

VTildeBitensor /: CovariantSeriesCoefficient[VTildeBitensor, n_]:= 
	VTildeBitensor /: CovariantSeriesCoefficient[VTildeBitensor, n] = 
		Expand[-1/(n+1) (1/2 CovariantSeriesCoefficient[SqrtDeltaInvDalSqrtDeltaBitensor,n] +
			Binomial[n, 2] g CovariantSeriesCoefficient[SqrtDeltaInvDalWBitensor, n - 2])
		];

(************************************ V ****************************************)
VBitensor /: CovariantSeries[VBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[VBitensor, i],{i,0,n}]

VBitensor /: CovariantSeriesCoefficient[VBitensor, n_]:= 
	VBitensor /: CovariantSeriesCoefficient[VBitensor, n] = 
		Expand[Sum[Binomial[n, k] CovariantSeriesCoefficient[SqrtDeltaBitensor, k]*
				CovariantSeriesCoefficient[VTildeBitensor,n-k], {k, 0, n}]
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
				- 1/2 (CovariantSeriesCoefficient[Dal[SqrtDeltaBitensor], n] 
				- m^2 CovariantSeriesCoefficient[SqrtDeltaBitensor, n]))
		];

(************************************ Vl ***************************************)
VBitensor /: CovariantSeries[VBitensor[l_Integer?Positive], n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[VBitensor[l], i],{i,0,n}]

VBitensor /: CovariantSeriesCoefficient[VBitensor[l_Integer?Positive], n_]:= 
	VBitensor /: CovariantSeriesCoefficient[VBitensor[l], n] = 
		Expand[1/(n+l+1)( Sum[Binomial[n, k] CovariantSeriesCoefficient[VBitensor[l], k]*
				CovariantSeriesCoefficient[TauPBitensor,n-k], {k, 0, n-2}]
				- 1/(2 l) (CovariantSeriesCoefficient[Dal[VBitensor[l-1]], n] 
				- m^2 CovariantSeriesCoefficient[VBitensor[l-1], n]))
		];
		
End[]

EndPackage[]
