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

(* Format Kappa_n and R_n nicely *)
Format[\[ScriptCapitalK][n_]] := Subscript[\[ScriptCapitalK], n];
Format[\[ScriptCapitalR][n_]] := Subscript[\[ScriptCapitalR], n];

(* We can calculate covariant series expansions for the following bitensors *)
Bitensors::usage = "Bitensors is a list of bitensors which can be expanded using CovariantSeries. The list is currently:\n"<>
	ToString[#[[1]] &/@ Bitensors]

Bitensors = { 
	{GammaBitensor, "\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"\[Mu]'\"], \[Nu]]\)"},
	{EtaBitensor,"\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"\[Mu]'\"], \[Nu]]\)"},
	{LambdaBitensor,"\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"\[Mu]\"], \[Nu]]\)"},
	{ZBitensor,"\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"\[Mu]'\"], \"\[Nu]'\"]\)"},
	{ABitensor,"\!\(\*SubscriptBox[SuperscriptBox[g,\"\[Mu]'\"], \"\[Nu] ; \[Tau]\"]\)"},
	{BBitensor,"\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"\[Mu]'\"], \"\[Nu] ; \[Tau]'\"]\)"},
	{ZetaBitensor,"\!\(\*SubscriptBox[SuperscriptBox[\[Sigma],\"\[Mu]'\"], \"\[Mu]'\"]\)"},
	{SqrtDeltaBitensor,"\!\(\*SuperscriptBox[\[CapitalDelta],\"1/2\"]\)"}
}

(* Usage messages for all the bitensors we support *)
(Evaluate[#[[1]]]::"usage" = ToString[#[[1]]] <> " is the bitensor " <> #[[2]] <> ".") & /@ Bitensors

AvramidiD::usage = "AvramidiD[x] is the derivative operator D on x."
AvramidiDPlus::usage = "AvramidiDPlus[x] is the part of the derivative operator D on x that increases the order in \!\(\*SuperscriptBox[\[Sigma], \"\[Mu]\"]\)."
AvramidiDSame::usage = "AvramidiDSame[x] is the part of the derivative operator D on x that preserves the order in \!\(\*SuperscriptBox[\[Sigma], \"\[Mu]\"]\)."

Begin["`Private`"]
(* Implementation of the package *)

(************************************** D *************************************)
(* Part that keeps the order the same *)
AvramidiDSame[AbstractDot[x_, y_]] :=  AbstractDot[AvramidiDSame[x], y] + AbstractDot[x, AvramidiDSame[y]]
AvramidiDSame[\[ScriptCapitalK][n_]] := n \[ScriptCapitalK][n]
AvramidiDSame[a_?NumericQ x_] := a AvramidiDSame[x]
e: AvramidiDSame[_Plus] :=  Distribute[Unevaluated[e]]
AvramidiDSame[a_?NumericQ] := a

(* Part that increases the order *)
AvramidiDPlus[AbstractDot[x_, y_]] :=  AbstractDot[AvramidiDPlus[x], y] + AbstractDot[x, AvramidiDPlus[y]]
AvramidiDPlus[\[ScriptCapitalK][n_]] := \[ScriptCapitalK][n + 1]
AvramidiDPlus[a_?NumericQ x_] := a AvramidiDPlus[x]
e: AvramidiDPlus[_Plus] := Distribute[Unevaluated[e]]
AvramidiDPlus[a_?NumericQ] := a

(* Total *)
AvramidiD[x_] := AvramidiDSame[x] + AvramidiDPlus[x]

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
	Expand[Sum[Binomial[n, k]*AbstractDot[ (n-k)*AvramidiDPlus[CovariantSeriesCoefficient[EtaBitensor,n-k-1]]
		 - AvramidiDSame[CovariantSeriesCoefficient[EtaBitensor,n-k]],
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

ABitensor /: CovariantSeriesCoefficient[ABitensor, 0] = -1;
ABitensor /: CovariantSeriesCoefficient[ABitensor, 1] = 0;

ABitensor /: CovariantSeriesCoefficient[ABitensor, n_]:= 
	ABitensor /: CovariantSeriesCoefficient[ABitensor, n] = 
		Expand[(1/(n + 1))*(n*\[ScriptCapitalR][n] 
			- Sum[Binomial[n, k]*AbstractDot[CovariantSeriesCoefficient[ABitensor,k],
				CovariantSeriesCoefficient[ZBitensor, n-k]], {k, 0, n - 2}])];

(************************************** B *************************************)
BBitensor /: CovariantSeries[BBitensor, n_]:= Sum[(-1)^i / i! CovariantSeriesCoefficient[BBitensor, i],{i,0,n}]

BBitensor /: CovariantSeriesCoefficient[BBitensor, 0] = -1;
BBitensor /: CovariantSeriesCoefficient[BBitensor, 1] = 0;

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

End[]

EndPackage[]

