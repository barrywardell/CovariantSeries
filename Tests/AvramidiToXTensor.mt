(* Mathematica Test File *)

Test[
	RiemannPart[\[ScriptCapitalK][3], b, -d, IndexList[c, a, e]] 
	,
	CD[-e][RiemannCD[b, -c, -d, -a]]
	,
	TestID->"RiemannPart"
]