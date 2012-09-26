(* Mathematica Test File *)

(******************************* AbstractDot *********************************)
(* Numbers pull through *)
Test[
	AbstractDot[5, k ]
	,
	5 k
	,
	TestID->"Numeric 1"
]

Test[
	AbstractDot[k, 5 ]
	,
	5 k
	,
	TestID->"Numeric 2"
]

Test[
	AbstractDot[5 k[1], 10 k[2]]
	,
	50 AbstractDot[k[1], k[2]]
	,
	TestID->"Numeric 3"
]

Test[
	AbstractDot[5 k[1], 10 k[2], 5 k[3]]
	,
	250 AbstractDot[k[1], k[2], k[3]]
	,
	TestID->"Numeric 4"
]

(* Non-commutativity *)
Test[
	AbstractDot[k[1], k[2]] + AbstractDot[k[2], k[1]]
	,
	AbstractDot[k[1], k[2]] + AbstractDot[k[2], k[1]]
	,
	TestID->"Non-commutativity 1"
]

Test[
	AbstractDot[k[1], k[2]] + AbstractDot[k[1], k[2]]
	,
	2 AbstractDot[k[1], k[2]]
	,
	TestID->"Non-commutativity 2"
]

(* Distributivity *)
Test[
	AbstractDot[k[1] + k[3], k[2]]
	,
	AbstractDot[k[1], k[2]] + AbstractDot[k[3], k[2]]
	,
	TestID->"Distributivity 1"
]

Test[
	AbstractDot[k[1], k[2] + k[3]]
	,
	AbstractDot[k[1], k[2]] + AbstractDot[k[1], k[3]]
	,
	TestID->"Distributivity 2"
]

Test[
	Expand[AbstractDot[4k[1] + 3k[3], 2k[2]]]
	,
	8 AbstractDot[k[1], k[2]] + 6 AbstractDot[k[3], k[2]]
	,
	TestID->"Distributivity 3"
]

Test[
	Expand[AbstractDot[4 k[1], 3 k[2] + 2 k[3]]]
	,
	12 AbstractDot[k[1], k[2]] + 8 AbstractDot[k[1], k[3]]
	,
	TestID->"Distributivity 4"
]

(****************************** AbstractTrace ********************************)
Test[
	AbstractTrace[5 AbstractDot[k[1], k[2]]]
	,
	5 AbstractTrace[AbstractDot[k[1], k[2]]]
	,
	TestID->"AbstractTrace numeric factor"
]

Test[
	AbstractTrace[5]
	,
	5
	,
	TestID->"AbstractTrace numeric"
]

Test[
	AbstractTrace[AbstractDot[k[1], k[2]] + AbstractDot[k[2], k[1]]]
	,
	AbstractTrace[AbstractDot[k[1], k[2]]] + AbstractTrace[AbstractDot[k[2], k[1]]]
	,
	TestID->"AbstractTrace distributivity"
]

(****************************** SimplifyTrace ********************************)
Test[
	SimplifyTrace[AbstractTrace[AbstractDot[k[1], k[2]] + AbstractDot[k[2], k[1]]]]
	,
	2 AbstractTrace[AbstractDot[k[1], k[2]]]
	,
	TestID->"SimplifyTrace 1"
]

Test[
	SimplifyTrace[AbstractTrace[AbstractDot[k[3], k[2], k[6]] + AbstractDot[k[2], k[6], k[3]]]]
	,
	2 SimplifyTrace[AbstractTrace[AbstractDot[k[3], k[2], k[6]]]]
	,
	TestID->"SimplifyTrace 2"
]