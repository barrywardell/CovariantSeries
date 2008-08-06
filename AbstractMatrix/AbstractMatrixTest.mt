(* Mathematica Test File *)

(******************************* AbstractDot *********************************)
(* Numbers pull through *)
Test[
	AbstractDot[5, k ]
	,
	5 k
	,
	TestID->"AbstractMatrixTest-20080716-O7C7F2"
]

Test[
	AbstractDot[k, 5 ]
	,
	5 k
	,
	TestID->"AbstractMatrixTest-20080716-Z6M0O9"
]

Test[
	AbstractDot[5 k[1], 10 k[2]]
	,
	50 AbstractDot[k[1], k[2]]
	,
	TestID->"AbstractMatrixTest-20080716-Y5M9B5"
]

Test[
	AbstractDot[5 k[1], 10 k[2], 5 k[3]]
	,
	250 AbstractDot[k[1], k[2], k[3]]
	,
	TestID->"AbstractMatrixTest-20080806-J2L7S9"
]

(* Non-commutativity *)
Test[
	AbstractDot[k[1], k[2]] + AbstractDot[k[2], k[1]]
	,
	AbstractDot[k[1], k[2]] + AbstractDot[k[2], k[1]]
	,
	TestID->"AbstractMatrixTest-20080716-J9X4L5"
]

Test[
	AbstractDot[k[1], k[2]] + AbstractDot[k[1], k[2]]
	,
	2 AbstractDot[k[1], k[2]]
	,
	TestID->"AbstractMatrixTest-20080716-I5G2E1"
]

(* Distributivity *)
Test[
	AbstractDot[k[1] + k[3], k[2]]
	,
	AbstractDot[k[1], k[2]] + AbstractDot[k[3], k[2]]
	,
	TestID->"AbstractMatrixTest-20080716-R4K6X3"
]

Test[
	AbstractDot[k[1], k[2] + k[3]]
	,
	AbstractDot[k[1], k[2]] + AbstractDot[k[1], k[3]]
	,
	TestID->"AbstractMatrixTest-20080716-Y4X6R9"
]

Test[
	Expand[AbstractDot[4k[1] + 3k[3], 2k[2]]]
	,
	8 AbstractDot[k[1], k[2]] + 6 AbstractDot[k[3], k[2]]
	,
	TestID->"AbstractMatrixTest-20080717-P4Q1V8"
]

Test[
	Expand[AbstractDot[4 k[1], 3 k[2] + 2 k[3]]]
	,
	12 AbstractDot[k[1], k[2]] + 8 AbstractDot[k[1], k[3]]
	,
	TestID->"AbstractMatrixTest-20080717-Y8M4M0"
]

(****************************** AbstractTrace ********************************)
Test[
	AbstractTrace[5 AbstractDot[k[1], k[2]]]
	,
	5 AbstractTrace[AbstractDot[k[1], k[2]]]
	,
	TestID->"AbstractMatrixTest-20080716-N5I1P0"
]

Test[
	AbstractTrace[5]
	,
	5
	,
	TestID->"AbstractMatrixTest-20080717-Q6D2G8"
]

Test[
	AbstractTrace[AbstractDot[k[1], k[2]] + AbstractDot[k[2], k[1]]]
	,
	AbstractTrace[AbstractDot[k[1], k[2]]] + AbstractTrace[AbstractDot[k[2], k[1]]]
	,
	TestID->"AbstractMatrixTest-20080717-Y0U9T3"
]

(****************************** SimplifyTrace ********************************)
Test[
	SimplifyTrace[AbstractTrace[AbstractDot[k[1], k[2]] + AbstractDot[k[2], k[1]]]]
	,
	2 AbstractTrace[AbstractDot[k[1], k[2]]]
	,
	TestID->"AbstractMatrixTest-20080717-W0H7W2"
]

Test[
	SimplifyTrace[AbstractTrace[AbstractDot[k[3], k[2], k[6]] + AbstractDot[k[2], k[6], k[3]]]]
	,
	2 SimplifyTrace[AbstractTrace[AbstractDot[k[3], k[2], k[6]]]]
	,
	TestID->"AbstractMatrixTest-20080806-X7D0O5"
]