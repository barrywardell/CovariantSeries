(* Mathematica Test File *)

(* Checks against Decanini and Folacci (Phys.Rev.D73:044027,2006) *)

(* Gamma up to 11th order *)
Test[
	CovariantSeries[GammaBitensor,11]
	,
	-1 + \[ScriptCapitalK][2]/6 - \[ScriptCapitalK][3]/12 + 
 (-AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]/5 + (3*\[ScriptCapitalK][4])/5)/24 + 
 (AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]/3 + 
   (2*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2]])/3 - (2*\[ScriptCapitalK][5])/3)/120 + 
 ((-3*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4]])/7 - 
   (10*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3]])/7 - 
   (10*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][2]])/7 + 
   AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]/7 + (5*\[ScriptCapitalK][6])/7)/720 + 
 (AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5]]/2 + 
   (9*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4]])/4 + 
   (15*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][3]])/4 + 
   (5*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][2]])/2 - 
   AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]/4 - 
   AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2]]/2 - 
   (3*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/4 - (3*\[ScriptCapitalK][7])/4)/
  5040 + ((-5*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][6]])/9 - 
   (28*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][5]])/9 - 
   7*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][4]] - 
   (70*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][3]])/9 - 
   (35*AbstractDot[\[ScriptCapitalK][6], \[ScriptCapitalK][2]])/9 + 
   AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4]]/3 + 
   (10*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3]])/9 + 
   (10*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][2]])/9 + 
   (14*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][3]])/9 + 
   (28*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/9 + 
   (7*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/3 - 
   AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]/9 + 
   (7*\[ScriptCapitalK][8])/9)/40320 + ((3*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][7]])/5 + 
   4*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][6]] + 
   (56*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][5]])/5 + 
   (84*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][4]])/5 + 
   14*AbstractDot[\[ScriptCapitalK][6], \[ScriptCapitalK][3]] + 
   (28*AbstractDot[\[ScriptCapitalK][7], \[ScriptCapitalK][2]])/5 - 
   (2*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][5]])/5 - 
   (9*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][4]])/5 - 
   3*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][3]] - 
   2*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5], \[ScriptCapitalK][2]] - 
   (12*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][4]])/5 - 
   8*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][3]] - 
   8*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4], \[ScriptCapitalK][2]] - 
   (28*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][2], \[ScriptCapitalK][3]])/5 - 
   (56*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/5 - 
   (28*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/5 + 
   AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]/5 + 
   (2*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/5 + 
   (3*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/5 + 
   (4*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/5 - 
   (4*\[ScriptCapitalK][9])/5)/362880 + 
 ((-7*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][8]])/11 - 
   (54*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][7]])/11 - 
   (180*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][6]])/11 - 
   (336*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][5]])/11 - 
   (378*AbstractDot[\[ScriptCapitalK][6], \[ScriptCapitalK][4]])/11 - 
   (252*AbstractDot[\[ScriptCapitalK][7], \[ScriptCapitalK][3]])/11 - 
   (84*AbstractDot[\[ScriptCapitalK][8], \[ScriptCapitalK][2]])/11 + 
   (5*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][6]])/11 + 
   (28*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][5]])/11 + 
   (63*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][4]])/11 + 
   (70*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5], \[ScriptCapitalK][3]])/11 + 
   (35*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][6], \[ScriptCapitalK][2]])/11 + 
   (36*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][5]])/11 + 
   (162*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][4]])/11 + 
   (270*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4], \[ScriptCapitalK][3]])/11 + 
   (180*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][5], \[ScriptCapitalK][2]])/11 + 
   (108*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][2], \[ScriptCapitalK][4]])/11 + 
   (360*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][3], \[ScriptCapitalK][3]])/11 + 
   (360*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][4], \[ScriptCapitalK][2]])/11 + 
   (168*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][2], \[ScriptCapitalK][3]])/11 + 
   (336*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/11 + 
   (126*AbstractDot[\[ScriptCapitalK][6], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/11 - 
   (3*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4]])/11 - 
   (10*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3]])/11 - 
   (10*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][2]])/11 - 
   (14*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][3]])/11 - 
   (28*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/11 - 
   (21*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/11 - 
   (18*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]])/11 - 
   (36*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/11 - 
   (54*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/11 - 
   (36*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/11 + 
   AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]/11 + 
   (9*\[ScriptCapitalK][10])/11)/3628800 + 
 ((2*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][9]])/3 + 
   (35*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][8]])/6 + 
   (45*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][7]])/2 + 
   50*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][6]] + 
   70*AbstractDot[\[ScriptCapitalK][6], \[ScriptCapitalK][5]] + 
   63*AbstractDot[\[ScriptCapitalK][7], \[ScriptCapitalK][4]] + 
   35*AbstractDot[\[ScriptCapitalK][8], \[ScriptCapitalK][3]] + 
   10*AbstractDot[\[ScriptCapitalK][9], \[ScriptCapitalK][2]] - 
   AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][7]]/2 - 
   (10*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][6]])/3 - 
   (28*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][5]])/3 - 
   14*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5], \[ScriptCapitalK][4]] - 
   (35*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][6], \[ScriptCapitalK][3]])/3 - 
   (14*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][7], \[ScriptCapitalK][2]])/3 - 
   (25*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][6]])/6 - 
   (70*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][5]])/3 - 
   (105*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4], \[ScriptCapitalK][4]])/2 - 
   (175*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][5], \[ScriptCapitalK][3]])/3 - 
   (175*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][6], \[ScriptCapitalK][2]])/6 - 
   15*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][2], \[ScriptCapitalK][5]] - 
   (135*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][3], \[ScriptCapitalK][4]])/2 - 
   (225*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][4], \[ScriptCapitalK][3]])/2 - 
   75*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][5], \[ScriptCapitalK][2]] - 
   30*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][2], \[ScriptCapitalK][4]] - 
   100*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][3], \[ScriptCapitalK][3]] - 
   100*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][4], \[ScriptCapitalK][2]] - 
   35*AbstractDot[\[ScriptCapitalK][6], \[ScriptCapitalK][2], \[ScriptCapitalK][3]] - 
   70*AbstractDot[\[ScriptCapitalK][6], \[ScriptCapitalK][3], \[ScriptCapitalK][2]] - 
   21*AbstractDot[\[ScriptCapitalK][7], \[ScriptCapitalK][2], \[ScriptCapitalK][2]] + 
   AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][5]]/3 + 
   (3*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][4]])/2 + 
   (5*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][3]])/2 + 
   (5*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][5], \[ScriptCapitalK][2]])/3 + 
   2*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][4]] + 
   (20*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][3]])/3 + 
   (20*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][4], \[ScriptCapitalK][2]])/3 + 
   (14*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][2], \[ScriptCapitalK][3]])/3 + 
   (28*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/3 + 
   (14*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/3 + 
   (5*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4]])/2 + 
   (25*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3]])/3 + 
   (25*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][2]])/3 + 
   (35*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][3]])/3 + 
   (70*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/3 + 
   (35*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/2 + 
   (15*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]])/2 + 
   15*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2]] + 
   (45*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/2 + 
   10*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]] - 
   AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]/6 - 
   AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2]]/3 - 
   AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]/2 - 
   (2*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/3 - 
   (5*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/6 - 
   (5*\[ScriptCapitalK][11])/6)/39916800
	,
	TestID->"CovariantSeriesTest-20080806-K7X5V4"
]

(* Eta up to 9th order *)
Test[
	CovariantSeries[EtaBitensor, 9]
	,
	-1 - \[ScriptCapitalK][2]/6 + \[ScriptCapitalK][3]/12 + ((-7*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]])/15 - (3*\[ScriptCapitalK][4])/5)/24 + 
 ((4*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]])/3 + AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2]] + (2*\[ScriptCapitalK][5])/3)/120 + 
 ((-18*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4]])/7 - (25*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3]])/7 - 
   (11*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][2]])/7 - (31*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/21 - (5*\[ScriptCapitalK][6])/7)/
  720 + ((25*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5]])/6 + (33*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4]])/4 + 
   (27*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][3]])/4 + (13*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][2]])/6 + 
   (73*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]])/12 + (31*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/6 + 
   (17*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/4 + (3*\[ScriptCapitalK][7])/4)/5040 + 
 ((-55*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][6]])/9 - (140*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][5]])/9 - 
   (91*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][4]])/5 - (98*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][3]])/9 - 
   (25*AbstractDot[\[ScriptCapitalK][6], \[ScriptCapitalK][2]])/9 - (239*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4]])/15 - 
   (226*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3]])/9 - (106*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][2]])/9 - 
   (182*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][3]])/9 - (160*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/9 - 
   (43*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/5 - (127*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/15 - 
   (7*\[ScriptCapitalK][8])/9)/40320 + ((42*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][7]])/5 + 26*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][6]] + 
   (196*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][5]])/5 + (168*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][4]])/5 + 
   16*AbstractDot[\[ScriptCapitalK][6], \[ScriptCapitalK][3]] + (17*AbstractDot[\[ScriptCapitalK][7], \[ScriptCapitalK][2]])/5 + 
   (168*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][5]])/5 + (378*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][4]])/5 + 
   66*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][3]] + 22*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5], \[ScriptCapitalK][2]] + 
   60*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][4]] + 98*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][3]] + 
   47*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4], \[ScriptCapitalK][2]] + (232*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][2], \[ScriptCapitalK][3]])/5 + 
   (209*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/5 + (74*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/5 + 
   (226*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]])/5 + 
   (197*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/5 + 
   (184*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/5 + 
   31*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]] + (4*\[ScriptCapitalK][9])/5)/362880
	,
	TestID->"CovariantSeriesTest-20080806-E6O8H5"
]

(* Lambda up to 9th order *)
Test[
	CovariantSeries[LambdaBitensor, 9]
	,
	-1 - \[ScriptCapitalK][2]/3 + \[ScriptCapitalK][3]/12 + ((-8*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]])/15 - (2*\[ScriptCapitalK][4])/5)/24 + 
 (AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]] + AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2]] + \[ScriptCapitalK][5]/3)/120 + 
 ((-10*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4]])/7 - (17*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3]])/7 - 
   (10*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][2]])/7 - (32*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/21 - (2*\[ScriptCapitalK][6])/7)/
  720 + ((11*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5]])/6 + (17*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4]])/4 + 
   (17*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][3]])/4 + (11*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][2]])/6 + 
   (17*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]])/4 + (29*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/6 + 
   (17*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/4 + \[ScriptCapitalK][7]/4)/5040 + 
 ((-20*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][6]])/9 - (58*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][5]])/9 - 
   (44*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][4]])/5 - (58*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][3]])/9 - 
   (20*AbstractDot[\[ScriptCapitalK][6], \[ScriptCapitalK][2]])/9 - (42*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4]])/5 - 
   (146*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3]])/9 - (92*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][2]])/9 - 
   (124*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][3]])/9 - (146*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/9 - 
   (42*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/5 - (128*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/15 - 
   (2*\[ScriptCapitalK][8])/9)/40320 + ((13*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][7]])/5 + 9*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][6]] + 
   (77*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][5]])/5 + (77*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][4]])/5 + 
   9*AbstractDot[\[ScriptCapitalK][6], \[ScriptCapitalK][3]] + (13*AbstractDot[\[ScriptCapitalK][7], \[ScriptCapitalK][2]])/5 + 
   (71*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][5]])/5 + (187*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][4]])/5 + 
   40*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][3]] + 18*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5], \[ScriptCapitalK][2]] + 
   31*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][4]] + 62*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][3]] + 
   40*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4], \[ScriptCapitalK][2]] + 31*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][2], \[ScriptCapitalK][3]] + 
   (187*AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/5 + (71*AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/5 + 
   31*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]] + (181*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2]])/
    5 + (181*AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2]])/5 + 
   31*AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]] + \[ScriptCapitalK][9]/5)/362880
	,
	TestID->"CovariantSeriesTest-20080806-G6A4C9"
]

(* Zeta up to 11th order *)
Test[
	CovariantSeries[ZetaBitensor, 11]
	,
	AbstractTrace[\[ScriptCapitalK][2]]/12 - AbstractTrace[\[ScriptCapitalK][3]]/24 + 
 (AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]/15 + (3*AbstractTrace[\[ScriptCapitalK][4]])/10)/24 + 
 (-AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]/3 - AbstractTrace[\[ScriptCapitalK][5]]/3)/120 + 
 ((4*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4]]])/7 + (15*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3]]])/
    28 + (8*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]])/63 + (5*AbstractTrace[\[ScriptCapitalK][6]])/14)/
  720 + ((-5*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5]]])/6 - 
   (9*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4]]])/4 - 
   (4*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]])/3 - (3*AbstractTrace[\[ScriptCapitalK][7]])/8)/5040 + 
 ((10*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][6]]])/9 + (35*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][5]]])/
    9 + (14*AbstractTrace[AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][4]]])/5 + 
   (136*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4]]])/45 + 
   (50*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3]]])/9 + 
   (8*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]])/15 + (7*AbstractTrace[\[ScriptCapitalK][8]])/18)/
  40320 + ((-7*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][7]]])/5 - 
   6*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][6]]] - (56*AbstractTrace[AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][5]]])/5 - 
   (28*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][5]]])/5 - 
   (73*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][4]]])/5 - 
   (73*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][3]]])/5 - 
   9*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][3]]] - 
   (48*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]])/5 - (2*AbstractTrace[\[ScriptCapitalK][9]])/5)/
  362880 + ((56*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][8]]])/33 + 
   (189*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][7]]])/22 + 
   (216*AbstractTrace[AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][6]]])/11 + 
   (140*AbstractTrace[AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][5]]])/11 + 
   (304*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][6]]])/33 + 
   (1015*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][5]]])/33 + 
   (480*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][4]]])/11 + 
   (1015*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5], \[ScriptCapitalK][3]]])/33 + 
   81*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][4]]] + 
   (896*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4]]])/33 + 
   (149*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3]]])/3 + 
   (805*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]])/33 + 
   (128*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]])/33 + 
   (9*AbstractTrace[\[ScriptCapitalK][10]])/22)/3628800 + (-2*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][9]]] - 
   (35*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][8]]])/3 - (63*AbstractTrace[AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][7]]])/
    2 - 50*AbstractTrace[AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][6]]] - 
   14*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][7]]] - 
   (170*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][6]]])/3 - 
   103*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][5]]] - 
   103*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5], \[ScriptCapitalK][4]]] - 
   (95*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][6], \[ScriptCapitalK][3]]])/6 - 
   (245*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][6]]])/6 - 
   (575*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][5]]])/3 - 
   273*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4], \[ScriptCapitalK][4]]] - 
   (184*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][5]]])/3 - 
   (317*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][4]]])/2 - 
   (317*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][3]]])/2 - 
   (461*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][4]]])/3 - 
   (860*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][3]]])/3 - 
   (320*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]])/3 - 
   (5*AbstractTrace[\[ScriptCapitalK][11]])/12)/39916800
	,
	TestID->"CovariantSeriesTest-20080806-U7L1X0"
]

(* SqrtDelta up to 11th order *)
Test[
	CovariantSeries[SqrtDeltaBitensor, 11]
	,
	1 + AbstractTrace[\[ScriptCapitalK][2]]/12 - AbstractTrace[\[ScriptCapitalK][3]]/24 + 
 (AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]/15 + AbstractTrace[\[ScriptCapitalK][2]]^2/12 + 
   (3*AbstractTrace[\[ScriptCapitalK][4]])/10)/24 + (-AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]/3 - 
   (5*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][3]])/12 - AbstractTrace[\[ScriptCapitalK][5]]/3)/120 + 
 ((4*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4]]])/7 + (15*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3]]])/
    28 + (8*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]])/63 + 
   (AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][2]])/6 + 
   (5*AbstractTrace[\[ScriptCapitalK][2]]^3)/72 + (5*AbstractTrace[\[ScriptCapitalK][3]]^2)/8 + 
   (3*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][4]])/4 + (5*AbstractTrace[\[ScriptCapitalK][6]])/14)/720 + 
 ((-5*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5]]])/6 - (9*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4]]])/
    4 - (4*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]])/3 - 
   (7*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][2]])/6 - 
   (7*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][3]])/12 - 
   (35*AbstractTrace[\[ScriptCapitalK][2]]^2*AbstractTrace[\[ScriptCapitalK][3]])/48 - 
   (21*AbstractTrace[\[ScriptCapitalK][3]]*AbstractTrace[\[ScriptCapitalK][4]])/8 - 
   (7*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][5]])/6 - (3*AbstractTrace[\[ScriptCapitalK][7]])/8)/5040 + 
 ((7*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]^2)/45 + 
   (10*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][6]]])/9 + (35*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][5]]])/
    9 + (14*AbstractTrace[AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][4]]])/5 + 
   (136*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4]]])/45 + 
   (50*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3]]])/9 + 
   (8*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]])/15 + 
   (8*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4]]]*AbstractTrace[\[ScriptCapitalK][2]])/3 + 
   (5*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][2]])/2 + 
   (16*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][2]])/27 + 
   (7*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][2]]^2)/18 + 
   (35*AbstractTrace[\[ScriptCapitalK][2]]^4)/432 + (14*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*
     AbstractTrace[\[ScriptCapitalK][3]])/3 + (35*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][3]]^2)/12 + 
   (7*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][4]])/5 + 
   (7*AbstractTrace[\[ScriptCapitalK][2]]^2*AbstractTrace[\[ScriptCapitalK][4]])/4 + (63*AbstractTrace[\[ScriptCapitalK][4]]^2)/20 + 
   (14*AbstractTrace[\[ScriptCapitalK][3]]*AbstractTrace[\[ScriptCapitalK][5]])/3 + 
   (5*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][6]])/3 + (7*AbstractTrace[\[ScriptCapitalK][8]])/18)/40320 + 
 ((-14*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]])/5 - 
   (7*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][7]]])/5 - 6*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][6]]] - 
   (56*AbstractTrace[AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][5]]])/5 - 
   (28*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][5]]])/5 - 
   (73*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][4]]])/5 - 
   (73*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][3]]])/5 - 
   9*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][3]]] - 
   (48*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]])/5 - 
   5*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5]]]*AbstractTrace[\[ScriptCapitalK][2]] - 
   (27*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4]]]*AbstractTrace[\[ScriptCapitalK][2]])/2 - 
   8*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][2]] - 
   (7*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][2]]^2)/2 - 
   12*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4]]]*AbstractTrace[\[ScriptCapitalK][3]] - 
   (45*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][3]])/4 - 
   (8*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][3]])/3 - 
   (7*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][3]])/2 - 
   (35*AbstractTrace[\[ScriptCapitalK][2]]^3*AbstractTrace[\[ScriptCapitalK][3]])/24 - (35*AbstractTrace[\[ScriptCapitalK][3]]^3)/8 - 
   (63*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][4]])/5 - 
   (63*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][3]]*AbstractTrace[\[ScriptCapitalK][4]])/4 - 
   (14*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][5]])/5 - 
   (7*AbstractTrace[\[ScriptCapitalK][2]]^2*AbstractTrace[\[ScriptCapitalK][5]])/2 - 
   (63*AbstractTrace[\[ScriptCapitalK][4]]*AbstractTrace[\[ScriptCapitalK][5]])/5 - 
   (15*AbstractTrace[\[ScriptCapitalK][3]]*AbstractTrace[\[ScriptCapitalK][6]])/2 - 
   (9*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][7]])/4 - (2*AbstractTrace[\[ScriptCapitalK][9]])/5)/362880 + 
 (14*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]^2 + 8*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*
    AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4]]] + (56*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][8]]])/33 + 
   (15*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3]]])/2 + 
   (189*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][7]]])/22 + 
   (216*AbstractTrace[AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][6]]])/11 + 
   (140*AbstractTrace[AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][5]]])/11 + 
   (16*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]])/
    9 + (304*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][6]]])/33 + 
   (1015*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][5]]])/33 + 
   (480*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][4]]])/11 + 
   (1015*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5], \[ScriptCapitalK][3]]])/33 + 
   81*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][4]]] + 
   (896*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4]]])/33 + 
   (149*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3]]])/3 + 
   (805*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]])/33 + 
   (128*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]])/33 + 
   (7*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]^2*AbstractTrace[\[ScriptCapitalK][2]])/6 + 
   (25*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][6]]]*AbstractTrace[\[ScriptCapitalK][2]])/3 + 
   (175*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][5]]]*AbstractTrace[\[ScriptCapitalK][2]])/6 + 
   21*AbstractTrace[AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][4]]]*AbstractTrace[\[ScriptCapitalK][2]] + 
   (68*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4]]]*AbstractTrace[\[ScriptCapitalK][2]])/3 + 
   (125*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][2]])/3 + 
   4*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][2]] + 
   10*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4]]]*AbstractTrace[\[ScriptCapitalK][2]]^2 + 
   (75*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][2]]^2)/8 + 
   (20*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][2]]^2)/9 + 
   (35*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][2]]^3)/36 + 
   (35*AbstractTrace[\[ScriptCapitalK][2]]^5)/288 + 25*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5]]]*
    AbstractTrace[\[ScriptCapitalK][3]] + (135*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4]]]*AbstractTrace[\[ScriptCapitalK][3]])/2 + 
   40*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][3]] + 
   35*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][3]] + 
   (35*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][3]]^2)/4 + 
   (175*AbstractTrace[\[ScriptCapitalK][2]]^2*AbstractTrace[\[ScriptCapitalK][3]]^2)/16 + 
   36*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4]]]*AbstractTrace[\[ScriptCapitalK][4]] + 
   (135*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][4]])/4 + 
   8*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][4]] + 
   (21*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][4]])/2 + 
   (35*AbstractTrace[\[ScriptCapitalK][2]]^3*AbstractTrace[\[ScriptCapitalK][4]])/8 + 
   (315*AbstractTrace[\[ScriptCapitalK][3]]^2*AbstractTrace[\[ScriptCapitalK][4]])/8 + 
   (189*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][4]]^2)/8 + 
   28*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][5]] + 
   35*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][3]]*AbstractTrace[\[ScriptCapitalK][5]] + 14*AbstractTrace[\[ScriptCapitalK][5]]^2 + 
   5*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][6]] + 
   (25*AbstractTrace[\[ScriptCapitalK][2]]^2*AbstractTrace[\[ScriptCapitalK][6]])/4 + 
   (45*AbstractTrace[\[ScriptCapitalK][4]]*AbstractTrace[\[ScriptCapitalK][6]])/2 + 
   (45*AbstractTrace[\[ScriptCapitalK][3]]*AbstractTrace[\[ScriptCapitalK][7]])/4 + 
   (35*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][8]])/12 + (9*AbstractTrace[\[ScriptCapitalK][10]])/22)/3628800 + 
 (-88*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4]]] - 
   (55*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5]]])/3 - 
   2*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][9]]] - (165*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*
     AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3]]])/2 - 
   (99*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4]]])/2 - 
   (35*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][8]]])/3 - (63*AbstractTrace[AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][7]]])/
    2 - 50*AbstractTrace[AbstractDot[\[ScriptCapitalK][5], \[ScriptCapitalK][6]]] - 
   (176*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]])/
    9 - (88*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*
     AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]])/3 - 
   14*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][7]]] - 
   (170*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][6]]])/3 - 
   103*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][5]]] - 
   103*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5], \[ScriptCapitalK][4]]] - 
   (95*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][6], \[ScriptCapitalK][3]]])/6 - 
   (245*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][6]]])/6 - 
   (575*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][5]]])/3 - 
   273*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4], \[ScriptCapitalK][4]]] - 
   (184*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][5]]])/3 - 
   (317*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][4]]])/2 - 
   (317*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][3]]])/2 - 
   (461*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][2], \[ScriptCapitalK][4]]])/3 - 
   (860*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][3]]])/3 - 
   (320*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]])/3 - 
   (77*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*
     AbstractTrace[\[ScriptCapitalK][2]])/3 - (77*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][7]]]*AbstractTrace[\[ScriptCapitalK][2]])/
    6 - 55*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][6]]]*AbstractTrace[\[ScriptCapitalK][2]] - 
   (308*AbstractTrace[AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][5]]]*AbstractTrace[\[ScriptCapitalK][2]])/3 - 
   (154*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][5]]]*AbstractTrace[\[ScriptCapitalK][2]])/3 - 
   (803*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][4]]]*AbstractTrace[\[ScriptCapitalK][2]])/6 - 
   (803*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][2]])/6 - 
   (165*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][2]])/2 - 
   88*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][2]] - 
   (275*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5]]]*AbstractTrace[\[ScriptCapitalK][2]]^2)/12 - 
   (495*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4]]]*AbstractTrace[\[ScriptCapitalK][2]]^2)/8 - 
   (110*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][2]]^2)/3 - 
   (385*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][2]]^3)/36 - 
   (77*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]^2*AbstractTrace[\[ScriptCapitalK][3]])/12 - 
   (275*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][6]]]*AbstractTrace[\[ScriptCapitalK][3]])/6 - 
   (1925*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][5]]]*AbstractTrace[\[ScriptCapitalK][3]])/12 - 
   (231*AbstractTrace[AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][4]]]*AbstractTrace[\[ScriptCapitalK][3]])/2 - 
   (374*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4]]]*AbstractTrace[\[ScriptCapitalK][3]])/3 - 
   (1375*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][3]])/6 - 
   22*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][3]] - 
   110*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4]]]*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][3]] - 
   (825*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][3]])/8 - 
   (220*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][2]]*
     AbstractTrace[\[ScriptCapitalK][3]])/9 - (385*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*
     AbstractTrace[\[ScriptCapitalK][2]]^2*AbstractTrace[\[ScriptCapitalK][3]])/24 - 
   (1925*AbstractTrace[\[ScriptCapitalK][2]]^4*AbstractTrace[\[ScriptCapitalK][3]])/576 - 
   (385*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][3]]^2)/4 - 
   (1925*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][3]]^3)/48 - 
   (165*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5]]]*AbstractTrace[\[ScriptCapitalK][4]])/2 - 
   (891*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4]]]*AbstractTrace[\[ScriptCapitalK][4]])/4 - 
   132*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][4]] - 
   (231*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][4]])/2 - 
   (231*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][3]]*AbstractTrace[\[ScriptCapitalK][4]])/4 - 
   (1155*AbstractTrace[\[ScriptCapitalK][2]]^2*AbstractTrace[\[ScriptCapitalK][3]]*AbstractTrace[\[ScriptCapitalK][4]])/16 - 
   (2079*AbstractTrace[\[ScriptCapitalK][3]]*AbstractTrace[\[ScriptCapitalK][4]]^2)/16 - 
   88*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4]]]*AbstractTrace[\[ScriptCapitalK][5]] - 
   (165*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][5]])/2 - 
   (176*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][5]])/9 - 
   (77*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][5]])/3 - 
   (385*AbstractTrace[\[ScriptCapitalK][2]]^3*AbstractTrace[\[ScriptCapitalK][5]])/36 - 
   (385*AbstractTrace[\[ScriptCapitalK][3]]^2*AbstractTrace[\[ScriptCapitalK][5]])/4 - 
   (231*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][4]]*AbstractTrace[\[ScriptCapitalK][5]])/2 - 
   55*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]]*AbstractTrace[\[ScriptCapitalK][6]] - 
   (275*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][3]]*AbstractTrace[\[ScriptCapitalK][6]])/4 - 
   55*AbstractTrace[\[ScriptCapitalK][5]]*AbstractTrace[\[ScriptCapitalK][6]] - 
   (33*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]]*AbstractTrace[\[ScriptCapitalK][7]])/4 - 
   (165*AbstractTrace[\[ScriptCapitalK][2]]^2*AbstractTrace[\[ScriptCapitalK][7]])/16 - 
   (297*AbstractTrace[\[ScriptCapitalK][4]]*AbstractTrace[\[ScriptCapitalK][7]])/8 - 
   (385*AbstractTrace[\[ScriptCapitalK][3]]*AbstractTrace[\[ScriptCapitalK][8]])/24 - 
   (11*AbstractTrace[\[ScriptCapitalK][2]]*AbstractTrace[\[ScriptCapitalK][9]])/3 - (5*AbstractTrace[\[ScriptCapitalK][11]])/12)/39916800
	,
	TestID->"CovariantSeriesTest-20080806-H5R4P8"
]

(* Tau up to 9th order *)
Test[
	CovariantSeries[TauBitensor, 9]
	,
	AbstractTrace[\[ScriptCapitalK][2]]/6 - AbstractTrace[\[ScriptCapitalK][3]]/24 + 
 ((4*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2]]])/15 + AbstractTrace[\[ScriptCapitalK][4]]/5)/24 + 
 (-AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3]]] - AbstractTrace[\[ScriptCapitalK][5]]/6)/120 + 
 ((10*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4]]])/7 + (17*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3]]])/
    14 + (16*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]])/21 + AbstractTrace[\[ScriptCapitalK][6]]/7)/720 + 
 ((-11*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][5]]])/6 - (17*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][4]]])/
    4 - (20*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]])/3 - AbstractTrace[\[ScriptCapitalK][7]]/8)/5040 + 
 ((20*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][6]]])/9 + (58*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][5]]])/
    9 + (22*AbstractTrace[AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][4]]])/5 + 
   (608*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][4]]])/45 + 
   (208*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][3]]])/9 + 
   (64*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2]]])/15 + AbstractTrace[\[ScriptCapitalK][8]]/9)/
  40320 + ((-13*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][7]]])/5 - 
   9*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][6]]] - (77*AbstractTrace[AbstractDot[\[ScriptCapitalK][4], \[ScriptCapitalK][5]]])/5 - 
   (116*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][5]]])/5 - 
   (271*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][3], \[ScriptCapitalK][4]]])/5 - 
   (271*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][4], \[ScriptCapitalK][3]]])/5 - 
   31*AbstractTrace[AbstractDot[\[ScriptCapitalK][3], \[ScriptCapitalK][3], \[ScriptCapitalK][3]]] - 
   (336*AbstractTrace[AbstractDot[\[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][2], \[ScriptCapitalK][3]]])/5 - AbstractTrace[\[ScriptCapitalK][9]]/10)/
  362880
	,
	TestID->"CovariantSeriesTest-20080806-O0F3S8"
]

(* Checks against Christensen (PRD 17, 4 (1978)) *)

(* A up to order 4 *)
Test[
	CovariantSeries[ABitensor, 3]
	,
	-\[ScriptCapitalR][1]/2 + \[ScriptCapitalR][2]/3 + (-AbstractDot[\[ScriptCapitalR][1], \[ScriptCapitalK][2]]/4 - (3*\[ScriptCapitalR][3])/4)/6
	,
	TestID->"CovariantSeriesTest-20080806-M1D5V5"
]

(* B up to order 3 *)
Test[
	CovariantSeries[BBitensor, 3]
	,
	-\[ScriptCapitalR][1]/2 + \[ScriptCapitalR][2]/6 + (-AbstractDot[\[ScriptCapitalR][1], \[ScriptCapitalK][2]]/4 - \[ScriptCapitalR][3]/4)/6
	,
	TestID->"CovariantSeriesTest-20080806-U6B5M2"
]

(* CDSqrtDelta up to order 3 *)
Test[
	CovariantSeries[CDSqrtDeltaBitensor, 3]
	,
	-AddFreeIndex[AbstractTrace[\[ScriptCapitalK][2]], 1]/6 + 
 	AddFreeIndex[AbstractTrace[\[ScriptCapitalK][3]], 1]/8 + 
 	(-AbstractDot[AddFreeIndex[AbstractTrace[\[ScriptCapitalK][2]], 1], \[ScriptCapitalK][2]]/
    6 - AddFreeIndex[AbstractTrace[AbstractDot[\[ScriptCapitalK][2], 
       \[ScriptCapitalK][2]]], 1]/15 - AddFreeIndex[AbstractTrace[\[ScriptCapitalK][2]]^2, 
     1]/12 - (3*AddFreeIndex[AbstractTrace[\[ScriptCapitalK][4]], 1])/10)/6
	,
	TestID->"CovariantSeriesTest-20080806-A2R8M1"
]

(* CovD[SqrtDelta] up to order 3 *)
Test[
	CovariantSeries[CovD[SqrtDeltaBitensor],3]
	,
	Evaluate[CovariantSeries[CDSqrtDeltaBitensor, 3]]
	,
	TestID->"CovariantSeriesTest-20080806-A2R8M1"
]