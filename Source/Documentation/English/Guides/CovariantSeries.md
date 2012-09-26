{
 "Application" -> "CovariantSeries",
 "Package" -> "CovariantSeries",
 "Title" -> "CovariantSeries",
 "Summary" -> 
   "A package for computing covariant series expansions of many fundamental bi-tensors.",
 "Description" -> 
   {"The CovariantSeries package and the associated AbstractMatrix package provide ",
    "symbolic Mathematica code for calculating covariant series expansions of ",
    "fundamental bi-tensors including the world-function, the van Vleck determinant, ",
    "the bi-vector of parallel transport, the bi-scalar V(x,x') appearing in the ",
    "Hadamard form of the scalar Green function and the DeWitt coefficients. ",
    "A full description of the algorithm used by the package is available in the paper ",
    "\"Transport equation approach to calculations of Hadamard Green functions and ",
    "non-coincident DeWitt coefficients\", Adrian C. Ottewill and Barry Wardell, ",
    "Phys. Rev. D 84, 104039 (2011). If you find this package useful, please consider ",
    "citing that paper."},
 "Keywords" -> {"CovariantSeries", "Bi-tensor"},
 "Label" -> "CovariantSeries Application",
 "Synonyms" -> {"CovariantSeries"},
 "URL" -> "http://github.com/barrywardell/CovariantSeries" ,
 "Packages" -> {
   {"Title" -> "Computing Covariant Series Expansions",
    "DetailedFunctions" -> {
      {"CovariantSeries", "compute covariant series expansions"},
      {"CovariantSeriesCoefficient", "compute coefficients in covariant series expansions."},
      {"SetRicciFlat", "introduce the assumption that the Ricci tensor vanishes."}
    }
   },

   {"Title" -> "Supported Bitensors",
    "Functions" -> $Bitensors
   }
 },
 "Tutorials" -> {
   "CovariantSeries"
 } 
}