(* ::Package:: *)

BeginPackage["Chaometer`"];


Superoperator::usage = "Superoperator[t, EnvironentInitialState, eigenvaluesH, eigenvectosH, L] calculates the superoperator \!\(\*OverscriptBox[\(\[ScriptCapitalE]\), \(^\)]\).";


Begin["`Private`"];


Superoperator[t_,\[Psi]0E_,eigVals_,eigVecs_,L_]:=
Module[{\[Psi]f,\[Rho]fSE(*SE final state*)},
Transpose[
\[Psi]f=StateEvolution[t,KroneckerVectorProduct[VectorFromKetInComputationalBasis[#],\[Psi]0E],eigVals,eigVecs]&/@{{0},{1}};
\[Rho]fSE=Dyad[\[Psi]f[[#1]],\[Psi]f[[#2]]]&@@@Tuples[{1,2},2];
Map[Flatten[MatrixPartialTrace[#,2,{2,2^(L-1)}]]&,\[Rho]fSE]
]
]


End[];


EndPackage[];
