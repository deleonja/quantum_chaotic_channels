(* ::Package:: *)

BeginPackage["Chaometer`"];


Superoperator::usage = "Superoperator[t, EnvironentInitialState, eigenvaluesH, eigenvectosH, L] calculates the superoperator \!\(\*OverscriptBox[\(\[ScriptCapitalE]\), \(^\)]\).";


QubitChannelSuperoperator::usage = "QubitChannelSuperoperator[\[Eta],\[Kappa]] calculates the superoperator of single-qubit quantum channel transforming the Bloch sphere as \!\(\*OverscriptBox[\(r\), \(\[RightVector]\)]\) \[RightTeeArrow] \!\(\*OverscriptBox[\(\[Eta]\), \(\[RightVector]\)]\) . \!\(\*OverscriptBox[\(r\), \(\[RightVector]\)]\) + \!\(\*OverscriptBox[\(\[Kappa]\), \(\[RightVector]\)]\).";


Begin["`Private`"];


Get["QMB.wl"]


Superoperator[t_,\[Psi]0E_,eigVals_,eigVecs_,L_]:=
Module[{\[Psi]f,\[Rho]fSE(*SE final state*)},
Transpose[
\[Psi]f=StateEvolution[t,KroneckerVectorProduct[VectorFromKetInComputationalBasis[#],\[Psi]0E],eigVals,eigVecs]&/@{{0},{1}};
\[Rho]fSE=Dyad[\[Psi]f[[#1]],\[Psi]f[[#2]]]&@@@Tuples[{1,2},2];
Map[Flatten[MatrixPartialTrace[#,2,{2,2^(L-1)}]]&,\[Rho]fSE]
]
]


QubitChannelSuperoperator[\[Eta]_,\[Kappa]_]:=1/2*Reshuffle[{{1+\[Eta][[3]]+\[Kappa][[3]],0,\[Kappa][[1]]-I \[Kappa][[2]],\[Eta][[1]]+\[Eta][[2]]},{0,1-\[Eta][[3]]+\[Kappa][[3]],\[Eta][[1]]-\[Eta][[2]],\[Kappa][[1]]+I \[Kappa][[2]]},{\[Kappa][[1]]-I \[Kappa][[2]],\[Eta][[1]]-\[Eta][[2]],1-\[Eta][[3]]-\[Kappa][[3]],0},{\[Eta][[1]]+\[Eta][[2]],\[Kappa][[1]]-I \[Kappa][[2]],0,1+\[Eta][[3]]-\[Kappa][[3]]}}]


End[];


EndPackage[];
