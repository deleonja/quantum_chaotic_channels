(* ::Package:: *)

BeginPackage["Chaometer`"];


Superoperator::usage = "Superoperator[t, EnvironentInitialState, eigenvaluesH, eigenvectosH, L] calculates the superoperator \!\(\*OverscriptBox[\(\[ScriptCapitalE]\), \(^\)]\). ";


ChoiMatrix::usage = "Superoperator[U,A,L] computes the channel superoperator, where U is the unitary evolution operator of the chain of length L, and A=Id(2^(L-1))\[CircleTimes]\!\(\*TemplateBox[{SubscriptBox[\"\[Psi]\", \"E\"]},\n\"Ket\"]\)\!\(\*TemplateBox[{SubscriptBox[\"\[Psi]\", \"E\"]},\n\"Bra\"]\).";


QubitChannelSuperoperator::usage = "QubitChannelSuperoperator[\[Eta],\[Kappa]] calculates the superoperator of single-qubit quantum channel transforming the Bloch sphere as \!\(\*OverscriptBox[\(r\), \(\[RightVector]\)]\) \[RightTeeArrow] \!\(\*OverscriptBox[\(\[Eta]\), \(\[RightVector]\)]\) . \!\(\*OverscriptBox[\(r\), \(\[RightVector]\)]\) + \!\(\*OverscriptBox[\(\[Kappa]\), \(\[RightVector]\)]\).";


KrausOperatorsFromSuperoperator::usage = "KrausOperatorsFromSuperoperator[\!\(\*
StyleBox[\"superoperator\",\nFontSlant->\"Italic\"]\)] calculates the Kraus operators of \!\(\*
StyleBox[\"superoperator\",\nFontSlant->\"Italic\"]\).";


Begin["`Private`"];


Get["QMB.wl"]


Superoperator[t_,\[Psi]0E_,eigVals_,eigVecs_,L_]:=
	Module[{\[Psi]f,\[Rho]fSE(*SE final state*)},
		\[Psi]f=StateEvolution[t,KroneckerVectorProduct[VectorFromKetInComputationalBasis[#],\[Psi]0E],eigVals,eigVecs]&/@{{0},{1}};
		\[Rho]fSE=Dyad[\[Psi]f[[#1]],\[Psi]f[[#2]]]&@@@Tuples[{1,2},2];
		Transpose[ Flatten[MatrixPartialTrace[#,2,{2,2^(L-1)}]]&/@\[Rho]fSE ]
	]


Superoperator[t_,\[Psi]0E_,eigVals_,eigVecs_,L_,chaometer_]:=
	Module[{\[Psi]f,\[Rho]fSE(*SE final state*)},
		\[Psi]f=StateEvolution[t,KroneckerVectorProduct[\[Psi]0E[[1]],KroneckerVectorProduct[VectorFromKetInComputationalBasis[#],\[Psi]0E[[2]]]],eigVals,eigVecs]&/@{{0},{1}};
		\[Rho]fSE=Dyad[\[Psi]f[[#1]],\[Psi]f[[#2]]]&@@@Tuples[{1,2},2];
		Transpose[ Flatten[MatrixPartialTrace[#, Except[{chaometer}], ConstantArray[2,L] ]]&/@\[Rho]fSE ]
	]


ChoiMatrix[U_, A_, L_] := 
Module[{UR},
	UR = Reshuffle[U, 2, 2^(L-1)]; 
	Chop[UR . A . ConjugateTranspose[UR]](* A = Id(2^(L-1)\[CircleTimes])Dyad[Subscript[\[Psi], E]] *)
]


QubitChannelSuperoperator[\[Eta]_,\[Kappa]_]:=1/2*Reshuffle[{{1+\[Eta][[3]]+\[Kappa][[3]],0,\[Kappa][[1]]-I \[Kappa][[2]],\[Eta][[1]]+\[Eta][[2]]},{0,1-\[Eta][[3]]+\[Kappa][[3]],\[Eta][[1]]-\[Eta][[2]],\[Kappa][[1]]+I \[Kappa][[2]]},{\[Kappa][[1]]-I \[Kappa][[2]],\[Eta][[1]]-\[Eta][[2]],1-\[Eta][[3]]-\[Kappa][[3]],0},{\[Eta][[1]]+\[Eta][[2]],\[Kappa][[1]]-I \[Kappa][[2]],0,1+\[Eta][[3]]-\[Kappa][[3]]}}]


KrausOperatorsFromSuperoperator[superoperator_]:=ArrayReshape[Times[Sqrt[#1],#2]&@@Chop[Eigensystem[1/2Reshuffle[superoperator]]],{4,2,2}]


End[];


EndPackage[];
