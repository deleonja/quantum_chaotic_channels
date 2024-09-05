(* ::Package:: *)

BeginPackage["QuantumJA`"]


\[Omega]::usage="\[Omega][d] returns the first dth root of unity.
\[Omega][d,k] returns \[Omega][d]^k.";
WeylMatrix::usage="WeylMatrix[m, n, d] returns the multiparticle Weyl matrix U({m1,...,mN},{n1,...,nN}) of dimension {d1,...,dN}.";
ComplexToPolar::usage="ComplexToPolar[z] returns the polar form of a complex number.";
Commutator::usage="Commutator[A,B] returns A.B-B.A";
CommutationQ::usage="CommutationQ[A,B] yields True if A and B commute, and False otherwise.";
MutuallyCommutingSetQ::usage="MutuallyCommutingSetQ[ListOfMatrices] yields True if all matrices in the list mutually commute, and False otherwise.";
(*Funci\[OAcute]n prestada de *)
Pauli::usage= "Pauli[0-3] gives Pauli Matrices according to wikipedia, and Pauli[{i1,i2,...,in}] gives Pauli[i1] \[CircleTimes]Pauli[i2] \[CircleTimes] ... \[CircleTimes] Pauli[in]";
EigenvectorQ::usage="EigenvectorQ[M,v] returns True if vector v is an eigenvector of matrix M.";

(*Cosas espec\[IAcute]ficas de Weyl channels*)
PlotGroup::usage="PlotGroup[H_] plots a group H of the form Zd+Zd.";
A::usage="A[d_List] returns matrix A of Weyl channels of {d1,d2,...,dN}.";
GeneratorSubgroups::usage="GeneratorSubgroups[d] returns the generator subgroups of Zd+Zd";
ChoiMatrix::usage="ChoiMatrix[Tau, d] returns the Choi matrix of a map Tau={Tau1,...,Taud1^2...dN^2}";
PToverBofChoiMatrix::usage="PToverBofChoiMatrix[Tau, d] returns the partial transpose of the Choi matrix of a Weyl map Tau of a system of particles with dimensions d={d1,...,dN}.";
TauOfWCEFromSubgroup::usage="TauOfWCEFromSubgroup[Indices, TauValues, d] returns the vector Tau of a Weyl map given the TauValues different from zero and its Indices.";
TausOfWCEGs::usage="TausOfWCEGs[d] returns the vectors Tau of those single-qudit Weyl erasing channels associated with a generator subgroup.";


Begin["`Private`"]


\[Omega][d_]:=Exp[2*Pi*I/d]
\[Omega][d_,k_]:=Exp[2k*Pi*I/d]


WeylMatrix[m_Integer,n_Integer,d_Integer]:=Sum[SparseArray[Mod[{k,k+n},d]+1->\[Omega][d,k*m],{d,d}],{k,0,d-1}]
WeylMatrix[m_List,n_List,d_List]:=KroneckerProduct@@(WeylMatrix[#1,#2,#3]&@@@Transpose[Join[{m,n},{d}]])


ComplexToPolar[z_]/;z\[Element]Complexes:=Abs[z] Exp[I Arg[z]]


Commutator[A_,B_]:=A . B-B . A


ZeroMatrix[d_]:=ConstantArray[0,{d,d}]


CommutationQ[A_,B_]:=Commutator[A,B]==ZeroMatrix[Length[A]]


MutuallyCommutingSetQ[ListOfMatrices_]:=Module[{SetLength=Length[ListOfMatrices]},
AllTrue[Table[CommutationQ@@ListOfMatrices[[{i,j}]],{i,SetLength-1},{j,i+1,SetLength}],TrueQ,2]
]


(* {{{ Pauli matrices*)  
Pauli[0]=Pauli[{0}]=SparseArray[{{1,0}, {0,1}}]; 
Pauli[1]=Pauli[{1}]=SparseArray[{{0,1}, {1,0}}]; 
Pauli[2]=Pauli[{2}]=SparseArray[{{0,-I},{I,0}}]; 
Pauli[3]=Pauli[{3}]=SparseArray[{{1,0}, {0,-1}}];
Pauli[Indices_List] := KroneckerProduct @@ (Pauli /@ Indices)


EigenvectorQ[M_,v_]:=MatrixRank[{FullSimplify/@M . v,v}]==1


(*Cosas espec\[IAcute]ficas de los Weyl channels*)
PlotGroup[H_List]:=Module[{d=Max[H]+1},ArrayPlot[SparseArray[#+1->ConstantArray[1,Length[#]],{d,d}],ImageSize->100]&[H]]
PlotGroup[H_List,options_]:=Module[{d=Max[H]+1},ArrayPlot[SparseArray[#+1->ConstantArray[1,Length[#]],{d,d}],options]&[H]]

(* Matrix A of Overscript[\[Lambda], \[RightVector]] = A.Overscript[\[Tau], \[RightVector]], with Overscript[\[Tau], \[RightVector]] the map multipliers and Overscript[\[Lambda], \[RightVector]] the corresponding Choi matrix's eigenvalues*)
A[d_List]:=KroneckerProduct@@Flatten[({#,Conjugate[#]}&[FourierMatrix[#]]&/@d),1]

(* Computes subgroup generators of Subscript[\[DoubleStruckCapitalZ], d] \[CirclePlus] Subscript[\[DoubleStruckCapitalZ], d] using Tom\[AAcute]s' theorem *)
(* Notes: currently working only for d = {d} *)
GeneratorSubgroups[d_]:=DeleteDuplicates[If[PrimePowerQ[d^2/Count[#,1/d]],IntegerDigits[Flatten[Position[N[#],N[1/d]]-1],d,2],Nothing]&/@A[{d}]];


(* Computes the Choi matrix of a Weyl map of N qudits *)
ChoiMatrix[Tau_,d_]:=
Module[{AllIndices,WeylOperatorBasis},
AllIndices=Tuples[Range[0,d-1],2];
WeylOperatorBasis=WeylMatrix[#1,#2,d]&@@@AllIndices;
Tau . Normal[KroneckerProduct[#,Conjugate[#]]&/@WeylOperatorBasis]
]

ChoiMatrix[Tau_,d_List]:=
Module[{AllIndices,WeylOperatorBasis},
AllIndices=Tuples[Tuples[Range[0,#-1],2]&/@d];
WeylOperatorBasis=WeylMatrix[#1,#2,d]&@@@AllIndices;
Tau . Normal[KroneckerProduct[#,Conjugate[#]]&/@WeylOperatorBasis]
]


(* Computes the partial transpose of the Choi matrix of a Weyl map of N qudits *)
PToverBofChoiMatrix[Tau_,d_]:=
Module[{AllIndices,WeylOperatorBasis},
AllIndices=Tuples[Range[0,d-1],2];
WeylOperatorBasis=WeylMatrix[#1,#2,d]&@@@AllIndices;
Tau . Normal[KroneckerProduct[#,ConjugateTranspose[#]]&/@WeylOperatorBasis]
]

PToverBofChoiMatrix[Tau_,d_List]:=
Module[{AllIndices,WeylOperatorBasis},
AllIndices=Transpose/@(Tuples[Tuples[Range[0,#-1],2]&/@d]);
WeylOperatorBasis=WeylMatrix[#1,#2,d]&@@@AllIndices;
Tau . Normal[KroneckerProduct[#,ConjugateTranspose[#]]&/@WeylOperatorBasis]
]


TauOfWCEFromSubgroup[Indices_,TauValues_,d_]:=Normal[Flatten[SparseArray[Thread[Indices+1->TauValues],{d,d}]]]


TausOfWCEGs[d_]:=Module[{\[Phi]},
(*homomorphisms*)
\[Phi]=Tuples[Range[0,d-1],2];DeleteDuplicates[Flatten[Table[TauOfWCEFromSubgroup[subgroup,\[Omega][d,\[Phi]i . #&/@subgroup],d],{subgroup,GeneratorSubgroups[d]},{\[Phi]i,\[Phi]}],1]]]


End[];


EndPackage[];
