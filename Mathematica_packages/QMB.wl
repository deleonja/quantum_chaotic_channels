(* ::Package:: *)

(* If ForScience paclet not installed, install it. See https://github.com/MMA-ForScience/ForScience *)
If[Length[PacletFind["ForScience"]]==0, PacletInstall[FileNameJoin[{DirectoryName[$InputFileName], "ForScience-0.88.45.paclet"}]]];


(* ::Section:: *)
(*Begin package*)


BeginPackage["QMB`"];


(* For nice formatting of usage messages, see https://github.com/MMA-ForScience/ForScience *)
<<ForScience`;


(* ::Section:: *)
(*Notas*)


(* ::Text:: *)
(*Hay cosas en IAA_model.nb que 1) hay que migrar para ac\[AAcute] y 2) que hay que revisar si deber\[IAcute]a de poner ac\[AAcute]*)


(* ::Text:: *)
(*Hay cosas en los cuadernos del caometro donde hay rutinas para la secci\[OAcute]n de quantum chaos, como el unfolding etc*)


(* ::Text:: *)
(*Hay cosas de Heisenberg meets fuzzy que tambi\[EAcute]n tengo que pasar para ac\[AAcute]*)


(* ::Section:: *)
(*Usage definitions*)


(* ::Subsection::Closed:: *)
(*General quantum mechanics*)


(* All usage messages are evaluated quietly as FormatUsage[] requires FrontEnd. Therefore, if 
   QMB.wl is loaded in a .wls no error about FrontEnd pops up. *)


Quiet[
DensityMatrix::usage = FormatUsage["DensityMatrix[\[Psi]] returns the density matrix of state vector ```\[Psi]```."];
, {FrontEndObject::notavail, First::normal}];


Quiet[
Pauli::usage = FormatUsage["Pauli[0-3] gives the Pauli matrices. 
Pauli[{i_1,...,i_N}] returns the ```N```-qubit Pauli string '''Pauli[```i_1```]''' \[CircleTimes] '''...''' \[CircleTimes] '''Pauli[```i_N```]'''."];
, {FrontEndObject::notavail, First::normal}];


MatrixPartialTrace::usage = "MatrixPartialTrace[mat, n, d] calculates the partial trace of mat over the nth subspace, where all subspaces have dimension d.
MatrixPartialTrace[mat, n, {\!\(\*SubscriptBox[\(d\), \(1\)]\),\!\(\*SubscriptBox[\(d\), \(1\)]\),\[Ellipsis]}] calculates the partial trace of matrix mat over the nth subspace, where mat is assumed to lie in a space constructed as a tensor product of subspaces with dimensions {d1,d2,\[Ellipsis]}.";


VectorFromKetInComputationalBasis::usage = "VectorFromKetInComputationalBasis[ket] returns the matrix representation of ket.";


KetInComputationalBasisFromVector::usage = "KetInComputationalBasisFromVector[vector] returns the ket representation in computational basis of vector.";


Quiet[
RandomQubitState::usage = FormatUsage["RandomQubitState[] returns a Haar random qubit state."];
, {FrontEndObject::notavail, First::normal}];


Quiet[
RandomChainProductState::usage = "RandomChainProductState[L] returns a random ```L```-qubit product state.";
, {FrontEndObject::notavail, First::normal}];


Quiet[
Dyad::usage = FormatUsage["Dyad[\[Psi]] returns | \[Psi] \[RightAngleBracket]\[LeftAngleBracket] \[Psi] |.
Dyad[\[Psi],\[Phi]] returns | \[Psi] \[RightAngleBracket]\[LeftAngleBracket] \[Phi] |."];
, {FrontEndObject::notavail, First::normal}];


Quiet[
Commutator::usage = FormatUsage["Commutator[A,B] returns AB - BA."];
, {FrontEndObject::notavail, First::normal}];


Quiet[
CommutationQ::usage = FormatUsage["CommutationQ[A,B] yields True if ```A``` and ```B``` commute, and False otherwise."];
, {FrontEndObject::notavail, First::normal}];


Quiet[
MutuallyCommutingSetQ::usage=FormatUsage["MutuallyCommutingSetQ[{A,B,...}] yields True if all matrices ```{A,B,...}``` mutually commute, and False otherwise."];
, {FrontEndObject::notavail, First::normal}];


Braket::usage = "Braket[a,b] gives \!\(\*TemplateBox[{RowBox[{\"a\", \" \"}], RowBox[{\" \", \"b\"}]},\n\"BraKet\"]\).";


FixCkForStateEvoultion::usage = "FixCkForStateEvoultion[\!\(\*SubscriptBox[\(\[Psi]\), \(0\)]\), { \!\(\*TemplateBox[{SubscriptBox[\"E\", \"k\"]},\n\"Ket\"]\) }] fixes \!\(\*SubscriptBox[\(c\), \(k\)]\) = \!\(\*TemplateBox[{RowBox[{SubscriptBox[\"E\", \"k\"], \" \"}], RowBox[{\" \", SubscriptBox[\"\[Psi]\", \"0\"]}]},\n\"BraKet\"]\) for StateEvolution[]";


StateEvolution::usage = "StateEvolution[t, \!\(\*SubscriptBox[\(\[Psi]\), \(0\)]\), {E_i}, {\!\(\*TemplateBox[{\"E_i\"},\n\"Ket\"]\)} ] returns \!\(\*TemplateBox[{RowBox[{\"\[Psi]\", RowBox[{\"(\", \"t\", \")\"}]}]},\n\"Ket\"]\) = \!\(\*SubscriptBox[\(\[Sum]\), \(\(\\ \)\(i\)\)]\) \!\(\*SuperscriptBox[\(\[ExponentialE]\), \(\(-\[ImaginaryI]\)\\  \*SubscriptBox[\(E\), \(i\)]\\  t\)]\)\!\(\*TemplateBox[{SubscriptBox[\"E\", \"i\"], RowBox[{\" \", SubscriptBox[\"\[Psi]\", \"0\"]}]},\n\"BraKet\"]\)\!\(\*TemplateBox[{SubscriptBox[\"E\", \"i\"]},\n\"Ket\"]\).
StateEvolution[t, {\!\(\*SubscriptBox[\(E\), \(k\)]\)}] calculates \!\(\*TemplateBox[{RowBox[{\"\[Psi]\", RowBox[{\"(\", \"t\", \")\"}]}]},\"Ket\"]\) = \!\(\*SubscriptBox[\(\[Sum]\), \(\(\\\\\)\(i\)\)]\) \!\(\*SuperscriptBox[\(\[ExponentialE]\), \(\(-\[ImaginaryI]\)\\\\\*SubscriptBox[\(E\), \(i\)]\\\\t\)]\)\!\(\*TemplateBox[{SubscriptBox[\"E\", \"i\"], RowBox[{\" \", SubscriptBox[\"\[Psi]\", \"0\"]}]},\"BraKet\"]\)\!\(\*TemplateBox[{SubscriptBox[\"E\", \"i\"]},\"Ket\"]\) having fixed the \!\(\*SubscriptBox[\(c\), \(k\)]\)'s with FixCkForStateEvoultion[\!\(\*SubscriptBox[\(\[Psi]\), \(0\)]\), { \!\(\*TemplateBox[{SubscriptBox[\"E\", \"k\"]},\n\"Ket\"]\) }].";


Quiet[
BlochVector::usage = FormatUsage["BlochVector[\[Rho]] returns the Bloch vector of a single-qubit density matrix \[Rho]."];
, {FrontEndObject::notavail, First::normal}];


KroneckerVectorProduct::usage = "KroneckerVectorProduct[a,b] calculates \!\(\*TemplateBox[{\"a\"},\n\"Ket\"]\)\[CircleTimes]\!\(\*TemplateBox[{\"b\"},\n\"Ket\"]\).";


Purity::usage = "Purity[\[Rho]] calculates the purity of \[Rho].";


qubit::usage = "Generates a state with the parametrization of the Bloch sphere (\[Theta],\[Phi])";


coherentstate::usage = "coherentstate[state,L] Generates a spin coherent state of L spins given a general single qubit state";


(* ::Subsection::Closed:: *)
(*Quantum chaos*)


(*buscar la rutina del unfolding para meterla aqu\[IAcute]. Quiz\[AAcute]s tambi\[EAcute]n las cosas de wigner dyson y poisson*)


MeanLevelSpacingRatio::usage = "MeanLevelSpacingRatio[\!\(\*
StyleBox[\"eigenvalues\",\nFontSlant->\"Italic\"]\)] gives \[LeftAngleBracket]\!\(\*SubscriptBox[\(r\), \(n\)]\)\[RightAngleBracket] of \!\(\*
StyleBox[\"eigenvalues\",\nFontSlant->\"Italic\"]\).";


Quiet[
	IPR::usage = 
		FormatUsage["IPR[\[Psi]] computes the Inverse Participation Ratio of  ```\[Psi]``` in computational basis."];
, {FrontEndObject::notavail, First::normal}];


(* ::Subsection::Closed:: *)
(*Quantum channels*)


Reshuffle::usage = "Reshuffle[m] applies the reshuffle transformation to the matrix m with dimension \!\(\*SuperscriptBox[\(d\), \(2\)]\)\[Times]\!\(\*SuperscriptBox[\(d\), \(2\)]\).
Reshuffle[A,m,n] reshuffles matrix A, where dim(A) = mn.";


(* ::Subsection::Closed:: *)
(*Bose-Hubbard*)


BoseHubbardHamiltonian::usage = "BoseHubbardHamiltonian[N, L, J, U] returns the BH Hamiltonian for N bosons and L sites with hopping parameter J, and interaction parameter U.";


BosonEscapeKrausOperators::usage = "BosonEscapeKrausOperators[N, L]: bosons escape to nearest neighbouring sites. N: bosons, L: site.";


BosonEscapeKrausOperators2::usage = "sdfa";


HilbertSpaceDim::usage = "HilbertSpaceDim[N, L] returns the dimension of Hilbert space of a Bose Hubbard system of N bosons and L sites.";


FockBasis::usage = "FockBasis[N, L] returns the lexicographical-sorted Fock basis of N bosons and L sites.";


SortFockBasis::usage = "SortFockBasis[fockBasis] returns fockBasis in ascending-order according to the tag of Fock states.";


Tag::usage = "Tag[ { \!\(\*SubscriptBox[\(k\), \(1\)]\),\!\(\*SubscriptBox[\(k\), \(2\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(k\), \(L\)]\) } ] returns the tag \!\(\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(L\)]\)\!\(\*SqrtBox[\(100  i + 3\)]\)\!\(\*SubscriptBox[\(k\), \(i\)]\) of Fock state \!\(\*TemplateBox[{RowBox[{SubscriptBox[\"k\", \"1\"], \",\", SubscriptBox[\"k\", \"2\"], \",\", \"\[Ellipsis]\", \",\", SubscriptBox[\"k\", \"L\"]}]},\n\"Ket\"]\).";


(*Bosons*)


FockBasisStateAsColumnVector::usage = "FockBasisStateAsColumnVector[FockState, N, L] returns matrix representation of fockState={\!\(\*SubscriptBox[\(i\), \(1\)]\),\!\(\*SubscriptBox[\(i\), \(2\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(i\), \(L\)]\)}. N: bosons, L: sites.";


FockBasisIndex::usage = "FockBasisIndex[fockState, sortedTagsFockBasis] returns the position of fockState in the tag-sorted Fock basis with tags sortedTagsFockBasis.";


RenyiEntropy::usage = "RenyiEntropy[\[Alpha], \[Rho]] computes the \[Alpha]-th order Renyi entropy of density matrix \[Rho].";


BosonicPartialTrace::usage = "BosonicPartialTrace[\[Rho]] calculates the partial trace of \[Rho]. Requires initialization.";


InitializationBosonicPartialTrace::usage = "InitializationBosonicPartialTrace[{\!\(\*SubscriptBox[\(i\), \(1\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(i\), \(k\)]\)}, N, L] initializes variables for BosonicPartialTrace[] to calculate the reduced density matrix of sites {\!\(\*SubscriptBox[\(i\), \(1\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(i\), \(k\)]\)}.";


(* ::Subsubsection:: *)
(*Fuzzy measurements in bosonic systems*)


InitializeVariables::usage = "InitializeVariables[n, L, boundaries, FMmodel] sets up the necessary variables for correct running of FuzzyMeasurement[\[Psi], \!\(\*SubscriptBox[\(p\), \(fuzzy\)]\)]; boundaries: 'open' or 'closed'; FMmodel: '#NN'.";


FuzzyMeasurement::usage = "FuzzyMeasurement[\[Psi], \!\(\*SubscriptBox[\(p\), \(fuzzy\)]\)] gives \[ScriptCapitalF](\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Ket\"]\)\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Bra\"]\)) = (1 - \!\(\*SubscriptBox[\(p\), \(fuzzy\)]\))\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Ket\"]\)\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Bra\"]\) + \!\(\*SubscriptBox[\(p\), \(fuzzy\)]\) \!\(\*UnderscriptBox[\(\[Sum]\), \(i\)]\) \!\(\*SubscriptBox[\(S\), \(i\)]\)\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Ket\"]\)\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Bra\"]\)\!\(\*SubsuperscriptBox[\(S\), \(i\), \(\[Dagger]\)]\), where \!\(\*SubscriptBox[\(S\), \(i\)]\) must be initizalized runnning InitializeVariables[n, L, boundaries, FMmodel].";


(* ::Subsection:: *)
(*Spin chains*)


(* ::Subsubsection::Closed:: *)
(*Symmetries*)


SpinParityEigenvectors::usage = "SpinParityEigenvectors[L] gives a list of {even, odd} eigenvectors of the L-spin system parity operator P; P\!\(\*TemplateBox[{RowBox[{SubscriptBox[\"k\", \"1\"], \",\", \"\[Ellipsis]\", \",\", SubscriptBox[\"k\", \"L\"]}]},\n\"Ket\"]\) = \!\(\*TemplateBox[{RowBox[{SubscriptBox[\"k\", \"L\"], \",\", \"\[Ellipsis]\", \",\", SubscriptBox[\"k\", \"1\"]}]},\n\"Ket\"]\), \!\(\*SubscriptBox[\(k\), \(i\)]\)=0,1.";


(* ::Subsubsection::Closed:: *)
(*Hamiltonians*)


Quiet[
IsingNNOpenHamiltonian::usage = FormatUsage["IsingNNOpenHamiltonian[h_x,h_z,J,L] returns the Hamiltonian H = \[Sum]_{*i=1*}^L (```h_x```\[Sigma]_i^x + ```h_z```\[Sigma]_i^z) - ```J```\[Sum]_{*i=1*}^{*L-1*} \[Sigma]^z_i \[Sigma]^z_{*i+1*}.
IsingNNOpenHamiltonian[h_x,h_z,{J_1,...,J_L},L] returns the Hamiltonian H = \[Sum]_{*i=1*}^L (```h_x```\[Sigma]_i^x + ```h_z```\[Sigma]_i^z) - \[Sum]_{*i=1*}^{*L-1*} ```J_i``` \[Sigma]^z_i \[Sigma]^z_{*i+1*}."];
, {FrontEndObject::notavail, First::normal}];


IsingNNClosedHamiltonian::usage = "IsingNNClosedHamiltonian[\!\(\*
StyleBox[SubscriptBox[\"h\", \"x\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"h\", \"z\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"J\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)] returns the Hamiltonian \!\(\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(L\)]\)(\!\(\*SubscriptBox[\(h\), \(x\)]\) \!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(x\)]\) + \!\(\*SubscriptBox[\(h\), \(z\)]\) \!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(z\)]\)) + \!\(\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), L]\) \!\(\*SubscriptBox[\(J\), \(i\)]\) \!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(z\)]\)\!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i + 1\), \(z\)]\) with \!\(\*SubscriptBox[\(\[Sigma]\), \(L + 1\)]\) = \!\(\*SubscriptBox[\(\[Sigma]\), \(1\)]\).";


ClosedXXZHamiltonian::usage = "ClosedXXZHamiltonian[\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[CapitalDelta]\",\nFontSlant->\"Italic\"]\)] returns the closed XXZ 1/2-spin chain as in appendix A.1 of Quantum 8, 1510 (2024).";


OpenXXZHamiltonian::usage= "OpenXXZHamiltonian[\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[CapitalDelta]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"h1\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"h2\",\nFontSlant->\"Italic\"]\)] returns the open XXZ 1/2-spin chain as in appendix A.2 of Quantum 8, 1510 (2024).";


Quiet[
LeaSpinChainHamiltonian::usage = FormatUsage["LeaSpinChainHamiltonian[J_{*xy*},J_z,\[Omega],\[Epsilon]_d,L,d] returns the spin-1/2 chain H = \[Sum]_{*i=1*}^{*L-1*} ```J_{*xy*}```(S^x_i S^x_{*i+1*} + S^y_i S^y_{*i+1*}) + ```J_z```S^z_i S^z_{*i+1*} + \[Sum]_{*i=1*}^{*L*} ```\[Omega]``` S^z_i + \[Epsilon]_d S^z_d. [Eq. (1) in Am. J. Phys. 80, 246\[Dash]251 (2012)]."];
, {FrontEndObject::notavail, First::normal}];


HeisenbergXXXwNoise::usage="HeisenbergXXXwNoise[hz,L] returns the Heisenberg XXX spin 1/2 chain with noise: \!\(\*FormBox[\(H\\\  = \\\ \*FractionBox[\(1\), \(4\)]\\\ \(\*SubsuperscriptBox[\(\[Sum]\), \(i = 1\), \(L - 1\)]\\\ \((\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(x\)]\\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i + 1\), \(x\)]\\\  + \\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(y\)]\\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i + 1\), \(y\)]\\\  + \\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(z\)]\\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i + 1\), \(z\)])\)\)\\\  + \\\ \*FractionBox[\(1\), \(2\)]\\\ \(\*SubsuperscriptBox[\(\[Sum]\), \(i = 1\), \(L\)]\*SubsuperscriptBox[\(h\), \(i\), \(z\)]\\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(z\)]\\\ \(\((open\\\ boundaries)\)\(.\)\)\)\),
TraditionalForm]\)";


(* ::Subsection::Closed:: *)
(*Fuzzy measurement channels*)


SwapMatrix::usage = "SwapMatrix[targetSite, wrongSite, N, L] returns the swap matrix that exchanges site targetSite with wrongSite for a system of N bosons and L sites.";


FuzzyMeasurementChannel::usage = 
"FuzzyMeasurementChannel[\[Rho], p, PermMatrices] returns (1-p) \[Rho] + \!\(\*FractionBox[\(p\), \(N - 1\)]\) \!\(\*SubscriptBox[\(\[Sum]\), \(i\)]\)\!\(\*SubscriptBox[\(S\), \(i, i + 1\)]\)\!\(\*SubscriptBox[\(\[Rho]S\), \(i, i + 1\)]\).
FuzzyMeasurementChannel[\[Rho], {\!\(\*SubscriptBox[\(p\), \(totalError\)]\), \!\(\*SubscriptBox[\(p\), \(NN\)]\), \!\(\*SubscriptBox[\(p\), \(SNN\)]\)}, {{\!\(\*SubscriptBox[\(S\), \(i, i + 1\)]\)}, {\!\(\*SubscriptBox[\(S\), \(i, i + 2\)]\)}}] returns \[ScriptCapitalE](\[Rho])=(1-\!\(\*SubscriptBox[\(p\), \(totalError\)]\))\[Rho] + \!\(\*SubscriptBox[\(p\), \(totalError\)]\)(\!\(\*FractionBox[SubscriptBox[\(p\), \(NN\)], \(L - 1\)]\)) \!\(\*SuperscriptBox[SubscriptBox[\(\[Sum]\), \(i = 1\)], \(L - 1\)]\) \!\(\*SubscriptBox[\(S\), \(i, i + 1\)]\) \[Rho] \!\(\*SubscriptBox[\(S\), \(i, i + 1\)]\) + \!\(\*SubscriptBox[\(p\), \(totalError\)]\)(\!\(\*FractionBox[SubscriptBox[\(p\), \(SNN\)], \(L - 2\)]\)) \!\(\*SuperscriptBox[SubscriptBox[\(\[Sum]\), \(i = 1\)], \(L - 2\)]\) \!\(\*SubscriptBox[\(S\), \(i, i + 2\)]\) \[Rho] \!\(\*SubscriptBox[\(S\), \(i, i + 2\)]\).";


(* ::Section:: *)
(*Beginning of Package*)


Begin["`Private`"];


(* ::Section:: *)
(*Routine definitions*)


(*no poner los nombres de funciones p\[UAcute]blicas porque se joden la definici\[OAcute]n de uso*)
ClearAll[SigmaPlusSigmaMinus,SigmaMinusSigmaPlus,SigmaPlusSigmaMinus2,SigmaMinusSigmaPlus2];


(* ::Subsection::Closed:: *)
(*General quantum mechanics*)


DensityMatrix[\[Psi]_] := Outer[Times, \[Psi], Conjugate[\[Psi]]]


Pauli[0]=Pauli[{0}]=SparseArray[{{1,0}, {0,1}}]; 
Pauli[1]=Pauli[{1}]=SparseArray[{{0,1}, {1,0}}]; 
Pauli[2]=Pauli[{2}]=SparseArray[{{0,-I},{I,0}}]; 
Pauli[3]=Pauli[{3}]=SparseArray[{{1,0}, {0,-1}}];
Pauli[Indices_List] := KroneckerProduct @@ (Pauli /@ Indices)


MatrixPartialTrace=ResourceFunction["MatrixPartialTrace"];


VectorFromKetInComputationalBasis[ket_]:=Normal[SparseArray[FromDigits[ket,2]+1->1,Power[2,Length[ket]]]]


KetInComputationalBasisFromVector[vector_]:=IntegerDigits[Position[vector,1][[1,1]]-1,2,Log[2,Length[vector]]]


RandomQubitState[] := 
Module[{x,y,z,\[Theta],\[Phi]},
	{x,y,z} = RandomPoint[Sphere[]];
	{\[Theta],\[Phi]}={ArcCos[z],Sign[y]ArcCos[x/Sqrt[x^2+y^2]]};
	{Cos[\[Theta]/2],Exp[I \[Phi]]Sin[\[Theta]/2]}
]


RandomChainProductState[0] := {1}
RandomChainProductState[1] := RandomQubitState[]
RandomChainProductState[L_] := Flatten[KroneckerProduct@@Table[RandomQubitState[],L]]


Dyad[a_]:=Outer[Times,a,Conjugate[a]]
Dyad[a_,b_]:=Outer[Times,a,Conjugate[b]]


Commutator[A_,B_]:=A . B-B . A


ZeroMatrix[d_]:=ConstantArray[0,{d,d}]


CommutationQ[A_,B_]:=Commutator[A,B]==ZeroMatrix[Length[A]]


MutuallyCommutingSetQ[ListOfMatrices_]:=Module[{SetLength=Length[ListOfMatrices]},
AllTrue[Table[CommutationQ@@ListOfMatrices[[{i,j}]],{i,SetLength-1},{j,i+1,SetLength}],TrueQ,2]
]


Braket[a_,b_]:=Conjugate[a] . b


StateEvolution[t_,psi0_List,eigenvals_List,eigenvecs_List]:=
(*|\[Psi](t)\[RightAngleBracket] = Underscript[\[Sum], k] Subscript[c, k]\[ExponentialE]^(-Subscript[\[ImaginaryI]E, k]t)|Subscript[E, k]\[RightAngleBracket], Subscript[c, k]=\[LeftAngleBracket]Subscript[E, k]\[VerticalSeparator] Subscript[\[Psi], 0]\[RightAngleBracket]*)
	Module[{ck},
		ck = N[Chop[Conjugate[eigenvecs] . psi0]];
		N[Chop[Total[ ck * Exp[-I*eigenvals*N[t]] * eigenvecs]]]
	]


FixCkForStateEvoultion[\[Psi]0_, eigenvecs_] :=
	Module[{},
		ck = N[ Chop[ Conjugate[eigenvecs] . \[Psi]0 ] ];
		Heigenvecs = eigenvecs;
	]


StateEvolution[t_,eigenvals_List]:=
(*|\[Psi](t)\[RightAngleBracket] = Underscript[\[Sum], k] Subscript[c, k]\[ExponentialE]^(-Subscript[\[ImaginaryI]E, k]t)|Subscript[E, k]\[RightAngleBracket], Subscript[c, k]=\[LeftAngleBracket]Subscript[E, k]\[VerticalSeparator] Subscript[\[Psi], 0]\[RightAngleBracket]*)
	N[Chop[Total[ ck * Exp[-I*eigenvals*N[t]] * Heigenvecs]]]


BlochVector[\[Rho]_]:=Chop[Tr[Pauli[#] . \[Rho]]&/@Range[3]]


KroneckerVectorProduct[a_,b_]:=Flatten[KroneckerProduct[a,b]]


Purity[\[Rho]_]:=Tr[\[Rho] . \[Rho]] 


qubit[\[Theta]_,\[Phi]_]:=FullSimplify[Normalize[{Cos[\[Theta]/2],Exp[I \[Phi]]Sin[\[Theta]/2]}]]


coherentstate[state_,L_]:=Flatten[KroneckerProduct@@Table[state,L]]


(* ::Subsection::Closed:: *)
(*Quantum chaos*)


MeanLevelSpacingRatio[eigenvalues_]:=Mean[Min/@Transpose[{#,1/#}]&[Ratios[Differences[Sort[eigenvalues]]]]]


IPR[\[Psi]_] := Total[\[Psi]^4]


(* ::Subsection:: *)
(*Quantum channels*)


(* ::Subsubsection:: *)
(*Reshuffle*)


Reshuffle[m_] := ArrayFlatten[ArrayFlatten/@Partition[Partition[ArrayReshape[#,{Sqrt[Dimensions[m][[1]]],Sqrt[Dimensions[m][[1]]]}]&/@m,Sqrt[Dimensions[m][[1]]]],Sqrt[Dimensions[m][[1]]]],1];


Reshuffle[A_,m_,n_] := ArrayFlatten[ArrayReshape[A, {m, n, m, n}]]


(* ::Subsection::Closed:: *)
(*Bosons*)


FockBasisStateAsColumnVector[FockBasisState_,N_,L_]:=Normal[SparseArray[Position[SortFockBasis[Normal[FockBasis[N,L]]][[2]],FockBasisState]->1,Binomial[N+L-1,L]]]


(*FockBasis[N,L] computes the Fock basis for N bosons and L sites. 
It requires the subroutines:
 * Assignationk[M,N,fockState]
* HilbertSpaceDim[N,M]*)
FockBasis[N_,M_]:=Module[{k,fockState},
k=1;
Normal[Join[
(*First lexycographical Fock state*)
{fockState=SparseArray[{1->N},{M}]},
(*Rest of Fock states*)
Table[
(* With \[Eta] the new Fock state and n the previous one, assign Subscript[\[Eta], i]=Subscript[n, i] (1<=i<=k-1), Subscript[\[Eta], k]=Subscript[n, k]-1 y Subscript[\[Eta], i]=0 (i>=k+2) *)
fockState=SparseArray[Join[Table[i->fockState[[i]],{i,k-1}],{k->fockState[[k]]-1}],{M}];
(* With \[Eta] the new Fock state and n the previous one, assign Subscript[\[Eta], k+1]=N-\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(k\)]
\*SubscriptBox[\(\[Eta]\), \(i\)]\) *)
fockState[[k+1]]=N-Total[fockState[[1;;k]]];
(* Compute next value of k *)
k=Assignationk[M,N,fockState];
fockState
,HilbertSpaceDim[N,M]-1]
]]
]


SortFockBasis[fockBasis_]:=Transpose[Sort[{Tag[#],#}&/@fockBasis]]


(* ::Subsection::Closed:: *)
(*Bose Hubbard*)


(* ::Subsubsection:: *)
(*General*)


(* Define messages for incorrect types *)
BoseHubbardHamiltonian::int = "The first argument (`1`) and second argument (`2`) are expected to be integers.";
BoseHubbardHamiltonian::real = "The third argument (`1`) and fourth argument (`2`) are expected to be real numbers.";

(*BoseHubbardH[n_Integer,M_Integer,t_Real,U_Real] computes the Bose-Hubbard Hamiltonian of n particles, M sites with hopping parameter t and interaction parameter U.*)
BoseHubbardHamiltonian[N_Integer,L_Integer,J_Real,U_Real]:=Module[{sortedBasis,sortedTags},
{sortedTags,sortedBasis}=SortFockBasis[FockBasis[N,L]];
-J KineticEnergyOfBoseHubbardHamiltonian[sortedBasis,sortedTags]+
U/2PotentialEnergyOfBoseHubbardHamiltonian[sortedBasis]]

(* Handle cases where arguments don't match the expected types *)
BoseHubbardHamiltonian[N_, L_, J_, U_] := Module[{},
  If[!IntegerQ[N] || !IntegerQ[L],
    Message[BoseHubbardHamiltonian::int, N, L];
    Return[$Failed];
  ];
  If[Head[J] =!= Real || Head[U] =!= Real,
    Message[BoseHubbardHamiltonian::real, J, U];
    Return[$Failed];
  ];
];


(*Computes the tag FockBasisElement following Tag[FockBasisElement]=\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(M\)]\(
\*SqrtBox[\(100 i + 3\)]\*
StyleBox[
RowBox[{"FockBasisElement", "[", "i", "]"}],
FontSlant->"Italic"]\)\)*)
Tag[FockBasisElement_]:=N[Round[Sum[Sqrt[100 i+3]#[[i+1]],{i,0,Length[#]-1}]&[FockBasisElement],10^-8]]
Tag[Nothing]:=Nothing


InitializationBosonicPartialTrace[SitesOfSubsystem_, n_, L_] :=
Module[{SubsystemSize,SubsystemComplementSize,SystemFockBasisTags,SystemFockBasis,SubsystemFockBasisTags,SubsystemFockBasis,SubsystemComplementFockBasis,RulesForOrderedSystemBasis,RulesForOrderedSubsystemBasis,FockIndicesInRho,FockIndicesInReducedRho},
SubsystemSize=Length[SitesOfSubsystem];
SubsystemComplementSize=L-SubsystemSize;
(*System's Fock basis*)
{SystemFockBasisTags,SystemFockBasis}=SortFockBasis[FockBasis[n,L]];
(*Subsystem's Fock basis*)
{SubsystemFockBasisTags,SubsystemFockBasis}=SortFockBasis[Flatten[Table[FockBasis[k,SubsystemSize],{k,0,n}],1]];
(*Complement subsystem's Fock basis*)
SubsystemComplementFockBasis=Map[ReplacePart[ConstantArray[_,L],Thread[Complement[Range[L],SitesOfSubsystem]->#]]&,Flatten[Table[SortFockBasis[FockBasis[k,SubsystemComplementSize]][[2]],{k,0,n}],1]];(*<<<*)
RulesForOrderedSystemBasis=Thread[Rule[SystemFockBasisTags,Range[HilbertSpaceDim[n,L]]]];
SubsystemHilbertSpaceDim=Length[SubsystemFockBasis];
RulesForOrderedSubsystemBasis=Thread[Rule[SubsystemFockBasisTags,Range[SubsystemHilbertSpaceDim]]];
FockIndicesInRho=Map[Tuples[{#,#}]&[Extract[SystemFockBasis,Position[SystemFockBasis,#]]]&,SubsystemComplementFockBasis];(*<<<*)
FockIndicesInReducedRho=Map[#[[All,All,SitesOfSubsystem]]&,FockIndicesInRho];(*<<<*)
ComputationalIndicesInRho=ReplaceAll[Map[Tag,FockIndicesInRho,{3}],RulesForOrderedSystemBasis];(*<<<*)
ComputationalIndicesInReducedRho=ReplaceAll[Map[Tag,FockIndicesInReducedRho,{3}],RulesForOrderedSubsystemBasis];(*<<<*)];


BosonicPartialTrace[Rho_] := 
	Module[{MatrixElementsOfRho,rules},
	
		MatrixElementsOfRho=Extract[Rho,#]&/@ComputationalIndicesInRho;
		rules=MapThread[Thread[Rule[#1,#2]]&,{ComputationalIndicesInReducedRho,MatrixElementsOfRho}];
		Total[Map[SparseArray[#,{SubsystemHilbertSpaceDim,SubsystemHilbertSpaceDim}]&,rules]]
	]


(* ::Subsubsection:: *)
(*Fuzzy measurements in bosonic systems*)


(* Initialization function *)
InitializeVariables[n_, L_, boundaries_, FMmodel_] := 
 Module[{basis, SwapAppliedToBasis},
  
  indices = Which[
    boundaries == "open",
    Table[ReplacePart[Range[L], {i -> Mod[i + 1, L, 1], Mod[i + 1, L, 1] -> i}], {i, L - ToExpression[StringTake[FMmodel, 1]]}],
    boundaries == "closed" || boundaries == "close",
    Print["Yet to come"]
  ];
  
  permutedBasisIndices = Table[
    basis = SortFockBasis[FockBasis[n, L]][[2]];
    SwapAppliedToBasis = #[[i]] & /@ basis;
    Flatten[Position[basis, #] & /@ SwapAppliedToBasis],
    {i, indices}
  ];
  
  numberOfPermutations = Length[permutedBasisIndices];
];


FuzzyMeasurement[\[Psi]_,pFuzzy_] := 
(1 - pFuzzy) Dyad[\[Psi]] + (pFuzzy/numberOfPermutations) * Total[ Table[ Dyad[\[Psi][[i]]], {i,permutedBasisIndices}]]


(* ::Subsubsection:: *)
(*Private routines*)


(*Hkin[sortedBasis_,sortedTags_] computes the kinetic term T of Bose hubbard model.*)
KineticEnergyOfBoseHubbardHamiltonian[sortedBasis_,sortedTags_]:=Module[{HKin,d,k,i,j,tagu,l,n,M,u},
d=Length[sortedBasis];n=Plus@@sortedBasis[[1]];M=Length[sortedBasis[[1]]];
HKin=SparseArray[ConstantArray[0,{d,d}]];
k=1;
Do[
Do[
i=ij[[1]];j=ij[[2]];
u=CreationOp[i,AnnihilationOp[j,{1,v}],n];
If[u[[1]]==0,Continue,tagu=Tag[u[[2]]];
l=SearchTagPosition[tagu,sortedTags];
HKin+=SparseArray[{k,l}->N[u[[1]]],{d,d}]];,
{ij,NearestNeighborIndices[M][[;;-2]]}];(*NearestNeighborIndices[M] para condiciones peri\[OAcute]dicas*)
k=k+1;,
{v,sortedBasis}];HKin+ConjugateTranspose[HKin]]


(* AnnihilationOp[i_,fockState_] implements the action of anihiliation operator acting over fockState, with fockState of the form {Subscript[c, i],Subscript[v, i]}, with Subscript[c, i] the constant of vector Subscript[v, i] *)
AnnihilationOp[i_,fockState_]:=If[#[[2,i]]==0,{0,Nothing},{#[[1]] Sqrt[#[[2,i]]],ReplacePart[#[[2]],i->#[[2,i]]-1]}]&[fockState]
AnnihilationOp[i_,{0,Nothing}]:={0,Nothing}


(* CreationOp[i_,fockState_] implements the action of creation operator acting over fockState, with fockState of the form {Subscript[c, i],Subscript[v, i]}, with Subscript[c, i] the constant of vector Subscript[v, i] *)
CreationOp[i_,fockState_,N_]:=If[#[[2,i]]==N,{0,Nothing},{#[[1]] Sqrt[#[[2,i]]+1],ReplacePart[#[[2]],i->#[[2,i]]+1]}]&[fockState]
CreationOp[i_,{0,Nothing},N_]:={0,Nothing}


(*SearchTagPosition[tagv_,sortedTags_] searches for the position of tagv in array sortedTags. *)
SearchTagPosition[tagv_,sortedTags_]:=FromDigits[Flatten[Position[sortedTags,tagv]]]


(*NearestNeighborIndices[M_] computes the set of indices of nearest neighbours of M sites.*)
NearestNeighborIndices[M_]:=Transpose[Join[{#},{RotateLeft[#]}]]&[Range[M]]


(*Hint[sortedBasis_] computes the interaction term T of Bose hubbard model.*)
PotentialEnergyOfBoseHubbardHamiltonian[sortedBasis_]:=Module[{d},d=Length[sortedBasis];
SparseArray[Table[{i,i},{i,d}]->Table[N[ Plus@@(#^2-#&/@v)],{v,sortedBasis}],{d,d}]]


(* ::Subsection::Closed:: *)
(*BosonEscapeKrausOperators[n, L]*)


BosonEscapeKrausOperators[n_,L_] := Module[{dimH=HilbertSpaceDim[n,L],fockBasisTags,fockBasisStates,KrausOperators},
{fockBasisTags,fockBasisStates}=SortFockBasis[FockBasis[n,L]];

KrausOperators=1/Sqrt[2(n-1)]*Join[
Table[Normal[SparseArray[Thread[
Select[Table[Module[{fockState=fockBasisStates[[k]]},
{FockBasisIndex[SigmaPlusSigmaMinus[fockState,i,L],fockBasisTags],k}],
{k,dimH}],AllTrue[#,Positive]&]->1.],{dimH,dimH}]],
{i,n-1}],
Table[Normal[SparseArray[Thread[
Select[Table[Module[{fockState=fockBasisStates[[k]]},
{FockBasisIndex[SigmaMinusSigmaPlus[fockState,i,L],fockBasisTags],k}],
{k,dimH}],AllTrue[#,Positive]&]->1.],{dimH,dimH}]],
{i,n-1}]];
Append[KrausOperators,Sqrt[IdentityMatrix[dimH]-Sum[ConjugateTranspose[K] . K,{K,KrausOperators}]]]
]


SigmaPlusSigmaMinus[fockState_,i_,L_] := 
If[fockState[[i]]==L||fockState[[i+1]]==0,
Nothing,
ReplaceAt[ReplaceAt[fockState,x_:>x+1,i],x_:>x-1,i+1]
]


SigmaMinusSigmaPlus[fockState_,i_,L_] := 
If[fockState[[i]]==0||fockState[[i+1]]==L,
Nothing,
ReplaceAt[ReplaceAt[fockState,x_:>x-1,i],x_:>x+1,i+1]
]


(*Secondary routines*)


HilbertSpaceDim[n_,L_]:=(n+L-1)!/(n!(L-1)!)


FockBasisIndex[fockState_,sortedTagsFockBasis_]:=
FromDigits[Flatten[Position[sortedTagsFockBasis,Tag[fockState]],1]]


(*Check definitions*)


Assignationk[M_,N_,n_]:=If[n[[1;;M-1]]==ConstantArray[0,M-1],M-1,FromDigits[Last[Position[Normal[n[[1;;M-1]]],x_ /;x!=0]]]]


RenyiEntropy[\[Alpha]_,\[Rho]_]:=1/(1-\[Alpha]) Log[Tr[MatrixPower[\[Rho],\[Alpha]]]]


(* ::Subsection::Closed:: *)
(*Spins*)


SpinParityEigenvectors[L_]:=Module[{tuples,nonPalindromes,palindromes},
tuples=Tuples[{0,1},L];
nonPalindromes=Select[tuples,#!=Reverse[#]&];
palindromes=Complement[tuples,nonPalindromes];
nonPalindromes=DeleteDuplicatesBy[nonPalindromes,Sort[{#,Reverse[#]}]&];
Normal[
{
Join[SparseArray[FromDigits[#,2]+1->1.,2^L]&/@palindromes,Normalize[SparseArray[{FromDigits[#,2]+1->1.,FromDigits[Reverse[#],2]+1->1.},2^L]]&/@nonPalindromes],
Normalize[SparseArray[{FromDigits[#,2]+1->-1.,FromDigits[Reverse[#],2]+1->1.},2^L]]&/@nonPalindromes
}]
]


(* ::Subsubsection::Closed:: *)
(*Spin chains*)


IsingNNOpenHamiltonian[hx_,hz_,J_,L_] := Module[{NNIndices},
	NNIndices=Normal[SparseArray[Thread[{#,#+1}->3],{L}]&/@Range[L-1]];
	N[Normal[Total[{hx*Pauli[#]+hz*Pauli[3#]&/@IdentityMatrix[L],-J*(Pauli/@NNIndices)},2]]]]


IsingNNClosedHamiltonian[hx_,hz_,J_,L_] := Module[{NNIndices},
	NNIndices=Normal[SparseArray[Thread[{#,Mod[#+1,L,1]}->3],{L}]&/@Range[L]];
	N[Normal[Total[{hx*Pauli[#]+hz*Pauli[3#]&/@IdentityMatrix[L],-J*(Pauli/@NNIndices)},2]]]]


ClosedXXZHamiltonian[L_,\[CapitalDelta]_]:=
	Module[{NNindices},
		NNindices = Normal[ SparseArray[Thread[{#, Mod[# + 1, L, 1]}->1], {L}] &/@ Range[L] ];
		N[Normal[-1/2*Total[Join[Pauli/@NNindices,Pauli/@(2NNindices),\[CapitalDelta] (Pauli[#]-IdentityMatrix[2^L])&/@(3NNindices)]]]]
	]


OpenXXZHamiltonian[L_,\[CapitalDelta]_,h1_,h2_]:=
	Module[{NNindices},
		NNindices = Normal[ SparseArray[Thread[{#, Mod[# + 1, L, 1]}->1], {L}] &/@ Range[L-1] ];
		N[Normal[-1/2*Total[Join[Pauli/@NNindices,Pauli/@(2NNindices),\[CapitalDelta] (Pauli[#]-IdentityMatrix[2^L])&/@(3NNindices)]]  
		- 1/2*(h1 Pauli[Join[{1},ConstantArray[0,L-1]]] + h2*Pauli[Join[ConstantArray[0,L-1],{1}]])+ 1/2*(h1 + h2)IdentityMatrix[2^L]]]
	]


HamiltonianNN[Jxy_,Jz_,L_]:=
	Module[{NNindices},
		NNindices = Normal[ SparseArray[Thread[{#, Mod[# + 1, L, 1]}->1], {L}] &/@ Range[L-1] ];
		N[Normal[(1/4)*Total[Join[Jxy*(Pauli/@NNindices),Jxy*(Pauli/@(2NNindices)),Jz*(Pauli[#]&/@(3NNindices))]]]]
	]

HamiltonianZ[\[Omega]_,\[Epsilon]d_,L_,d_]:=N[(1/2)*(\[Omega]*Total[Pauli/@(3*IdentityMatrix[L])]+\[Epsilon]d*Pauli[Normal[SparseArray[d->3,L]]])]

LeaSpinChainHamiltonian[Jxy_,Jz_,\[Omega]_,\[Epsilon]d_,L_,d_]:=HamiltonianNN[Jxy,Jz,L]+HamiltonianZ[\[Omega],\[Epsilon]d,L,d]


HeisenbergXXXwNoise[h_List,L_]:=
Module[{NNIndices,firstSum,secondSum},
(* \sum_{k=1}^{L-1} S_k^xS_{k+1}^x + S_k^zS_{k+1}^z + S_k^zS_{k+1}^z *)
NNIndices=Normal[SparseArray[Thread[{#,#+1}->1],{L}]&/@Range[L-1]];
firstSum=1/4*Total[Table[Pauli[i*#]&/@NNIndices,{i,3}],2];

(* \sum_{k=1}^{L} h_k^z S_k^z *)
secondSum=1/2*h . (Pauli/@DiagonalMatrix[ConstantArray[3,L]]);

firstSum+secondSum
]


(* ::Subsection::Closed:: *)
(*Fuzzy measurement channels*)


SwapMatrix[targetSite_, wrongSite_, N_, L_] := Module[
    {
        tagsOfSwappedFockBasis, (* List to hold tags of swapped Fock basis *)
        sortedTags,             (* List to hold sorted tags *)
        sortedFockBasis         (* List to hold sorted Fock basis *)
    },
    
    (* Step 1: Generate and sort the Fock basis (according to its tags) *)
    {sortedTags, sortedFockBasis} = Normal[SortFockBasis[FockBasis[N, L]]];
    
    (* Step 2: Generate tags for the swapped Fock basis *)
    tagsOfSwappedFockBasis = Tag /@ Map[SwappedFockState[#, targetSite, wrongSite] &, sortedFockBasis];
    
    (* Step 3: Create the sparse matrix corresponding to the swap matrix *)
    SparseArray[
        Table[{i, Position[sortedTags, tagsOfSwappedFockBasis[[i]]][[1, 1]]},
        {i, Length[tagsOfSwappedFockBasis]}] -> 1]
]


FuzzyMeasurementChannel[\[Rho]_, pTotalError_, SwapNNmatrices_]:=(1 - pTotalError) \[Rho] +
	pTotalError 1/Length[SwapNNmatrices] Total[# . \[Rho] . # & /@ SwapNNmatrices]


FuzzyMeasurementChannel[\[Rho]_, pErrors_List, SwapMatrices_] := Module[
	{
		{pTotalError, pNN, pSNN} = pErrors,
		{SwapNN, SwapSNN} = SwapMatrices
	},
	
	(1 - pTotalError) \[Rho] + pNN ( pTotalError/Length[SwapNN] ) Total[# . \[Rho] . # & /@ SwapNN]
	+ pNN ( pTotalError/Length[SwapSNN] ) Total[# . \[Rho] . # & /@ SwapSNN]
]


(* ::Subsubsection:: *)
(*Private functions*)


SwappedFockState[fockState_, targetSite_, wrongSite_] := 
  ReplacePart[fockState, 
    Thread[# -> Part[fockState, Reverse[#]] ] & [{targetSite, wrongSite}]
  ]


(* ::Section:: *)
(*End of Package*)


End[];


EndPackage[];
