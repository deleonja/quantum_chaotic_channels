#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory["/home/jadeleon/Documents/quantum_chaotic_channels/ja_files/spin_chains/"];
Get["../../Mathematica_packages/QMB.wl"];
Get["../../Mathematica_packages/Chaometer.wl"];


myTicksXBelow=Join[Table[{10.^i,Superscript[10,i],{0.015,0}},{i,-2,2}](*major ticks*),Flatten[Table[{10^i+j*10.^(i),"",{.005,0}},{i,-2,1},{j,1,8}],1](*minor ticks*)];


myTicksXAbove=Join[Table[{10.^i,"",{0.015,0}},{i,-2,2}](*major ticks*),Flatten[Table[{10^i+j*10.^(i),"",{.005,0}},{i,-2,1},{j,1,8}],1](*minor ticks*)];


myTicksYLeft=Join[Table[{i,NumberForm[i,{3,2}],{.015,0}},{i,0.25,1,0.15}],Table[{i,"",{0.005,0}},{i,0.25,1,0.05}]];


myTicksYRight=Join[Table[{i,"",{.015,0}},{i,0.25,1,0.15}],Table[{i,"",{0.005,0}},{i,0.25,1,0.05}]];


(* ::Section:: *)
(*C\[AAcute]lculos*)


L=9;(* number of spins *)

t=Power[10,Range[-2,2,0.05]];

Do[
(* Define a list to save the superoperators *)
superoperators={};

Print["Delta = "<>ToString[\[CapitalDelta]]<>", "<>DateString[]];

(* For the chaometer being each site of the chain, compute the superoperators of the channel *)
Do[
(* Compute eigenvalues and eigenvectors of spin chain's Hamiltonian *)
{eigenvalsH,eigenvecsH}=Transpose[Sort[Transpose[Chop[Eigensystem[ClosedXXZHamiltonian[L,\[CapitalDelta]]]]]]];

Print["Ya diagonalic\[EAcute]"];

(* Set the environent's initial state as a product random state *)
SeedRandom[20349];(*set a seed to always get the same random state \[Psi]E *)
\[Psi]E={RandomChainProductState[chaometer-1],RandomChainProductState[L-chaometer]};

Print["Before"];

AppendTo[superoperators,Chop[Superoperator[#,\[Psi]E,eigenvalsH,eigenvecsH,L,chaometer]]&/@t];

Print["After"];

Print["Chaometer's site "<>ToString[chaometer]<>": done"<>", "<>DateString[]];
,{chaometer,L}];

(* Export data of superoperators *)
Export["superoperators_data/xxz/closed_L_9_Delta_"<>ToString[\[CapitalDelta]]<>".csv",superoperators,"CSV"];
Print["Export data of superoperators: done"];

(* Compute Choi matrices and their purities *)
chois=1/2*Map[Reshuffle,superoperators,{2}];
choiPurities=Chop[Map[Purity,chois,{2}]];
Print["Compute Chois matrices and their purities: done"];

(* Plot purities *)
fig=ListLogLogPlot[Transpose[{t,#}]&/@choiPurities,
PlotRange->{{0.01-0.0015,120},{0.24-0.01,1+0.1}},
Joined->True,PlotStyle->Directive[Dashed],PlotMarkers->Automatic,
Frame->True,
FrameLabel->{HoldForm[TraditionalForm[t]],HoldForm[TraditionalForm[Tr[\[ScriptCapitalD]^2]]]},
FrameStyle->Directive[Black,FontSize->22],
FrameTicks->{{myTicksYLeft,myTicksYRight},{myTicksXBelow,myTicksXAbove}},
FrameTicksStyle->Directive[Black,20],
GridLines->{Table[10^i,{i,-2,2}],Range[0.25,1,0.05]},
ImageSize->500,
PlotLabel->"XXZ Closed, "<>ToString[TraditionalForm[HoldForm[L]]]<>" = "<>ToString[L]<>", \[CapitalDelta] = "<>ToString[\[CapitalDelta]],
LabelStyle->Directive[Black,FontSize->19],
ImagePadding->{{Automatic,Automatic},{Automatic,5}},
PlotLegends->LineLegend[ToString[#]&/@Range[9],LegendLabel->Placed["Chaometer's site",Above]]
];

(* Export the corresponding figure *)
Export["../figs_ja/xxz/closed_L_9_Delta_"<>ToString[\[CapitalDelta]]<>".pdf",fig];
,{\[CapitalDelta], {0}}]
