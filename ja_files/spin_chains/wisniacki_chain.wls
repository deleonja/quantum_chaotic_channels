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


?ClosedXXZHamiltonian


L=9;(* number of spins *)
{J,hx}={1,1};

t=Power[10,Range[-2,2,0.05]];

Do[
(* Define a list to save the superoperators *)
superoperators={};

Print["hz = "<>ToString[hz]];

(* For the chaometer being each site of the chain, compute the superoperators of the channel *)
Do[
(* Compute eigenvalues and eigenvectors of spin chain's Hamiltonian *)
{eigenvalsH,eigenvecsH}=Transpose[Sort[Transpose[Chop[Eigensystem[IsingNNOpenHamiltonian[hx,hz,J,L]]]]]];

(* Set the environent's initial state as a product random state *)
SeedRandom[20349];(*set a seed to always get the same random state \[Psi]E *)
\[Psi]E={RandomChainProductState[chaometer-1],RandomChainProductState[L-chaometer]};

AppendTo[superoperators,Chop[Superoperator[#,\[Psi]E,eigenvalsH,eigenvecsH,L,chaometer]]&/@t];

Print["Chaometer's site "<>ToString[chaometer]<>": done"];
,{chaometer,L}];

(* Export data of superoperators *)
Export["superoperators_data/open/L_9_J_hx_1_hz_"<>ToString[hz]<>".csv",superoperators,"CSV"];
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
PlotLabel->"Open, "<>ToString[TraditionalForm[HoldForm[L]]]<>" = "<>ToString[L]<>", J = \!\(\*SubscriptBox[\(h\), \(x\)]\) = "<>ToString[J]<>", \!\(\*SubscriptBox[\(h\), \(z\)]\) = "<>ToString[hz],
LabelStyle->Directive[Black,FontSize->19],
ImagePadding->{{Automatic,Automatic},{Automatic,5}},
PlotLegends->LineLegend[ToString[#]&/@Range[9],LegendLabel->Placed["Chaometer's site",Above]]
];

(* Export the corresponding figure *)
Export["../figs_ja/wisniacki/varying_chaometers_site/open_L_9_J_hx_1_hz_"<>ToString[hz]<>".pdf",fig];
,{hz, 0, 3, 0.25}]
