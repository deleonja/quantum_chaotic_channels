#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Set directory to this notebook's directory *)
SetDirectory["/home/jadeleon/Documents/quantum_chaotic_channels/ja_files"];
Get["../Mathematica_packages/QMB.wl"];


(* Par\[AAcute]metros fijos L, Subscript[J, z], \[Omega], \[Epsilon], d (sitio del defecto)*)
L = 7; 
d = 4;
{Jz,\[Omega]}={1., 0.};
indices = Tuples[Range[0,3],L-1];

(* Lista de promedios temporales del valor esperado de Haar de la pureza de Choi *)
temporalAverageChoi = {};

Do[
Print[DateString[]];

(*Hamiltoniano*)
H=LeaSpinChainHamiltonian[Jxy,Jz,\[Omega],\[Epsilon],L,d];

(*Diagonalizaci\[OAcute]n*)
{eigenvalues,eigenvectors}=Chop[Eigensystem[H]];

(*Matriz de cambio de base*)
P = Transpose[eigenvectors];

(*Operador de evoluci\[OAcute]n*)
ClearAll[U];
U[t_] := P . DiagonalMatrix[Exp[-I*eigenvalues*t]] . ConjugateTranspose[P];

(* Calcular la pureza de Choi promedio en funci\[OAcute]n del tiempo *)
choiPurity = 
ParallelTable[
	u = U[t];
	{t, Chop[Total[(1/2)^2*(1/12)^(L-1)*Purity[MatrixPartialTrace[u . (Sqrt[Power[3,(L-1)-Count[#,x_/;x!=0]]]*Pauli[Join[{0},#]]) . ConjugateTranspose[u], 1, 2]] &/@ indices]]}
	,{t, 0., 100., 0.5}
, DistributedContexts -> Full];

(* Calcular el promedio temporal de la pureza promedio de Choi *)
fInterp = Interpolation[choiPurity]; 
choiPurityTempAvg = NIntegrate[fInterp[t],{t, 0, 100}] / 100.;
AppendTo[temporalAverageChoi, choiPurityTempAvg];


(* Graficar la pureza de Choi en el tiempo *)
(*
fig=ListPlot[choiPurity,
PlotRange->{All,{-0.01,1.01}},
Joined->True,
PlotMarkers->{False,0.04},
PlotStyle->Directive[Thickness[0.004]],
Frame->True,
FrameLabel->{TraditionalForm[HoldForm[t]],TraditionalForm[HoldForm[\[DoubleStruckCapitalE][Tr[\[ScriptCapitalD]^2]]]]},
FrameStyle->Directive[Black,FontSize->18],
PlotLabel->"\!\(\*SubscriptBox[\(J\), \(z\)]\) = "<>ToString[Jz]<>", d = "<>ToString[d]<>", \[Omega] = "<>ToString[\[Omega]]<>", \[Epsilon] = "<>ToString[\[Epsilon]]<>", \!\(\*SubscriptBox[\(J\), \(xy\)]\) = "<>ToString[NumberForm[Jxy,{Infinity,2}]],
LabelStyle->Directive[Black,FontSize->16],
ImageSize->400,
Epilog->{{Opacity[0.25],Thick,Red,Line[{{-10,0.25},{200,0.25}}]},{Opacity[0.25],Thick,Red,Line[{{-10,0.5},{200,0.5}}]},{Opacity[0.35],Thick,Green,Line[{{-10,choiPurityTempAvg},{200,choiPurityTempAvg}}]}},
GridLines->Automatic,
GridLinesStyle->Directive[Gray, Dashed,Thin]
];
*)

Print["Jxy = "<>ToString[Jxy]<>", epsilon = "<>ToString[\[Epsilon]]<>": done, "<>DateString[]];

(* Exportar figura y datos *)
(*Export["figs_ja/choi/haar/Lea_L_"<>ToString[L]<>"_Jz_"<>ToString[Jz]<>"_d_"<>ToString[d]<>"_omega_"<>ToString[\[Omega]]<>"_epsilon_"<>ToString[\[Epsilon]]<>"_Jxy_"<>ToString[NumberForm[Jxy,{Infinity,2}]]<>".pdf",fig,"PDF"];*)
Export["figs_ja/choi/haar/data/Lea_L_"<>ToString[L]<>"_Jz_"<>ToString[Jz]<>"_d_"<>ToString[d]<>"_omega_"<>ToString[\[Omega]]<>"_epsilon_"<>ToString[\[Epsilon]]<>"_Jxy_"<>ToString[NumberForm[Jxy,{Infinity,2}]]<>".csv",choiPurity,"CSV"];
,{Jxy,Subdivide[0., 2., 20]}, {\[Epsilon], Subdivide[0., 2., 20]}]; 


Export["data/temporal_haar_avg_choi_purity/xxz_w_defect_L_"<>ToString[L]<>"_Jz_"<>ToString[Jz]<>"_d_"<>ToString[d]<>"_omega_"<>ToString[\[Omega]]<>"/epsilon_"<>ToString[\[Epsilon]]<>"_Jxy_"<>ToString[NumberForm[Jxy,{Infinity,2}]]<>".csv",temporalAverageChoi,"CSV"];
