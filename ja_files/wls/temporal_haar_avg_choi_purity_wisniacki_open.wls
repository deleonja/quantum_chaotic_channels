#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Para correr en el servidor: nohup ./temporal_haar_avg_choi_purity_xxz_open.wls > log/xxz_open_L_7_Jz_1_d_3_omega_0.log 2>&1 & *)


(* Fijar directorio a mis archivos dentro del repositorio *)
SetDirectory["/home/jadeleon/Documents/chaos_meets_channels/ja_files"];


(* Cargar paquetes *)
Get["../Mathematica_packages/QMB.wl"];


(* Par\[AAcute]metros fijos L, hx *)
L = 7; 
hx = 1.;

(* Indices de las Pauli strings en la suma del valor esperado de Haar de la pureza de Choi *)
indices = Tuples[Range[0, 3], L - 1];

(* Lista de promedios temporales del valor esperado de Haar de la pureza de Choi *)
temporalAverageChoi = {};

(* Lista de tiempos por cada iteraci\[OAcute]n del Do *)
times={};

Do[
startTime = AbsoluteTime[];

(* Hamiltoniano *)
H=IsingNNOpenHamiltonian[hx, hz, J, L];

(*Diagonalizaci\[OAcute]n*)
{eigenvalues, eigenvectors}=Chop[Eigensystem[H]];

(*Matriz de cambio de base*)
P = Transpose[eigenvectors];

(*Operador de evoluci\[OAcute]n*)
ClearAll[U];
U[t_] := P . DiagonalMatrix[Exp[-I*eigenvalues*t]] . ConjugateTranspose[P];

(* Calcular valor promedio respecto a Haar de la pureza de Choi promedio *)
haarAvgChoiPurity = 
ParallelTable[
	u = U[t];
	{t, Chop[Total[(1/2)^2*(1/12)^(L-1)*Purity[MatrixPartialTrace[u . (Sqrt[Power[3,(L-1)-Count[#,x_/;x!=0]]]*Pauli[Join[{0},#]]) . ConjugateTranspose[u], 1, 2]] &/@ indices]]}
	,{t, 0., 100., 0.5}
, DistributedContexts -> Full];

(* Exportar datos de la pureza de Choi *)
Export["data/haar_avg_choi_purity/wisniacki/open/L_"<>ToString[L]<>"/hx_"<>ToString[hx]<>"hz_"<>ToString[NumberForm[hz,{Infinity,2}]]<>"_J_"<>ToString[NumberForm[J,{Infinity,2}]]<>".csv", 
	Prepend[haarAvgChoiPurity, {"t", "valor esperado de Haar de la pureza de Choi"}], "CSV"];

(* Calcular el promedio temporal de la pureza promedio de Choi *)
fInterp = Interpolation[haarAvgChoiPurity]; 
haarAvgChoiPurityTempAvg = NIntegrate[fInterp[t],{t, 0, 100}] / 100.;
AppendTo[temporalAverageChoi, {hz, J, haarAvgChoiPurityTempAvg}];

(* Graficar la pureza de Choi en el tiempo *)
(*
fig=ListPlot[haarAvgChoiPurity,
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
Epilog->{{Opacity[0.25],Thick,Red,Line[{{-10,0.25},{200,0.25}}]},{Opacity[0.25],Thick,Red,Line[{{-10,0.5},{200,0.5}}]},{Opacity[0.35],Thick,Green,Line[{{-10,haarAvgChoiPurityTempAvg},{200,haarAvgChoiPurityTempAvg}}]}},
GridLines->Automatic,
GridLinesStyle->Directive[Gray, Dashed,Thin]
];
*)

(* Exportar figura *)
(*Export["figs_ja/choi/haar/Lea_L_"<>ToString[L]<>"_Jz_"<>ToString[Jz]<>"_d_"<>ToString[d]<>"_omega_"<>ToString[\[Omega]]<>"_epsilon_"<>ToString[\[Epsilon]]<>"_Jxy_"<>ToString[NumberForm[Jxy,{Infinity,2}]]<>".pdf",fig,"PDF"];*)

endTime = AbsoluteTime[];
time= endTime - startTime;
remainingIterations = (100-(10*(((hz-2.1)/0.1+1)-1)+((J-2.1)/0.1+1)));
AppendTo[times,time];
Print["hz = ", ToString[NumberForm[hz,{Infinity,1}]], ", J = ", ToString[NumberForm[J,{Infinity,1}]], ": ", ToString[NumberForm[time,{Infinity,0}]], " s. Terminar\[AAcute]: ", DateString[DatePlus[Now,Quantity[remainingIterations*Mean[times]/3600,"Hours"]]]];

,{hz, Subdivide[2.1, 3., 9]}, {J, Subdivide[2.1, 3., 9]}]; 


(* Exportar el promedio temporal del valor esperado de Haar de la pureza de Choi *)
Export["data/temporal_haar_avg_choi_purity/wisniacki_open_L_"<>ToString[L]<>"_hx_"<>ToString[hx]<>"_2.csv", 
	Prepend[temporalAverageChoi, {"hz", "J", "promedio temporal de 0 a 100 del valor esperado de Haar de la pureza de Choi"}], "CSV"];
