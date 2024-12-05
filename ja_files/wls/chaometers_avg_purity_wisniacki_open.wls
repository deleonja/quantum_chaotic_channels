#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Set directory to this notebook's directory *)
SetDirectory["/home/jadeleon/Documents/chaos_meets_channels/ja_files"];
Get["../Mathematica_packages/QMB.wl"];
Print["Paquetes cargados. ", DateString[],"."];


(* Par\[AAcute]metros fijos L, hx *)
L = 6; 
hx = 1.;


(* Lista de tiempos por cada iteraci\[OAcute]n del Do *)
times = {};

Print["Hola"];
hzAndJvalues = Tuples[{#,#}]&[Subdivide[0., 3., 30]];
lenhzAndJvalues = Length[hzAndJvalues];


(* Encender todos los kernels *)
LaunchKernels[];
Print["Kernels iniciados. ", DateString[],"."];


Do[
startTime = AbsoluteTime[];

(* hz y J *)
{hz, J} = hzAndJvalues[[hzAndJ]];

(* Hamiltoniano *)
H=IsingNNOpenHamiltonian[hx, hz, J, L];

(*Diagonalizaci\[OAcute]n*)
{eigenvalues, eigenvectors} = Chop[Eigensystem[H]];

(* Calcular valor promedio de la pureza *)
meanPurity = Mean[
Table[
	\[Psi]0 = RandomChainProductState[L];
	
	purity = Chop[
	ParallelTable[
		\[Psi] = Chop[StateEvolution[t, \[Psi]0, eigenvalues, eigenvectors]];
		\[Rho]A = MatrixPartialTrace[Dyad[\[Psi]],2,{2,2^(L-1)}];
		{t, Purity[\[Rho]A]}
		, {t, 0, 100, 0.5}
	, DistributedContexts->Full]
	];

	(* Calcular el promedio temporal de la pureza del ca\[OAcute]metro *)
	fInterp = Interpolation[purity]; 
	tmpAvgPurity = NIntegrate[fInterp[t],{t, 0, 100}] / 100.
, 50]
];

(* Agregar {hz, J, meanPurity} al archivo CSV *)
stream=OpenAppend["data/chaometer_purity/wisniacki_open_L_"<>ToString[L]<>"_hx_1.csv"];
WriteString[stream,StringJoin[Riffle[ToString/@{hz, J, meanPurity},","]]<>"\n"];
Close[stream];

endTime = AbsoluteTime[];
time = endTime - startTime;
remainingIterations = (lenhzAndJvalues - hzAndJ);
AppendTo[times,time];
Print[DateString[], ". ",ToString[hzAndJ], "/", lenhzAndJvalues,": ", ToString[NumberForm[time,{Infinity,0}]], " s. Terminar\[AAcute]: ", DateString[DatePlus[Now,Quantity[remainingIterations*Mean[times]/3600,"Hours"]]]];

, {hzAndJ, lenhzAndJvalues}];


Print["Termin\[EAcute], ", DateString[],"."];
