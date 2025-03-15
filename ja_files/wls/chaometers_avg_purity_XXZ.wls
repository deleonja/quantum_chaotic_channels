#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Set directory to this notebook's directory *)
SetDirectory["/home/jadeleon/Documents/chaos_meets_channels/ja_files"];
Get["../Mathematica_packages/QMB.wl"];
Print["Paquetes cargados. ", DateString[],"."];


(* Par\[AAcute]metros fijos L, Jz, \[Omega], d *)
L = 7; d = 3;
tmax = 100;
ensembleSize = 50;(* n\[UAcute]mero de estados iniciales sobre los cuales se promedia *)
kernels = 24;
(*******)
Jz = 1.;
\[Omega] = 0.;


(* Lista de tiempos por cada iteraci\[OAcute]n del Do *)
times = {};

points = Tuples[{#,#}]&[Subdivide[0., 5., 29]];
lenpoints = Length[points];


(* Encender todos los kernels *)
LaunchKernels[kernels];
Print["N\[UAcute]mero de kernels disponibles: ", $KernelCount];


(* Crear archivo de datos *)
dataFilename = "data/chaometer_purity/XXZ_open_L_"<>ToString[L]<>"_d_"<>ToString[d]<>"_Jz_"<>ToString[NumberForm[Jz,{Infinity,2}]]<>"_omega_"<>ToString[NumberForm[\[Omega],{Infinity,2}]]<>"_epsilon_"<>ToString[NumberForm[\[Epsilon]d,{Infinity,2}]]<>"_Jxy_"<>ToString[NumberForm[Jxy,{Infinity,2}]]<>".csv";
dataFile = OpenAppend[dataFilename];
WriteString[dataFile, "Jz,Jxy,promedio temporal de 0 a 100 de la pureza promedio del caometro"<>"\n"];
Close[dataFile];


Do[
startTime = AbsoluteTime[];

(* hz y J *)
{\[Epsilon]d, Jxy} = points[[i]];

(* Hamiltoniano *)
H = LeaSpinChainHamiltonian[Jxy, Jz, \[Omega], \[Epsilon]d, L, d];

(*Diagonalizaci\[OAcute]n*)
{eigvals, eigvecs} = Chop[Eigensystem[H]];

(* Calcular valor promedio de la pureza *)
meanPurity = 
Mean[
	Table[
		\[Psi]0 = RandomChainProductState[L];
		
		purity = 
			Chop[
			ParallelTable[
				\[Psi] = Chop[StateEvolution[t, \[Psi]0, eigvals, eigvecs]];
				\[Rho]A = MatrixPartialTrace[Dyad[\[Psi]], 2, {2, 2^(L-1)}];
				{t, Purity[\[Rho]A]}
				, {t, 0, 100, 0.5}
			, DistributedContexts->Full]
			];
	
		(* Calcular el promedio temporal de la pureza del ca\[OAcute]metro *)
		purityInterp = Interpolation[purity]; 
		tmpAvgPurity = NIntegrate[purityInterp[t],{t, 0., tmax}] / tmax
	, ensembleSize]
];

(* Agregar {\[Epsilon], Jxy, meanPurity} al archivo CSV *)
dataFile = OpenAppend[dataFilename];
WriteString[dataFile, StringJoin[Riffle[ToString/@{hz, J, meanPurity}, ","]]<>"\n"];
Close[dataFile];

endTime = AbsoluteTime[];
time = endTime - startTime;
remainingIterations = (lenpoints - i);
AppendTo[times, time];
Print[DateString[], ". ",ToString[i], "/", lenpoints,": ", ToString[NumberForm[time,{Infinity,0}]], " s. Terminar\[AAcute]: ", DateString[DatePlus[Now,Quantity[remainingIterations*Mean[times]/3600,"Hours"]]]];

, {i, lenpoints}];


Print["Termin\[EAcute], ", DateString[],"."];
