#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Para correr en el servidor, dentro del directorio del wls: nohup ./<wls-script>.wls > log/<log-file>.log 2>&1 & *)


(* Fijar directorio a mis archivos dentro del repositorio *)
SetDirectory["/home/jadeleon/Documents/chaos_meets_channels/ja_files"];


(* Cargar paquetes *)
Get["../Mathematica_packages/QMB.wl"];


Print["Paquetes cargados. ", DateString[],"."];


(* Par\[AAcute]metros fijos L, hx *)
L = 7; 
hx = 1.;

(* Indices de las Pauli strings en la suma del valor esperado de Haar de la pureza de Choi *)
indices = Tuples[Range[0, 3], L - 1];

(* Lista de tiempos por cada iteraci\[OAcute]n del Do *)
times={};

hzAndJvalues = Tuples[{#,#}]&[Subdivide[0.05,2.95,29]];
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

(* Agregar {hz, J, haarAvgChoiPurityTempAvg} al archivo CSV *)
stream=OpenAppend["data/temporal_haar_avg_choi_purity/wisniacki_open_L_7_hx_1.csv"];
WriteString[stream,StringJoin[Riffle[ToString/@{hz, J, haarAvgChoiPurityTempAvg},","]]<>"\n"];
Close[stream];

endTime = AbsoluteTime[];
time = endTime - startTime;
remainingIterations = (lenhzAndJvalues - hzAndJ);
AppendTo[times,time];
Print[DateString[], ". ",ToString[hzAndJ], "/", lenhzAndJvalues,": ", ToString[NumberForm[time,{Infinity,0}]], " s. Terminar\[AAcute]: ", DateString[DatePlus[Now,Quantity[remainingIterations*Mean[times]/3600,"Hours"]]]];

,{hzAndJ, lenhzAndJvalues}]; 


Print["Termin\[EAcute], ", DateString[],"."];
