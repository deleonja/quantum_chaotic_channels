#!/usr/bin/env wolframscript
(* ::Package:: *)

L=10;(* espines *)


(* valores para (hz, J) *)
hx=1.;
{hz0,hzf,hzStep}={0.05,2.95,0.15};
{J0,Jf,Jstep}={1.,1.,0.};
hzAndJvalues=Tuples[{Range[hz0,hzf,hzStep],Range[J0,Jf,Jstep]}];


(* tiempos para la curva de Haar *)
{t0,tf,tStep}={0.,50.,1.};


(* Nombre de archivo para guardar el promedio temporal del valor esperado de Haar de la pureza de Choi *)
temporalAvgOfHaarAvgOfChoiPurityFilename="data/temporal_haar_avg_choi_purity/wisniacki_open_L_"<>ToString[NumberForm[L,{Infinity,0}]]<>"_hx_"<>ToString[NumberForm[hx,{Infinity,2}]]<>".csv";


(* n\[UAcute]mero de kernels a usar *)
kernels=64;


Print["Hola"];


(* Fijar directorio a mis archivos dentro del repositorio *)
SetDirectory["/home/miguel.gonzalez/Documents/quantum_chaotic_channels/ja_files"];
(* Cargar paquetes *)
Get["../Mathematica_packages/QMB.wl"];


(* Encender los kernels *)
LaunchKernels[kernels];


pauliIndices=Tuples[Range[0,3],L-1];

(* Pauli strings with fixed \[Sigma]_1 = id *)
pauliStrings=Pauli[Join[{0},#]]&/@pauliIndices;

constants=(1/2)^2*(1/12)^(L-1)*Power[3.,(L-1)-Count[#,x_/;x!=0]]&/@pauliIndices;

(* Lista de tiempos por cada iteraci\[OAcute]n del Do *)
times={};

lenhzAndJvalues = Length[hzAndJvalues];

dimH=2^L;
halfDimH=dimH/2;
halfDimHPlusOne=dimH/2+1;


(* Archivo para guardar promedio temporal del valor esperado de Haar de la pureza de Choi *)
dataFile=OpenAppend[temporalAvgOfHaarAvgOfChoiPurityFilename];
WriteString[dataFile,"hz,J,promedio temporal de 0 a 100 del valor esperado de Haar de la pureza de Choi"<>"\n"];
Close[dataFile];


Do[
startTime = AbsoluteTime[];

(* hz y J *)
{hz,J}=hzAndJvalues[[hzAndJ]];

(* Crear Hamiltoniano y diagonalizar *)
H=IsingNNOpenHamiltonian[hx,hz,J,L];
{eigenvalues, eigenvectors}=Chop[Eigensystem[H]];

(* Calcular la matriz de cambio de base de la de eigenenerg\[IAcute]as a la computacional *)
P = Transpose[eigenvectors];

(* Operador de evoluci\[OAcute]n *)
ClearAll[U];U[t_]:=U[t]=P . DiagonalMatrix[Exp[-I*eigenvalues*t]] . ConjugateTranspose[P];

(* Calcular valor esperado de Haar de la pureza de Choi para cada t *)
haarAvgChoiPurity = {};
Do[
startTimeTimePoint=AbsoluteTime[];
u=U[t];uDagger=ConjugateTranspose[u];
h=Total[
ParallelTable[
(* Compute partial trace *)
A=#[[1;;halfDimH,1;;halfDimH]]+#[[halfDimHPlusOne;;dimH,halfDimHPlusOne;;dimH]]&[u . pauliStrings[[i]] . uDagger,L];
(* Purity *)
Chop[constants[[i]]*Tr[A . A]]
,{i,4^(L-1)}
,DistributedContexts->Full
,ProgressReporting->False]
];
AppendTo[haarAvgChoiPurity,{t,h}];
endTimeTimePoint=AbsoluteTime[];
Print["(", ToString[NumberForm[t,{Infinity,0}]], ", ", ToString[NumberForm[h,{Infinity,0}]], "), ", ToString[NumberForm[startTimeTimePoint-endTimeTimePoint,{Infinity,0}]], " s."];
	,{t,t0,tf,tStep}];

(* Exportar valor esperado de la pureza de Choi para (hz,J) *)
HaarAvgChoiPurityFilename="data/haar_avg_choi_purity/wisniacki/open/L_"<>ToString[NumberForm[L,{Infinity,0}]]<>"/hx_"<>ToString[NumberForm[hx,{Infinity,2}]]<>"_hz_"<>ToString[NumberForm[hz,{Infinity,2}]]<>"_J_"<>ToString[NumberForm[J,{Infinity,2}]]<>".csv";
Export[HaarAvgChoiPurityFilename,Prepend[haarAvgChoiPurity, {"t", "valor esperado de Haar de la pureza de Choi"}], "CSV"];

(* Calcular el promedio temporal de la pureza promedio de Choi *)
fInterp = Interpolation[haarAvgChoiPurity]; 
haarAvgChoiPurityTempAvg = NIntegrate[fInterp[t],{t,t0,tf}]/tf;

(* Agregar {hz, J, haarAvgChoiPurityTempAvg} al archivo CSV *)
dataFile=OpenAppend[temporalAvgOfHaarAvgOfChoiPurityFilename];
WriteString[dataFile,StringJoin[Riffle[ToString/@{hz,J,haarAvgChoiPurityTempAvg},","]]<>"\n"];
Close[dataFile];

endTime = AbsoluteTime[];
time = endTime-startTime;
remainingIterations=(lenhzAndJvalues-hzAndJ);
AppendTo[times,time];
Print[DateString[], ". ",ToString[hzAndJ], "/", lenhzAndJvalues,": ", ToString[NumberForm[time,{Infinity,0}]], " s. Terminar\[AAcute]: ", DateString[DatePlus[Now,Quantity[remainingIterations*Mean[times]/3600,"Hours"]]]];

,{hzAndJ,lenhzAndJvalues}]; 
