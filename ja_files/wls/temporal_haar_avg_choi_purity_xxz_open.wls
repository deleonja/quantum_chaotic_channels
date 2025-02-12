#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Para correr en el servidor: nohup ./temporal_haar_avg_choi_purity_xxz_open.wls > log/xxz_open_L_7_Jz_1_d_3_omega_0.log 2>&1 & *)


(* Set directory to this notebook's directory *)
SetDirectory["/home/jadeleon/Documents/chaos_meets_channels/ja_files"];
Get["../Mathematica_packages/QMB.wl"];


(* Par\[AAcute]metros fijos L, Jz, \[Omega], \[Epsilon], d (sitio del defecto)*)
L = 7; 
d = 3;
{Jz, \[Omega]}={1., 0.};

step=0.1;
points = Flatten[Table[Flatten[Transpose/@{{Range[0,i,step],ConstantArray[i,IntegerPart[i/step]+1]},{ConstantArray[i,IntegerPart[i/step]],Range[i-step,0,-step]}},1],{i,2.1,5.,step}],1];
lenpoints = Length[points];

(* Indices de las Pauli strings en la suma del valor esperado de Haar de la pureza de Choi *)
indices = Tuples[Range[0,3],L-1];

(* Lista de tiempos por cada iteraci\[OAcute]n del Do *)
times={};

Do[
startTime = AbsoluteTime[];

{\[Epsilon], Jxy} = i;

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

(* Exportar datos de la pureza de Choi *)
Export["data/haar_avg_choi_purity/xxz/L_"<>ToString[L]<>"/Jz_"<>ToString[Jz]<>"_d_"<>ToString[d]<>"_omega_"<>ToString[\[Omega]]<>"/epsilon_"<>ToString[NumberForm[\[Epsilon],{Infinity,2}]]<>"_Jxy_"<>ToString[NumberForm[Jxy,{Infinity,2}]]<>".csv", 
	Prepend[choiPurity, {"t", "valor esperado de Haar de la pureza de Choi"}], "CSV"];

(* Calcular el promedio temporal de la pureza promedio de Choi *)
fInterp = Interpolation[choiPurity]; 
choiPurityTempAvg = NIntegrate[fInterp[t],{t, 0, 100}] / 100.;

(* Agregar {hz, J, haarAvgChoiPurityTempAvg} al archivo CSV *)
stream=OpenAppend["data/temporal_haar_avg_choi_purity/xxz_w_defect_L_"<>ToString[L]<>"_Jz_"<>ToString[Jz]<>"_d_"<>ToString[d]<>"_omega_"<>ToString[\[Omega]]<>".csv"];
WriteString[stream,StringJoin[Riffle[ToString/@{\[Epsilon], Jxy, choiPurityTempAvg},","]]<>"\n"];
Close[stream];

endTime = AbsoluteTime[];
time = endTime - startTime;
remainingIterations = (lenpoints - Position[points, i][[1,1]]);
AppendTo[times, time];
Print[DateString[], ". ",ToString[Position[points, i][[1,1]]], "/", lenpoints,": ", ToString[NumberForm[time,{Infinity,0}]], " s. Terminar\[AAcute]: ", DateString[DatePlus[Now,Quantity[remainingIterations*Mean[times]/3600,"Hours"]]]];

,{i, points}]; 


Print["Termin\[EAcute], ", DateString[],"."];
