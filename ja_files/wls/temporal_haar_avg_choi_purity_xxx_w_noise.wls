#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Fijar directorio a ja_files dentro del repositorio *)
SetDirectory["/home/jadeleon/Documents/chaos_meets_channels/ja_files"];


(* Cargar paquetes *)
Print["Empezando a cargar paquetes."]
Get["../Mathematica_packages/QMB.wl"];
Print["Paquetes cargados. ", DateString[],"."];


(* Par\[AAcute]metros fijos: L *)
L = 7;(* n\[UAcute]mero de espines *)
hz = Flatten[Table[#, 25]]&[Subdivide[0.1, 5., 39]];(*40 puntos para hz entre 0.1 y 5*)
hzLen = Length[hz];

(* Indices de las Pauli strings en la suma del valor esperado de Haar de la pureza de Choi *)
pauliIndices = Tuples[Range[0, 3], L - 1];

(* Lista para tiempos por cada iteraci\[OAcute]n del Do *)
times={};


(* Crear archivo de datos para guardar el promedio temporal del promedio de Haar de la pureza de Choi *)

(* Comentado porque ya existe el archivo *)
(*stream = OpenAppend["data/temporal_haar_avg_choi_purity/xxx_w_noise_L_"<>ToString[L]<>".csv"];
WriteString[stream, "hz,promedio temporal de 0 a 100 del valor esperado de Haar de la pureza de Choi"<>"\n"];
Close[stream];*)


(* Encender todos los kernels (en jungkook 24) *)
LaunchKernels[];
Print["Kernels iniciados. ", DateString[],"."];


(* Para etiquetar los archivos del temporal Haar average *)
n = 26;


Do[
(* Hora en la que comienza la iteraci\[OAcute]n *)
startTime = AbsoluteTime[];

(* Lista de valores aleatorios entre [-hz, hz] *)
hzRandom = RandomReal[{-hz[[i]], hz[[i]]}, L];

(* Hamiltoniano *)
H = HeisenbergXXXwNoise[hzRandom, L];

(*Diagonalizaci\[OAcute]n*)
{eigenvalues, eigenvectors} = Chop[Eigensystem[H]];

(*Matriz de cambio de base*)
P = Transpose[eigenvectors];

(*Operador de evoluci\[OAcute]n*)
ClearAll[U];
U[t_] := P . DiagonalMatrix[Exp[-I*eigenvalues*t]] . ConjugateTranspose[P];

(* Calcular valor esperado de Haar de la pureza de Choi promedio *)
haarAvgChoiPurity = 
ParallelTable[
	u = U[t];
	{t, Chop[Total[(1/2)^2*(1/12)^(L-1)*Purity[MatrixPartialTrace[u . (Sqrt[Power[3,(L-1)-Count[#,x_/;x!=0]]]*Pauli[Join[{0},#]]) . ConjugateTranspose[u], 1, 2]] &/@ pauliIndices]]}
	,{t, 0., 100., 0.5}
, DistributedContexts -> Full];

(* Exportar datos del valor esperado de Haar de la pureza de Choi *)
Export["data/haar_avg_choi_purity/xxx_w_noise/open/L_"<>ToString[L]<>"hz_"<>ToString[NumberForm[hz[[i]],{Infinity,2}]]<>"_"<>IntegerString[n, 10, 2]<>".csv", 
	Prepend[haarAvgChoiPurity, {"t", "valor esperado de Haar de la pureza de Choi"}], "CSV"];

(* Calcular el promedio temporal de la pureza promedio de Choi *)
fInterp = Interpolation[haarAvgChoiPurity]; 
haarAvgChoiPurityTempAvg = NIntegrate[fInterp[t],{t, 0, 100}] / 100.;

(* Agregar {hz, J, haarAvgChoiPurityTempAvg} al archivo CSV *)
stream = OpenAppend["data/temporal_haar_avg_choi_purity/xxx_w_noise_L_"<>ToString[L]<>".csv"];
WriteString[stream, StringJoin[Riffle[ToString/@{hz[[i]], haarAvgChoiPurityTempAvg},","]]<>"\n"];
Close[stream];

(* Imprimir cu\[AAcute]nto tiempo tard\[OAcute] la iteraci\[OAcute]n y cu\[AAcute]ndo va a terminar aprox *)
endTime = AbsoluteTime[];
time = endTime - startTime;
remainingIterations = hzLen - i;
totalIterations = hzLen;
AppendTo[times, time];
Print[DateString[], ". n = ", n, ". hz = ", hz[[i]], ". ", i, "/", totalIterations," iteraci\[OAcute]n tard\[OAcute] ", ToString[NumberForm[time, {Infinity,0}]], " s. Terminar\[AAcute] ", DateString[DatePlus[Now,Quantity[remainingIterations*Mean[times]/3600, "Hours"]]]];

(* Contador *)
If[hz[[i]]==5., n++];

,{i, hzLen}]; 


Print["Fin"];
