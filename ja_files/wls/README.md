### Para correr en el servidor

En este directorio ejecutar:

`nohup ./<wls-file> > log/<log-file> 2>&1 &`

El comando `nohup` es para ejecutar el comando de bash y éste se siga ejecutando después de cerrar la sesión de SSH. La parte de `> log/<log-file>.log 2>&1 &` guarda el output del script en `log/<log-file>.log` 

### Cadenas

- XXZ abierto
- XXZ cerrado 
- Wisniacki abierto 
- Wisniacki cerrado 
- Wisniacki en subespacio simétrico, sólo se puede con el espín de enmedio

### XXZ abierto

L=7, 441 puntos, barriendo Jxy y varepsilon de 0 a 2 en pasos de 0.1. **Comenzó** Wed 27 Nov 2024 23:17:11 y **terminó** Thu 28 Nov 2024 11:16:33, tardó **12 horas**. 

#### Comando para correr el script en el servidor

`nohup ./temporal_haar_avg_choi_purity_xxz_open.wls > log/xxz_open_L_7_Jz_1_d_3_omega_0.log 2>&1 &`
