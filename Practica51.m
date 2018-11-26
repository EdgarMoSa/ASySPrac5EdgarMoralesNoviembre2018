%% Practica5
%%   Portada
% *Intituto Politécnico Nacional*
%
% *Unidad Profecional Interdiciplinaria en Ingeniería y Tenologías Avanzadas*
%                         
% *Análisis de señales y sistemas*
%
% *Práctica 3: Señales en tiempo discreto.*
% Grupo: *2MV1*
%
% Profesor: *Dr. Rafael Martinez Martinez*
%  
% Integrantes:
%
% *Calva Lima Leonardo Ashley*
%
% *Escarcega Corona Luis*
%
% *Morales Sanabria Edgar Esteban*
%%
clear; clc; close all; 
%% Objetivos
%%
% 
% * Manipulación básica de MATLAB.
% * Gráficas de señales reales y complejas discretas.
% * Transformación de señales discretas (escalamientos y traslaciones).
% * Calculo de energía y potencia de señales discretas
% 
%% Ejercicio 6.5
% Datos para la función con 5 armonicos.
d0=0.504;
dn=@(n) 0.504/(1+4*n*1j);
t0=0;
tf=pi;
f=@(t) exp(-t/2);
armo=4;
a=-7;
b=9;
% Gráfica de la funcion Delta con 5 armonicos.
sfc(t0,tf,dn,d0,f,armo,a,b)
% Datos para la función con 15 armonicos.
d0=0.504;
dn=@(n) 0.504/(1+4*n*1j);
t0=0;
tf=pi;
f=@(t) exp(-t/2);
armo=4;
a=-7;
b=8;
% Gráfica de la funcion Delta con 15 armonicos.
sfc(t0,tf,dn,d0,f,armo,a,b)
%% Ejercicio 6.7
% Datos de la funcion Delta con 5 armonicos.
d0=0;
dn=@(n) 1/2;
t0=0;
tf=2;
f=@(t) 1;
armo=5;
a=-5;
b=5;
%Gráfica de la funcion Delta con 5 armonicos.
sfc(t0,tf,dn,d0,f,armo,a,b)

% Datos de la funcion Delta con 15 armonicos.
d0=0;
dn=@(n) 1/2;
t0=0;
tf=2;
f=@(t) 1;
armo=15;
a=-5;
b=5;
%Gráfica de la funcion Delta con 15 armonicos.
sfc(t0,tf,dn,d0,f,armo,a,b)