%% Practica5
%%   Portada
% *Intituto Polit�cnico Nacional*
%
% *Unidad Profecional Interdiciplinaria en Ingenier�a y Tenolog�as Avanzadas*
%                         
% *An�lisis de se�ales y sistemas*
%
% *Pr�ctica 3: Se�ales en tiempo discreto.*
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
% * Manipulaci�n b�sica de MATLAB.
% * Gr�ficas de se�ales reales y complejas discretas.
% * Transformaci�n de se�ales discretas (escalamientos y traslaciones).
% * Calculo de energ�a y potencia de se�ales discretas
% 
%% Ejercicio 6.5
% Datos para la funci�n con 5 armonicos.
d0=0.504;
dn=@(n) 0.504/(1+4*n*1j);
t0=0;
tf=pi;
f=@(t) exp(-t/2);
armo=4;
a=-7;
b=9;
% Gr�fica de la funcion Delta con 5 armonicos.
sfc(t0,tf,dn,d0,f,armo,a,b)
% Datos para la funci�n con 15 armonicos.
d0=0.504;
dn=@(n) 0.504/(1+4*n*1j);
t0=0;
tf=pi;
f=@(t) exp(-t/2);
armo=4;
a=-7;
b=8;
% Gr�fica de la funcion Delta con 15 armonicos.
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
%Gr�fica de la funcion Delta con 5 armonicos.
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
%Gr�fica de la funcion Delta con 15 armonicos.
sfc(t0,tf,dn,d0,f,armo,a,b)