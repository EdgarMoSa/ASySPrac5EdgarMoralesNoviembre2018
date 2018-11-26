%%   Portada
% *Intituto Polit�cnico Nacional*
%
% *Unidad Profecional Interdiciplinaria en Ingenier�a y Tenolog�as
% Avanzadas*
%                         
% *An�lisis de se�ales y sistemas*
%
% *Pr�ctica 5: Series de Fourier.*
%
% *Grupo: 2MV1*
%
% *Profesor: Dr. Rafael Martinez Martinez*
%  
% *Integrantes:*
%
% *Calva Lima Leonardo Ashley*
%
% *Escarcega Corona Luis*
%
% *Morales Sanabria Edgar Esteban*
%% Objetivo
%
% Los objetivos son los siguientes:
% 
% Realizar gr�ficas de series de Fourier exponenciales y trigonom�tricas en tiempo continuo
%
% Manipulaci�n de instrucciones en MATLAB
%
% Calculo n�merico de los coeficientes de Fourier
%
%% Introducci�n
%
% Existen tres maneras de representar una serie de Fourier, cada una con
% sus respectivas ventajas matem�ticas respecto a qu� tipo de funci�n, �sta
% pretende aproximar. 
%
% En este reporte de pr�ctica se analizar�n las tres; las cuales son:
% 
% * *Trigonom�trica* ; Dadas las propiedades peri�dicas de la funci�n seno
% y coseno es posible construir modelos muy semejantes a otras funciones
% con una sumando una cantidad infinita de estas dos y un par de arreglos 
% en el periodo.
% Inicialmente nos valemos de un par de constantes $$a_n$ y $$b_n$
% a quienes llamaremos coeficientes de Fourier, mismos que modear�n las
% amplitudes del seno y coseno, respectivamente. Al variar estos
% coeficientes dentro de la sumatoria se logra aproximar la se�al cada vez
% m�s.
% Esta aproximaci�n cuenta con la siguiente estructura.
%  
% $$s_f(t)= a_0/2 + \sum_{n=1}^{\infty} (a_n  cos (nt) + b_n  sin (nt))$$ 
% 
% * *Trigonom�trica Compacta* ; En ocasiones es m�s �til conocer la amplitud
% y la fase en t�rminos cosinusoidales en lugar de amplitudes cosinusoidales
% y sinusoidal. Otra forma de expresar la compleja forma de la serie de 
% Fourier. Esta aproximaci�n cuenta con la siguiente estructura.
%  
% $$s_f(t)= a_0/2 + \sum_{n=1}^{\infty} (a_n  cos (w_n t-\theta_n) )$$ 
% 
% * *Exponencial* ; Esta serie a diferencia de las anteriores nos permite
% modelar expresiones tanto reales como imaginarias dadas sus
% caracter�sticas exponenciales. 
% 
% $$s_f(t)= \sum_{-\infty}^{\infty} c_n e^{jntw_0}$$ 
% 
% Donde en este caso especial los coeficientes de Fourier son llamados
% $$c_n$$; los cuales se obtienen aplicando algunos despejes de la forma
% trigonom�tricas y finalmente se calculan de la siguiente manera
%
% $$c_n = 1/T_0 \int \limits_{t_0}^{t_0 + T_0} f(x) e^{-nw_0jt}\cdot dt$$
% 
% Como es posible apreciar en esta �ltima ecuaci�n el c�lculo de $$c_n$$
% est� en funci�n del elemento a calcular. Esto es que var�a en funci�n de
% n que es un par�metro dentro de la serie. 
% Otra forma de encontrar los coeficientes de Fourier de forma num�rica es
% a trav�s de software. En el cap�tulo 6.6 del libro Continuous-Time Signal 
% Analysis se mencionan algunos pasos para llegar a un m�todo anal�tico
%
% Aplicando sobre un intervalo de muestreo T con periodo $$T_0=1$ en la 
% se�al $$x(t)$ , tenemos que $$N_0=T_0/T$ donde $$N_0$ es el 
% n�mero de muestras en un peri�do $$T_0$$. Para encontrar la relaci�n
% entre $$C_n$ y $$x(t)$ consideramos la definici�n de los coeficientes de 
% Fourier; donde al aplicarle l�mites en lugar de integrar se llega a las siguientes dos relaciones
% 
% $$N_0=T_0/ T $$
%
% $$\Omega_0=T\omega_0=2 \pi/N_0 $$
%   
% Cabe destacar en esta s�ntesis del m�todo un concepto que traducido al
% espa�ol es llamado "Error de Solape"; el cu�l es b�sicamente una
% aproximaci�n cuando el l�mite de $$T_0->0$ y se indetermina el l�mite.
% Para este caso el m�todo num�rico tiende a valores muy peque�os cercando
% al cero, pero nunca siendo cero. Ello implica que la Serie de Fourier
% tenga una ligera diferencia a la funci�n original.
%
% $$C_n=1/N_0 \sum_{k=0}^{N_0-1} x(kT) e^{-jkn\omega_0}$$
% 
% Aplicando las dos ultimas relaciones en la ecuaci�n anterior, se tiene 
% una nueva igualdad entre los valores de $$C_n$$
%   
% $$C_{n+N_0}=C_n$
%   
% Conclusiones: La igualdad anterior establece una llamada propiedad
% de periocidad que significa que m�s all� de que $$n=N_0/2$$ n representa
% aquellos valores negativos. Por ejemplo si $$N_0=32, C_17=C_-15 , C_18=C_-14,..., 
% C_31=C_-1$ Y de esta manera el ciclo se repite desde $$n=32$
% 
% Se puede aprovechar la eficiencia de la FFT o Transformada R�pida de Fourier 
% Por sus siglas en ingl�s; que es un algoritmo utilizado en matlab en el que 
% solo se necasitan rampas $$x(t)$ sobre un periodo que inicie en $$t=0$,
% tambi�n es preferible que $$N_0=2^m$
%
%%
 T_0 = pi; N_0 = 256; T = T_0/N_0; t = (0:T:T*(N_0-1))'; M = 10;
 x = exp(-t/2); x(1) = (exp(-pi/2) + 1)/2;
 D_n = fft (x)/N_0; n = -N_0/2:N_0/2-1'; 
 clf; subplot (2, 2, 1); stem(n, abs(fftshift (D_n)),'k');
 axis ([-M M -.1 .6]); xlabel('n'); ylabel('|D_n|');
 subplot (2, 2, 2);
 stem(n, angle(fftshift(D_n)),'k');
 axis([-M M -pi pi]);
 xlabel ('n'); 
 ylabel('\angle D n [rad]')
 %%
%% Ejemplo 6.1
%
% Con serie y espectro trigonometrico, no es necesario entregar el c�digo, 
% solo la aplicaci�n al problema especifico, debe de indicar la funci�n y 
% los valores de sus coeficientes (sin incluir el procedimiento)
%
%
% 1) Realizar el programa de la serie que se indica
%% 
% 2) Gr�fica de la serie de Fourier en un intervalo que muestre 5 repeticiones
% 3) Gr�fica de la se�al y la serie de Fourier para 4 armonicos
% 4) Gr�fica del error
% 5) Gr�fica de la energ�a del error
% 6) Espectro de magnitud para 4 armonicos
% 7) Espectro de fase para 4 armonicos
% 8) Todo lo anterior para 15 armonicos
%%
%
% Utilizando el software visto en clase se define la funci�n de la
% siguiente manera:
%
%%
d0=0.504;
dn=@(n) 0.504/(1+4*n*j);
t0=0;
tf=pi;
f=@(t) exp(-t/2);
armo=4;
a=-7;
b=7;
sfcon(t0,tf,dn,d0,f,armo,a,b)
%% Ejemplo 6.2
% 
% Con serie y espectro exponencial y A=3, no es necesario entregar el c�digo,
% solo la aplicaci�n al problema especifico, debe de indicar la funci�n y los
% valores de sus coeficientes (sin incluir el procedimiento)
%
%%
%
d0=0.504;
dn=@(n) 0.504/(1+4*n*j);
t0=0;
tf=pi;
f=@(t) (2-2*t).*(t>=1/2&t<3/2);
armo=3;
a=-7;
b=7;
sfcon(t0,tf,dn,d0,f,armo,a,b)
%%
%% Ejemplo 6.4
%
% Con serie y espectro exponencial, no es necesario entregar el c�digo, solo 
% la aplicaci�n al problema especifico, debe de indicar la funci�n y los valores
% de sus coeficientes (sin incluir el procedimiento)
%
%%
% 
d0=0.504;
dn=@(n) 0.504/(1+4*n*j);
t0=0;
tf=pi;
f=@(t) 1.*(t>=-2*pi&t<2*pi);
armo=4;
a=-7;
b=7;
sfcon(t0,tf,dn,d0,f,armo,a,b)
armo=15;
sfcon(t0,tf,dn,d0,f,armo,a,b)
%%
%% Ejemplo 6.5
%
%% Ejemplo 6.7
% 
%% COMPUTER EXAMPLE C6.2 
%
% Elabore un c�digo similar al COMPUTER EXAMPLE C6.2 que se encuentra al final de la secci�n 6.2 de Lathi para el Ejempo 6.2 con los datos indicados anteriormente (no utilice inline)
%% Trapecio Compuesto
%
% Elabore un c�digo que implemente el algoritmo de trapecio compuesto para
% $n=15$, Utilice este c�digo para aproximar $D_0,...,D_4$ del ejemplo de 
% la pr�ctica. Ahora implemente el c�digo COMPUTER EXAMPLE C6.4 que se 
% encuentra al final de la secci�n 6.6 de Lathi, y calcule nuevamente el 
% los coeficientes $D_0,...,D_4$ del ejemplo propuesto. Muestre una tabla 
% que contenga los coeficientes mencionados calculados con los dos algorit-
% mos y de forma exacta, �Qu� algortmo aproxima mejor a los coeficientes?,
% para esto compare los coefientes con el valor absoluto de la resta.
%% Referencias
% 
% Dr. Rafael Martinez Martinez. (2016). Se ?nales y Sistemas. 2018, de creativecommons.orgSitio web: http://rafneta.github.io/
% 
% B P Lathi. (2005). Continuous-Time Signal Analysis-The Fourier Series. En Linear Systems and Signals(448-449). New York: Oxford University Press, Inc. .