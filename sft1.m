function sft1(t0,a,b,n)
 % a0=0; %Coeficiente de ao
 % an=0; %Coeficiente de an
 % w=1;
 % n=20; %Número de armónicos
  t=0:.1:10;

a0=;  % Dado por software
an=(2/t0)*((exp(t/-2)*(4*n*w*sin(n*w*t)-2*cos(n*w*t))/(4*n*2*w*2+1)));


sum=0;
 for k=1:n
   sum=sum+a0+an*cos(w*k*t)+((1/(pi*k))*((1-cos(k*pi))-(cos(k*pi)-cos(2*k*pi))))*sin(w*k*t); %Series de Fourier
 end
     plot(t,sum)
        grid on
end

%https://www.calculadora-de-integrales.com/
%https://la.mathworks.com/help/matlab/ref/integral.html
%https://soymarcelo88.wordpress.com/2010/07/27/series-trigonometricas-de-fourier-en-matlab-guide/
%http://rafneta.github.io/Prac5ASySSeriesFourierNov2018/Prac5ASySSeriesFourierNov2018.html
%https://classroom.google.com/c/MTU5NjcxMzI5OTla/t/MTQ5OTA4MDU4NTha