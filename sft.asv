function sft(T_0,f,n,a,b)
% T_0 Periodo para realizar la serie
% f función original
% n número de armonicos a utilizar en la gráfica
% a, b intevalo para realizar la grafica de la serie
t=a:0.0001:b;

 % a0=0; %Coeficiente de ao
 % an=0; %Coeficiente de an
 % w=1;
 % n=20; %Número de armónicos
 % t=0:.1:10;
n=20;
W_0=(2*pi)/T_0;

f0=@(t) f;
f1=@(t) f.*cos(W_0*n*t);
f2=@(t) f.*sin(W_0*n*t);

c=T_0/2;
d=-(T_0/2);
q=integral(f0,c,d);
q1=integral(f1,c,d);
q2=integral(f2,c,d);

a0=(1/T_0)*(q); 
an=(2/T_0)*(q1);
bn=(2/T_0)*(q2);

 sum=0;
 for k=1:n
   sum=sum+a0+an*cos(w*k*t)+((1/(pi*k))*((1-cos(k*pi))-(cos(k*pi)-cos(2*k*pi))))*sin(w*k*t); %Series de Fourier
 end
     plot(t,sum)
        grid on
end