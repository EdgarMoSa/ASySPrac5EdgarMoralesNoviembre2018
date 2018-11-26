function sff(n,T_0, f0, f1, f2)
% n Numero de armonicos
 W_0=(2*pi)/T_0;
 t=-10:0.0001:10;

 % f1 Es el producto de f.*cos(W_0*n*t);
 % f2 Es el producto de f.*sin(W_0*n*t);
 
c=T_0/2;
d=-(T_0/2);
q=@(t) integral(f0,d,c);
q1=@(t)integral(f1,d,c);
q2=@(t)integral(f2,d,c);

a0=(q); 
an=(q1);
bn=(q2);

 sum=0;
 for k=1:n
   % sum=sum+a0/2+an*cos(W_0*k*t)+bn*sin(W_0*k*t) +((1/(pi*k))*((1-cos(k*pi))-(cos(k*pi)-cos(2*k*pi))))*sin(w*k*t); %Series de Fourier
    sum = sum + a0 + an*cos(W_0*k*t) + bn*sin(W_0*k*t);
 end
     plot(t,sum)
        grid on