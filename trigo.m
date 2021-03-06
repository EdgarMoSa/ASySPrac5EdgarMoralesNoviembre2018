clear; clc; close all;
function trigo()
ab=input('Ingresa el intervalo del periodo [a b]:   ');
f=input('Ingresa la funcion periodica:   ');
T0=ab(2)-ab(1);
w0=2*pi/pi;
a0=0.5042;
an=@(n) 0.504*(2/(1+16.*(n.^2)));
bn=@(n) 0.504*((8.*n)/(1+16.*(n.^2)));
arm=input('Ingresa el n�mero de armonicos:   ');
a=-2*T0; b=3*T0;
% Serie de Fourier
sf=a0;
t=a:0.001:b;

for n=1:arm
    sf=sf+an(n)*cos(n*w0*t)+bn(n)*sin(n*w0*t);
end
figure
subplot(2,2,1)
plot(t,sf,'LineWidth',2)
grid on
title('Serie de Fourier');
legend('Serie de Fourier','Location','Best')
xlabel('t','FontWeight','bold','FontSize',16)
% Serie de fourier vs funcion original
sf=a0;
t1=ab(1):0.001:ab(2);

for n=1:arm
    sf=sf+an(n)*cos(n*w0*t1)+bn(n)*sin(n*w0*t1);
end
subplot(2,2,2)
plot(t1,sf,'LineWidth',2)
grid on; hold on
plot(t1,f(t1),'LineWidth',2)
legend('Funci�n original(e^{-t/2})','Serie de Fourier ','Location','Best')
title('Serie de Fourier vs Funci�n original','FontWeight','bold','FontSize',16)
xlabel('t','FontWeight','bold','FontSize',16)
axis auto
% Error
E=f(t1)-sf;
subplot(2,2,3); plot(t1,E,'LineWidth',2)
title('Error','FontWeight','bold','FontSize',16)
xlabel('t','FontWeight','bold','FontSize',16)
axis auto; grid on;
% Energia del Error
subplot(2,2,4)
area(t1,E.^2);
legend('Energia del error','Location','Best')
title('Energ�a del Error','FontWeight','bold','FontSize',16)
xlabel('t','FontWeight','bold','FontSize',16)
axis auto; grid on
end