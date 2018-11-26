close all; clear; clc;
ar=15;
sef=1/3;
t0=3;
t=-4:4;
for i=1:ar
   sef=sef+exp(2*pi*i*t*1j/3)/3+ exp(-2*pi*i*t.*1j/3)/3;
end
for i=1:length(t)
    dn(i)=1/3;
    On(i)=0;
    if i==length(t)
       On(2)=pi;
       On(8)=-pi;
    end
end
figure
subplot(3,1,1,'Color',[0,0.7,0.9])
stem(t,sef,'LineWidth',2,'Color',[1 0 0])
grid on
legend('Serie de Fourier','Location','Best')
xlabel('Serie de Fourier.','FontWeight','bold','FontSize',16)
axis([-4.5 4.5 0 14])

subplot(3,1,2)
stem(t,dn,'LineWidth',2,'Color',[0 0 1])
grid on
legend('Espectro de magnitud','Location','Best')
xlabel('Espectro de Fase.','FontWeight','bold','FontSize',16)

subplot(3,1,3)
stem(t,On,'LineWidth',2,'Color',[0 0 1])
grid on
legend('Espectro de fase','Location','Best')
xlabel('Espectro de Magnitud.','FontWeight','bold','FontSize',16)