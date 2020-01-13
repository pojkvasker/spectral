close all
clear all

phi1=0;
phi2=0;
N=64;
t=0:N-1;
e=randn(1);

y=10*sin(0.2*2*pi*t+phi1) + 5*sin((0.2+1/N)*2*pi*t+phi2) + e;


p1=periodogram(y);
p2=periodogram([y,zeros(1,N)]);
p3=periodogram([y,zeros(1,3*N)]);
p4=periodogram([y,zeros(1,5*N)]);
p5=periodogram([y,zeros(1,7*N)]);

subplot(3,2,1), plot(p1), title('Periodogram of y')
subplot(3,2,2), plot(p2), title('Periodogram of y zero-padded with N')
subplot(3,2,3), plot(p3), title('Periodogram of y zero-padded with 3N')
subplot(3,2,4), plot(p4), title('Periodogram of y zero-padded with 5N')
subplot(3,2,5), plot(p5), title('Periodogram of y zero-padded with 7N')