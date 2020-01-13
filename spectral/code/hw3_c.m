close all 
clear all

N = 64; % number of data points
t = 1:N;
p = 50; %nr of realizations
k = N; %nr of time lags
sig2 = 0; % error variance
e = sig2*randn(1,N);
phi_shift = zeros(1,p);
for i=1:p
    phi_shift(i) =  0 + (2*pi-0)*rand(1); % uniformly dist. between 0 and 2pi
end
y = 10*sin(0.24*pi*t + phi_shift(1)) + 5*sin(0.26*pi*t + phi_shift(2)) + e;
phi_y = periodogram(y);
omega = -pi:(2*pi/(N-1)):pi;
omega_y = 0:(pi/(length(phi_y)-1)):pi;

%figure(1)
%plot(omega_y,10*log10(phi_y)), title('True spectrum of y using periodogram'),xlabel('Frequency [rad]'),ylabel('Power spectral density estimate [dB]');

% true acs
p = 50; %nr of realizations
k = N; %nr of time lags
r = zeros(p,k);
alpha = [5 2.5];
for i=1:p %realizations "t"
    for j=1:k %time lags "k"
        r(i,j) = (alpha(1).^2)*(exp(1i*0.24*pi*(j-1)+phi_shift(i))-exp(-1i*0.24*pi*(j-1)+phi_shift(i)))+...
                 (alpha(2).^2)*(exp(1i*0.26*pi*(j-1)+phi_shift(i))-exp(-1i*0.26*pi*(j-1)+phi_shift(i))) + sig2*diracfunc(j,1);
    end
end

[yw4_a, yw4_sig2] = yulewalker_acs(r(1,:),4);
[yw12_a, yw12_sig2] = yulewalker_acs(r(1,:),12);
[myw44_a, myw44_gamma] = mywarma_acs(r(1,:),4,4,4); % M=n
[myw1212_a, myw1212_gamma] = mywarma_acs(r(1,:),12,12,12); % M=n

phi_yw4 = fftshift(armase(1,yw4_a,yw4_sig2,length(omega)));
phi_yw12 = fftshift(armase(1,yw12_a,yw12_sig2,length(omega)));
phi_myw44 = abs(fftshift(argamse(myw44_gamma,myw44_a,length(omega))));
phi_myw1212 = abs(fftshift(argamse(myw1212_gamma,myw1212_a,length(omega))));

figure(2)
subplot(2,4,1), plot(omega,phi_yw4),title('PSD est using YW for AR(4)')...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,4,2), plot(omega,phi_yw12),title('PSD est using YW for AR(12)')...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,4,3), plot(omega,phi_myw44),title('PSD est using MYW for ARMA(4,4)')...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,4,4), plot(omega,phi_myw1212),title('PSD est using MYW for ARMA(12,12)')...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
% only myw for arma(12,12) got peaks from both sinusoids.
% the others is provably underfitted.

sys_yw4 = tf(1,yw4_a');
sys_yw12 = tf(1,yw12_a');
sys_myw44 = tf(1,myw44_a');
sys_myw1212 = tf(1,myw1212_a');
w = 0:1/100:2*pi;
unitcircle = cos(w)+1j*sin(w);
subplot(2,4,5), zplane(1,yw4_a'),hold on,plot(unitcircle),title('PZ-plot YW for AR(4)');
subplot(2,4,6), zplane(1,yw12_a'),hold on,plot(unitcircle),title('PZ-plot YW for AR(12)');
subplot(2,4,7), zplane(myw44_gamma',myw44_a'),hold on,plot(unitcircle),title('PZ-plot MYW for ARMA(4,4)');
subplot(2,4,8), zplane(myw1212_gamma',myw1212_a'),hold on,plot(unitcircle),title('PZ-plot MYW for ARMA(12,12)');

figure(3)
zplane(1,yw4_a');
function a = lsmerge(coeff,omega,index)
    a = coeff*exp(-1j*omega*index);
end

function c = diracfunc(a,b)
    if a == b
        c = 1;
    else
        c = 0;
    end
end