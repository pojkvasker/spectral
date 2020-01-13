close all 
clear all

N = 1000; % number of data points
t = 1:N;
p = 50; %nr of realizations
k = N; %nr of time lags
sig2 = 1; % error variance
e = sig2*randn(1,N);
phi_shift1 = zeros(1,p);
phi_shift2 = zeros(1,p);
for i=1:p
    phi_shift1(i) =  0 + (2*pi-0)*rand(1); % uniformly dist. between 0 and 2pi
    phi_shift2(i) =  0 + (2*pi-0)*rand(1); % uniformly dist. between 0 and 2pi
end
y = 10*sin(0.24*pi*t + phi_shift1(1)) + 5*sin(0.26*pi*t + phi_shift1(2)) + e;
phi_y = periodogram(y);
phi_y = [fliplr(phi_y') phi_y'];
omega = -pi:(2*pi/(N-1)):pi;
omega_y = -pi:(2*pi/(length(phi_y)-1)):pi;

figure(1)
plot(omega_y,10*log10(phi_y')), title('True spectrum of y using periodogram'),xlabel('Frequency [rad]'),ylabel('Power spectral density estimate [dB]');

% true acs
p = 50; %nr of realizations
k = N; %nr of time lags
r = zeros(p,k);
alpha = [5 2.5];
for i=1:p %realizations "t"
    for j=1:k %time lags "k"
        r(i,j) = (alpha(1).^2)*(exp(1i*0.24*pi*(j-1)+phi_shift1(i))-exp(-1i*0.24*pi*(j-1)+phi_shift2(i)))+...
                 (alpha(2).^2)*(exp(1i*0.26*pi*(j-1)+phi_shift1(i))-exp(-1i*0.26*pi*(j-1)+phi_shift2(i))) + sig2*diracfunc(j,1);
    end
end

[yw4_a, yw4_sig2] = yulewalker_acs(r(1,:),4);
[yw12_a, yw12_sig2] = yulewalker_acs(r(1,:),12);
[myw44_a, myw44_gamma] = mywarma_acs(r(1,:),4,4,4); % M=n
[myw1212_a, myw1212_gamma] = mywarma_acs(r(1,:),12,12,12); % M=n

phi_yw4 = 10*log10(fftshift(armase(1,yw4_a,yw4_sig2,length(omega))));
phi_yw12 = 10*log10(fftshift(armase(1,yw12_a,yw12_sig2,length(omega))));
phi_myw44 = 10*log10(abs(fftshift(argamse(myw44_gamma,myw44_a,length(omega)))));
phi_myw1212 = 10*log10(abs(fftshift(argamse(myw1212_gamma,myw1212_a,length(omega)))));

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
subplot(2,4,7), zplane(1,myw44_a'),hold on,plot(unitcircle),title('PZ-plot MYW for ARMA(4,4)');
subplot(2,4,8), zplane(1,myw1212_a'),hold on,plot(unitcircle),title('PZ-plot MYW for ARMA(12,12)');

% c-------------------------------------------------------------------------------------------------------------------------
Nc = 64; % number of data points
t = 1:Nc;
sig2c = 0; % error variance
phi_shiftc =  0 + (2*pi-0)*rand(1,Nc); % uniformly dist. between 0 and 2pi
yc = zeros(p,Nc);
for i=1:p
    yc(i,:) = 10*sin(0.24*pi*t + phi_shift1(i)) + 5*sin(0.26*pi*t + phi_shift1(i));
end
yw4_ac = zeros(5,p);
yw4_sig2c = zeros(5,p);
yw12_ac = zeros(13,p);
yw12_sig2c = zeros(13,p);
myw44_ac = zeros(5,p);
myw44_gammac = zeros(5,p);
myw1212_ac = zeros(13,p);
myw1212_gammac = zeros(13,p);
myw44_acc = zeros(5,p);
myw44_gammacc = zeros(5,p);
myw1212_acc = zeros(13,p);
myw1212_gammacc = zeros(13,p);
phi_yw4c = zeros(length(omega),p); 
phi_yw12c = zeros(length(omega),p);
phi_myw44c = zeros(length(omega),p);
phi_myw1212c = zeros(length(omega),p);
phi_myw44cc = zeros(length(omega),p);
phi_myw1212cc = zeros(length(omega),p);
yw8_a = zeros(9,p);
yw8_sig2 = zeros(1,p);
phi_yw8 = zeros(length(omega),p);

ls4_a = zeros(5,p);
ls4_sig2 = zeros(1,p);
ls12_a = zeros(13,p);
ls12_sig2 = zeros(1,p);
ls44_a = zeros(5,p);
ls44_b = zeros(5,p);
ls44_sig2 = zeros(1,p);
ls1212_a = zeros(13,p);
ls1212_b = zeros(13,p);
ls1212_sig2 = zeros(1,p);
ls44_a2 = zeros(5,p);
ls44_b2 = zeros(5,p);
ls44_sig22 = zeros(1,p);
ls1212_a2 = zeros(13,p);
ls1212_b2 = zeros(13,p);
ls1212_sig22 = zeros(1,p);
phi_ls4 = zeros(length(omega),p); 
phi_ls12 = zeros(length(omega),p);
phi_ls44 = zeros(length(omega),p);
phi_ls44c = zeros(length(omega),p);
phi_ls1212 = zeros(length(omega),p);
phi_ls1212c = zeros(length(omega),p);
phi_ls44cc = zeros(length(omega),p);
phi_ls1212cc = zeros(length(omega),p);
ls8_a = zeros(9,p);
ls8_sig2 = zeros(1,p);
phi_ls8 = zeros(length(omega),p);

mult = 2;
ordtest = 8;
for i=1:p
    [yw4_ac(:,i), yw4_sig2c(i)] = yulewalker(yc(i,:),4);
    [yw8_a(:,i), yw8_sig2(i)] = yulewalker(yc(i,:),8); % order 8 for d
    [yw12_ac(:,i), yw12_sig2c(i)] = yulewalker(yc(i,:),12);
    [myw44_ac(:,i), myw44_gammac(:,i)] = mywarma(yc(i,:),4,4,4); % M=n
    [myw1212_ac(:,i), myw1212_gammac(:,i)] = mywarma(yc(i,:),12,12,12); % M=n
    [myw44_acc(:,i), myw44_gammacc(:,i)] = mywarma(yc(i,:),4,4,4*mult); % M=2n
    [myw1212_acc(:,i), myw1212_gammacc(:,i)] = mywarma(yc(i,:),12,12,12*mult); % M=2n
    
    phi_yw4c(:,i) = 10*log10(fftshift(armase(1,yw4_ac(:,i),yw4_sig2c(i),length(omega))));
    phi_yw8(:,i) = 10*log10(fftshift(armase(1,yw8_a(:,i),yw8_sig2(i),length(omega)))); % order 8
    phi_yw12c(:,i) = 10*log10(fftshift(armase(1,yw12_ac(:,i),yw12_sig2c(i),length(omega))));
    phi_myw44c(:,i) = 10*log10(abs(fftshift(argamse(myw44_gammac(:,i),myw44_ac(:,i),length(omega)))));
    phi_myw1212c(:,i) = 10*log10(abs(fftshift(argamse(myw1212_gammac(:,i),myw1212_ac(:,i),length(omega)))));
    phi_myw44cc(:,i) = 10*log10(abs(fftshift(argamse(myw44_gammacc(:,i),myw44_acc(:,i),length(omega)))));
    phi_myw1212cc(:,i) = 10*log10(abs(fftshift(argamse(myw1212_gammacc(:,i),myw1212_acc(:,i),length(omega)))));
    
    [ls4_a(:,i), ls4_sig2(i)] = lsar(yc(i,:),4);
    [ls8_a(:,i), ls8_sig2(i)] = lsar(yc(i,:),8); % order 8 for d
    [ls12_a(:,i), ls12_sig2(i)] = lsar(yc(i,:),12);
    [ls44_a(:,i),ls44_b(:,i), ls44_sig2(i)] = lsarma(yc(i,:),4,4,4); % m = K
    [ls1212_a(:,i),ls1212_b(:,i), ls1212_sig2(i)] = lsarma(yc(i,:),12,12,12); % m = K
    [ls44_a2(:,i),ls44_b2(:,i), ls44_sig22(i)] = lsarma(yc(i,:),4,4,4*mult); % m = 2K
    [ls1212_a2(:,i),ls1212_b2(:,i), ls1212_sig22(i)] = lsarma(yc(i,:),12,12,12*mult); % m = 2K
    
    phi_ls4(:,i) = 10*log10(fftshift(armase(1,ls4_a(:,i),ls4_sig2(i),length(omega))));
    phi_ls8(:,i) = 10*log10(fftshift(armase(1,ls8_a(:,i),ls8_sig2(i),length(omega)))); % order 8
    phi_ls12(:,i) = 10*log10(fftshift(armase(1,ls12_a(:,i),ls12_sig2(i),length(omega))));
    phi_ls44(:,i) = 10*log10(abs(fftshift(argamse(ls44_b(:,i),ls44_a(:,i),length(omega)))));
    phi_ls1212(:,i) = 10*log10(abs(fftshift(argamse(ls1212_b(:,i),ls1212_a(:,i),length(omega)))));
    phi_ls44c(:,i) = 10*log10(abs(fftshift(argamse(ls44_b2(:,i),ls44_a2(:,i),length(omega)))));
    phi_ls1212c(:,i) = 10*log10(abs(fftshift(argamse(ls1212_b2(:,i),ls1212_a2(:,i),length(omega)))));
end
% Yule-Walker
figure(3)
subplot(2,6,1), plot(omega,phi_yw4c),title({'PSD est using','YW for AR(4)'})...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,6,2), plot(omega,phi_yw12c),title({'PSD est using','YW for AR(12)'})...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,6,3), plot(omega,phi_myw44c),title({'PSD est using','MYW for ARMA(4,4), M=n'})...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,6,4), plot(omega,phi_myw1212c),title({'PSD est using','MYW for ARMA(12,12), M=n'})...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,6,5), plot(omega,phi_myw44cc),title({'PSD est using','MYW for ARMA(4,4), M=2n'})...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,6,6), plot(omega,phi_myw1212cc),title({'PSD est using','MYW for ARMA(12,12), M=2n'})...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,6,7), zplane(1,yw4_ac'),hold on,plot(unitcircle),title({'PZ-plot YW, for AR(4)'});
subplot(2,6,8), zplane(1,yw12_ac'),hold on,plot(unitcircle),title({'PZ-plot YW, for AR(12)'});
subplot(2,6,9), zplane(1,myw44_ac'),hold on,plot(unitcircle),title({'PZ-plot MYW, for ARMA(4,4)'});
subplot(2,6,10), zplane(1,myw1212_ac'),hold on,plot(unitcircle),title({'PZ-plot MYW for, ARMA(12,12)'});
subplot(2,6,11), zplane(1,myw44_acc'),hold on,plot(unitcircle),title({'PZ-plot MYW for, ARMA(4,4), M=2n'});
subplot(2,6,12), zplane(1,myw1212_acc'),hold on,plot(unitcircle),title({'PZ-plot MYW for, ARMA(12,12),M=2n'});

% Least Squares
figure(4)
subplot(2,6,1), plot(omega,phi_ls4),title({'PSD est using','LS for AR(4)'})...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,6,2), plot(omega,phi_ls12),title({'PSD est using','LS for AR(12)'})...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,6,3), plot(omega,phi_ls44),title({'PSD est using','LSAR for ARMA(4,4), m=K'})...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,6,4), plot(omega,phi_ls1212),title({'PSD est using','LSAR for ARMA(12,12), m=K'})...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,6,5), plot(omega,phi_ls44c),title({'PSD est using','LSAR for ARMA(4,4), m=2K'})...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,6,6), plot(omega,phi_ls1212c),title({'PSD est using','LSAR for ARMA(12,12), m=2K'})...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,6,7), zplane(1,ls4_a'),hold on,plot(unitcircle),title({'PZ-plot LS, for AR(4)'});
subplot(2,6,8), zplane(1,ls12_a'),hold on,plot(unitcircle),title({'PZ-plot ls, for AR(12)'});
subplot(2,6,9), zplane(1,ls44_a'),hold on,plot(unitcircle),title({'PZ-plot LSAR, for ARMA(4,4)'});
subplot(2,6,10), zplane(1,ls1212_a'),hold on,plot(unitcircle),title({'PZ-plot LSAR for, ARMA(12,12)'});
subplot(2,6,11), zplane(1,ls44_a2'),hold on,plot(unitcircle),title({'PZ-plot LSAR for, ARMA(4,4), m=2K'});
subplot(2,6,12), zplane(1,ls1212_a2'),hold on,plot(unitcircle),title({'PZ-plot LSAR for, ARMA(12,12),m=2K'});
% order 8
figure(5)
subplot(2,2,1), plot(omega,phi_yw8),title({'PSD est using','YW for AR(8)'})...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,2,2), plot(omega,phi_ls8),title({'PSD est using','LS for AR(8)'})...
    ,xlabel('Frequency [rad]'),ylabel('Spectrum [dB]');
subplot(2,2,3), zplane(1,yw8_a'),hold on,plot(unitcircle),title({'PZ-plot YW, for AR(8)'});
subplot(2,2,4), zplane(1,ls8_a'),hold on,plot(unitcircle),title({'PZ-plot LS, for AR(8)'});


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