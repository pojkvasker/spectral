clear all
close all

% import data from current folder
load('lynxdata.mat')
load('sunspotdata.mat')

% removing mean
loglynx = loglynx - mean(loglynx);
lynx = lynx - mean(lynx);
sunspot = sunspot - mean(sunspot);

len_logl = length(loglynx);
len_lynx = length(lynx);
len_sun = length(sunspot);

orders_n = 20;
orders_m = 20;

% vector for saving variances for different orders_n
vars_logl = zeros(orders_n);
vars_lynx = zeros(orders_n);
vars_sun = zeros(orders_n);

% vectors for saving aic, aicc, gic and bic values for each order
aic_logl = zeros(1,orders_n);
aic_lynx = zeros(1,orders_n);
aic_sun = zeros(1,orders_n);
aicc_logl = zeros(1,orders_n);
aicc_lynx = zeros(1,orders_n);
aicc_sun = zeros(1,orders_n);
gic_logl = zeros(1,orders_n);
gic_lynx = zeros(1,orders_n);
gic_sun = zeros(1,orders_n);
bic_logl = zeros(1,orders_n);
bic_lynx = zeros(1,orders_n);
bic_sun = zeros(1,orders_n);

aic_logl_ = zeros(1,orders_m);
aic_lynx_ = zeros(1,orders_m);
aic_sun_ = zeros(1,orders_m);
aicc_logl_ = zeros(1,orders_m);
aicc_lynx_ = zeros(1,orders_m);
aicc_sun_ = zeros(1,orders_m);
gic_logl_ = zeros(1,orders_m);
gic_lynx_ = zeros(1,orders_m);
gic_sun_ = zeros(1,orders_m);
bic_logl_ = zeros(1,orders_m);
bic_lynx_ = zeros(1,orders_m);
bic_sun_ = zeros(1,orders_m);

% match the index of the lowest criterion with this array
% i1 corresp to aic for logl
% i12 corresp to bic to sunspot
% etc. .. ... . .
i1 = zeros(1,orders_m);
i2 = zeros(1,orders_m);
i3 = zeros(1,orders_m);
i4 = zeros(1,orders_m);
i5 = zeros(1,orders_m);
i6 = zeros(1,orders_m);
i7 = zeros(1,orders_m);
i8 = zeros(1,orders_m);
i9 = zeros(1,orders_m);
i10 = zeros(1,orders_m);
i11 = zeros(1,orders_m);
i12 = zeros(1,orders_m);

for m=1:orders_m
    for n=1:orders_n
        k = 20;
        m_lynx = 11;
        n_lynx = 4;
        m_sun = 9;
        n_sun = 3;
        m_logl = 5;
        n_logl = 2;
        
%         m_lynx = m;
%         n_lynx = n;
%         m_sun = m;
%         n_sun = n;
%         m_logl = m;
%         n_logl = n;
        
        % Least Square ARMA coeff est
        [a_params_logl_lsarma,b_params_logl_lsarma,var_logl] = lsarma(loglynx,n_logl,m_logl,k);
        [a_params_lynx_lsarma,b_params_lynx_lsarma,var_lynx] = lsarma(lynx,n_lynx,m_lynx,k);
        [a_params_sun_lsarma,b_params_sun_lsarma,var_sun] = lsarma(sunspot,n_sun,m_sun,k);
        
        phi_logl_arma = armase(b_params_logl_lsarma,a_params_logl_lsarma,var_logl,len_logl);
        phi_lynx_arma = armase(b_params_lynx_lsarma,a_params_lynx_lsarma,var_lynx,len_lynx);
        phi_sun_arma = armase(b_params_sun_lsarma,a_params_sun_lsarma,var_sun,len_sun);
        
        vars_logl(n) = var_logl;
        vars_lynx(n) = var_lynx;
        vars_sun(n) = var_sun;
        
        aic_logl(n) = aic(var_logl,n,len_logl);
        aic_lynx(n) = aic(var_lynx,n,len_lynx);
        aic_sun(n) = aic(var_sun,n,len_sun);
        aicc_logl(n) = aicc(var_logl,n,len_logl);
        aicc_lynx(n) = aicc(var_lynx,n,len_lynx);
        aicc_sun(n) = aicc(var_sun,n,len_sun);
        gic_logl(n) = gic(var_logl,n,len_logl,4);
        gic_lynx(n) = gic(var_lynx,n,len_lynx,4);
        gic_sun(n) = gic(var_sun,n,len_sun,4);
        bic_logl(n) = bic(var_logl,n,len_logl);
        bic_lynx(n) = bic(var_lynx,n,len_lynx);
        bic_sun(n) = bic(var_sun,n,len_sun);
    end
    % i corresp to n = AR order = poles
    [aic_logl_(m),i1(m)] = min(aic_logl);
    [aic_lynx_(m),i2(m)] = min(aic_lynx);
    [aic_sun_(m),i3(m)] = min(aic_sun);
    [aicc_logl_(m),i4(m)] = min(aicc_logl);
    [aicc_lynx_(m),i5(m)] = min(aicc_lynx);
    [aicc_sun_(m),i6(m)] = min(aicc_sun);
    [gic_logl_(m),i7(m)] = min(gic_logl);
    [gic_lynx_(m),i8(m)] = min(gic_lynx);
    [gic_sun_(m),i9(m)] = min(gic_sun);
    [bic_logl_(m),i10(m)] = min(gic_logl);
    [bic_lynx_(m),i11(m)] = min(gic_lynx);
    [bic_sun_(m),i12(m)] = min(gic_sun);
end

omega_logl = -pi:(2*pi/(len_logl-1)):pi;
omega_lynx = -pi:(2*pi/(len_lynx-1)):pi;
omega_sun = -pi:(2*pi/(len_sun-1)):pi;
w = 0:1/100:2*pi;
unitcircle = cos(w)+1j*sin(w);

figure(1)
subplot(2,3,1), plot(omega_logl,fftshift(10*log10(phi_logl_arma))), title('Loglynx LSARMA'),xlabel('frequency [rad]'), ylabel('spectrum [dB]');
subplot(2,3,2), plot(omega_lynx,fftshift(10*log10(phi_lynx_arma))), title('Lynx LSARMA'),xlabel('frequency [rad]'), ylabel('spectrum [dB]');
subplot(2,3,3), plot(omega_sun,fftshift(10*log10(phi_sun_arma))), title('Sunspots LSARMA'),xlabel('frequency [rad]'), ylabel('spectrum [dB]');
sys_logl_lsar = tf(b_params_logl_lsarma',a_params_logl_lsarma');
sys_lynx_lsar = tf(b_params_lynx_lsarma',a_params_lynx_lsarma');
sys_sun_lsar = tf(b_params_sun_lsarma',a_params_sun_lsarma');
subplot(2,3,4), pzplot(sys_logl_lsar), hold on, plot(unitcircle), title('PZ-plot Loglynx');
subplot(2,3,5), pzplot(sys_lynx_lsar), hold on, plot(unitcircle), title('PZ-plot Lynx');
subplot(2,3,6), pzplot(sys_sun_lsar), hold on, plot(unitcircle), title('PZ-plot Sunspots');

figure(3)
subplot(1,3,1), plot(1:orders_n,vars_logl), title('variance vs order loglynx');
subplot(1,3,2), plot(1:orders_n,vars_lynx), title('variance vs order lynx');
subplot(1,3,3), plot(1:orders_n,vars_sun), title('variance vs order sunspot');

figure(4)
subplot(2,2,1), plot(1:orders_n,aic_logl_), title('AIC vs order loglynx');
subplot(2,2,2), plot(1:orders_n,aicc_logl_), title('AICC vs order loglynx');
subplot(2,2,3), plot(1:orders_n,gic_logl_), title('GIC vs order loglynx');
subplot(2,2,4), plot(1:orders_n,bic_logl_), title('BIC vs order loglynx');

figure(5)
subplot(2,2,1), plot(1:orders_n,aic_lynx_),grid on, title('AIC vs order lynx');
subplot(2,2,2), plot(1:orders_n,aicc_lynx_),grid on, title('AICC vs order lynx');
subplot(2,2,3), plot(1:orders_n,gic_lynx_), grid on,title('GIC vs order lynx');
subplot(2,2,4), plot(1:orders_n,bic_lynx_), grid on,title('BIC vs order lynx');

figure(6)
subplot(2,2,1), plot(1:orders_n,aic_sun_), title('AIC vs order sunspot');
subplot(2,2,2), plot(1:orders_n,aicc_sun_), title('AICC vs order sunspot');
subplot(2,2,3), plot(1:orders_n,gic_sun_), title('GIC vs order sunspot');
subplot(2,2,4), plot(1:orders_n,bic_sun_), title('BIC vs order sunspot');


function b = aic(sig2,order,N)
    b = N*log(sig2) + 2*order;
end

function b = aicc(sig2,order,N)
    b = N*log(sig2) + (2*order * N)./(N-order-1);
end

function b = gic(sig2,order,N,gic_param)
% gic param usually between [2 6], hence set to 4 :)
    b = N*log(sig2) +  order * gic_param;
end

function b = bic(sig2,order,N)
    b = N*log(sig2) +  order * log(N);
end