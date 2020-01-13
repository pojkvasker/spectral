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

orders = floor(len_logl/4);

% vector for saving variances for different orders
vars_logl = zeros(1,orders);
vars_lynx = zeros(1,orders);
vars_sun = zeros(1,orders);

% vectors for saving aic, aicc, gic and bic values for each order
aic_logl = zeros(1,orders);
aic_lynx = zeros(1,orders);
aic_sun = zeros(1,orders);
aicc_logl = zeros(1,orders);
aicc_lynx = zeros(1,orders);
aicc_sun = zeros(1,orders);
gic_logl = zeros(1,orders);
gic_lynx = zeros(1,orders);
gic_sun = zeros(1,orders);
bic_logl = zeros(1,orders);
bic_lynx = zeros(1,orders);
bic_sun = zeros(1,orders);
        
for n=1:orders
        
        % model orders
        
        N_logl_lsar = 12;
        N_lynx_lsar = 8;
        N_sun_lsar = 9;
%         N_logl_lsar = n;
%         N_lynx_lsar = n;
%         N_sun_lsar = n;
%         
        % K should not be large wrt N
        K_logl_lsarma = floor(len_logl/20);
        K_lynx_lsarma = floor(len_lynx/20);
        K_sun_lsarma = floor(len_sun/20);
        
        
        % Least Square AR coeff est
        [params_logl_lsar,var_logl_lsar] = lsar(loglynx,N_logl_lsar);
        [params_lynx_lsar,var_lynx_lsar] = lsar(lynx,N_lynx_lsar);
        [params_sun_lsar,var_sun_lsar] = lsar(sunspot,N_sun_lsar);
        
        omega_logl = -pi:(2*pi/(len_logl-1)):pi;
        omega_lynx = -pi:(2*pi/(len_lynx-1)):pi;
        omega_sun = -pi:(2*pi/(len_sun-1)):pi;
        
        est_logl_lsar = zeros(length(omega_logl),1);
        est_lynx_lsar = zeros(length(omega_lynx),1);
        est_sun_lsar = zeros(length(omega_sun),1);
        
        % Calculating spectrum manually
        
        % Loglynx
        k1 = 1:N_logl_lsar+1;
        for i=k1
            est_logl_lsar(:) = est_logl_lsar(:) + lsmerge(params_logl_lsar(i),omega_logl,i-1)';
        end
        
        % Lynx
        k2 = 1:N_lynx_lsar+1;
        for i=k2
            est_lynx_lsar(:) = est_lynx_lsar(:) + lsmerge(params_lynx_lsar(i),omega_lynx,i-1)';
        end
        
        % Sunspot
        k3 = 1:N_sun_lsar+1;
        for i=k3
            est_sun_lsar(:) = est_sun_lsar(:) + lsmerge(params_sun_lsar(i),omega_sun,i-1)';
        end
        
        phi_logl_lsar = var_logl_lsar./abs(est_logl_lsar).^2;
        phi_logl_lsar = 10*log10(phi_logl_lsar);
        
        phi_lynx_lsar = var_lynx_lsar./abs(est_lynx_lsar).^2;
        phi_lynx_lsar = 10*log10(phi_lynx_lsar);
        
        phi_sun_lsar = var_sun_lsar./abs(est_sun_lsar).^2;
        phi_sun_lsar = 10*log10(phi_sun_lsar);
        
        vars_logl(n) = var_logl_lsar;
        vars_lynx(n) = var_lynx_lsar;
        vars_sun(n) = var_sun_lsar;
        
        aic_logl(n) = aic(var_logl_lsar,n,len_logl);
        aic_lynx(n) = aic(var_lynx_lsar,n,len_lynx);
        aic_sun(n) = aic(var_sun_lsar,n,len_sun);
        aicc_logl(n) = aicc(var_logl_lsar,n,len_logl);
        aicc_lynx(n) = aicc(var_lynx_lsar,n,len_lynx);
        aicc_sun(n) = aicc(var_sun_lsar,n,len_sun);
        gic_logl(n) = gic(var_logl_lsar,n,len_logl,4);
        gic_lynx(n) = gic(var_lynx_lsar,n,len_lynx,4);
        gic_sun(n) = gic(var_sun_lsar,n,len_sun,4);
        bic_logl(n) = bic(var_logl_lsar,n,len_logl);
        bic_lynx(n) = bic(var_lynx_lsar,n,len_lynx);
        bic_sun(n) = bic(var_sun_lsar,n,len_sun);
        
end

figure(1)
subplot(2,3,1), plot(omega_logl,phi_logl_lsar), title('Loglynx LSAR'), ylabel('spectrum [dB]'),xlabel('frequency [rad]');
subplot(2,3,2), plot(omega_lynx,phi_lynx_lsar), title('Lynx LSAR'), ylabel('spectrum [dB]'),xlabel('frequency [rad]');
subplot(2,3,3), plot(omega_sun,phi_sun_lsar), title('Sunspots LSAR'), ylabel('spectrum [dB]'),xlabel('frequency [rad]');


sys_logl_lsar = tf(1,params_logl_lsar');
sys_lynx_lsar = tf(1,params_lynx_lsar');
sys_sun_lsar = tf(1,params_sun_lsar');

w = 0:1/100:2*pi;
unitcircle = cos(w)+1j*sin(w);


subplot(2,3,4), pzplot(sys_logl_lsar), hold on, plot(unitcircle), title('PZ-plot Loglynx');
subplot(2,3,5), pzplot(sys_lynx_lsar), hold on, plot(unitcircle), title('PZ-plot Lynx');
subplot(2,3,6), pzplot(sys_sun_lsar), hold on, plot(unitcircle), title('PZ-plot Loglynx');

figure(3)
subplot(1,3,1), plot(1:orders,vars_logl), grid on, title('variance vs order loglynx'), xlabel('order'),ylabel('white noise variance estimate');
subplot(1,3,2), plot(1:orders,vars_lynx), grid on, title('variance vs order lynx'),xlabel('order'),ylabel('white noise variance estimate');
subplot(1,3,3), plot(1:orders,vars_sun), grid on, title('variance vs order sunspot'),xlabel('order'),ylabel('white noise variance estimate');

figure(4)
subplot(2,2,1), plot(1:orders,aic_logl), title('AIC vs order loglynx');
subplot(2,2,2), plot(1:orders,aicc_logl), title('AICC vs order loglynx');
subplot(2,2,3), plot(1:orders,gic_logl), title('GIC vs order loglynx');
subplot(2,2,4), plot(1:orders,bic_logl), title('BIC vs order loglynx');

figure(5)
subplot(2,2,1), plot(1:orders,aic_lynx), title('AIC vs order lynx');
subplot(2,2,2), plot(1:orders,aicc_lynx), title('AICC vs order lynx');
subplot(2,2,3), plot(1:orders,gic_lynx), title('GIC vs order lynx');
subplot(2,2,4), plot(1:orders,bic_lynx), title('BIC vs order lynx');

figure(6)
subplot(2,2,1), plot(1:orders,aic_sun), title('AIC vs order sunspot');
subplot(2,2,2), plot(1:orders,aicc_sun), title('AICC vs order sunspot');
subplot(2,2,3), plot(1:orders,gic_sun), title('GIC vs order sunspot');
subplot(2,2,4), plot(1:orders,bic_sun), title('BIC vs order sunspot');



function a = lsmerge(coeff,omega,index)
    a = coeff*exp(-1j*omega*index);
end

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