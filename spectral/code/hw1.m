clear all
close all

% normal dist. noise. 50 realizations. 256 datapoints.
N=256;
e=randn(N,50);

% broadband arma-coeffs
A=[1 -1.3817 1.5632 -0.8843 0.4096];
B=[1 0.3544 0.3508 0.1736 0.2401];
e_filt=filter(B,A,e);

% narrowband arma-coeffs
A2=[1 -1.6408 2.2044 -1.4808 0.8145];
B2=[1 1.5857 0.9604];
e_filt_n=filter(B2,A2,e);

M_N4=N/8;
M_N16=N/16;

K_N4=M_N4/2;
K_N16=M_N16/2;

L_N4=N/M_N4;
L_N16=N/M_N16;

% bartlett-window
bw_N4=bartlett(M_N4);
bw_N16=bartlett(M_N16);
% rectangular window
rw_N4=rectwin(M_N4);
rw_N16=rectwin(M_N16);
% home made triangular window
win_N4=zeros(1,M_N4);
win_N16=zeros(1,M_N16);
for i=1:M_N4
    win_N4(i) = 1-(i-1)/(M_N4-1);
end
for i=1:M_N16
    win_N16(i)=1-(i-1)/(M_N16-1);
end

% broadband init
bt_N4=zeros(N,50);
bt_N16=zeros(N,50);
w_N4=zeros(N,50);
w_N16=zeros(N,50);
mean_w=zeros(2,256);
mean_bt=zeros(2,256);
std_w=zeros(2,256);
std_bt=zeros(2,256);
% narrowband init
bt_N4_n=zeros(N,50);
bt_N16_n=zeros(N,50);
w_N4_n=zeros(N,50);
w_N16_n=zeros(N,50);
mean_w_n=zeros(2,256);
mean_bt_n=zeros(2,256);
std_w_n=zeros(2,256);
std_bt_n=zeros(2,256);
for i=1:50
    % broadband
    bt_N4(:,i)=btse(e_filt(:,i),win_N4,N);
    bt_N16(:,i)=btse(e_filt(:,i),win_N16,N);
    w_N4(:,i)=welchse(e_filt(:,i),win_N4,K_N4,N);
    w_N16(:,i)=welchse(e_filt(:,i),win_N16,K_N16,N);
    % narrowband
    bt_N4_n(:,i)=btse(e_filt_n(:,i),win_N4,N);
    bt_N16_n(:,i)=btse(e_filt_n(:,i),win_N16,N);
    w_N4_n(:,i)=welchse(e_filt_n(:,i),win_N4,K_N4,N);
    w_N16_n(:,i)=welchse(e_filt_n(:,i),win_N16,K_N16,N);
end
% mean_bt(1,:)=mean(bt_N4,2);
% mean_bt(2,:)=mean(bt_N16,2);
% std_bt(1,:)=std(bt_N4,[],2);
% std_bt(2,:)=std(bt_N16,[],2);
% mean_bt_n(1,:)=mean(bt_N4_n,2);
% mean_bt_n(2,:)=mean(bt_N16_n,2);
% std_bt_n(1,:)=std(bt_N4_n,[],2);
% std_bt_n(2,:)=std(bt_N16_n,[],2);
% 
% mean_w(1,:)=mean(w_N4,2);
% mean_w(2,:)=mean(w_N16,2);
% std_w(1,:)=std(w_N4,0,2);
% std_w(2,:)=std(w_N16,0,2);
% mean_w_n(1,:)=mean(w_N4_n,2);
% mean_w_n(2,:)=mean(w_N16_n,2);
% std_w_n(1,:)=std(w_N4_n,0,2);
% std_w_n(2,:)=std(w_N16_n,0,2);

% converting to [-pi,pi]
mean_bt(1,1:N/2)=mean(bt_N4((N/2)+1:N,:),2);
mean_bt(1,(N/2)+1:N)=mean(bt_N4(1:N/2,:),2);
mean_bt(2,1:N/2)=mean(bt_N16((N/2)+1:N,:),2);
mean_bt(2,(N/2)+1:N)=mean(bt_N16(1:N/2,:),2);
std_bt(1,1:N/2)=std(bt_N4((N/2)+1:N,:),[],2);
std_bt(1,(N/2)+1:N)=std(bt_N4(1:N/2,:),[],2);
std_bt(2,1:N/2)=std(bt_N16((N/2)+1:N,:),[],2);
std_bt(2,(N/2)+1:N)=std(bt_N16(1:N/2,:),[],2);
mean_bt_n(1,1:N/2)=mean(bt_N4_n((N/2)+1:N,:),2);
mean_bt_n(1,(N/2)+1:N)=mean(bt_N4_n(1:N/2,:),2);
mean_bt_n(2,1:N/2)=mean(bt_N16_n((N/2)+1:N,:),2);
mean_bt_n(2,(N/2)+1:N)=mean(bt_N16_n(1:N/2,:),2);
std_bt_n(1,1:N/2)=std(bt_N4_n((N/2)+1:N,:),[],2);
std_bt_n(1,(N/2)+1:N)=std(bt_N4_n(1:N/2,:),[],2);
std_bt_n(2,1:N/2)=std(bt_N16_n((N/2)+1:N,:),[],2);
std_bt_n(2,(N/2)+1:N)=std(bt_N16_n(1:N/2,:),[],2);

mean_w(1,1:N/2)=mean(w_N4((N/2)+1:N,:),2);
mean_w(1,(N/2)+1:N)=mean(w_N4(1:N/2,:),2);
mean_w(2,1:N/2)=mean(w_N16((N/2)+1:N,:),2);
mean_w(2,(N/2)+1:N)=mean(w_N16(1:N/2,:),2);
std_w(1,1:N/2)=std(w_N4((N/2)+1:N,:),[],2);
std_w(1,(N/2)+1:N)=std(w_N4(1:N/2,:),[],2);
std_w(2,1:N/2)=std(w_N16((N/2)+1:N,:),[],2);
std_w(2,(N/2)+1:N)=std(w_N16(1:N/2,:),[],2);
mean_w_n(1,1:N/2)=mean(w_N4_n((N/2)+1:N,:),2);
mean_w_n(1,(N/2)+1:N)=mean(w_N4_n(1:N/2,:),2);
mean_w_n(2,1:N/2)=mean(w_N16_n((N/2)+1:N,:),2);
mean_w_n(2,(N/2)+1:N)=mean(w_N16_n(1:N/2,:),2);
std_w_n(1,1:N/2)=std(w_N4_n((N/2)+1:N,:),[],2);
std_w_n(1,(N/2)+1:N)=std(w_N4_n(1:N/2,:),[],2);
std_w_n(2,1:N/2)=std(w_N16_n((N/2)+1:N,:),[],2);
std_w_n(2,(N/2)+1:N)=std(w_N16_n(1:N/2,:),[],2);

p=[-pi:2*pi/255:pi];
%dB = mag2db(abs(h));

% plots

figure(1)
subplot(2,2,1), plot(p,mag2db(abs(mean_bt(1,:)))),hold on,...
    plot(p,mag2db(abs(mean_bt(1,:)+std_bt(1,:)))),...
    plot(p,mag2db(abs(mean_bt(1,:)-std_bt(1,:)))),
grid on, title('Broadband BT sample mean for M=N/8'),
xlabel('Frequency [rad]'),
ylabel('Amplitude [dB]');
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

subplot(2,2,2), plot(p,mag2db(abs(mean_bt(2,:)))),hold on,...
    plot(p,mag2db(abs(mean_bt(2,:)+std_bt(2,:)))),...
    plot(p,mag2db(abs(mean_bt(2,:)-std_bt(2,:)))),
grid on, title('Broadband BT sample mean for M=N/16'),
xlabel('Frequency [rad]'),
ylabel('Amplitude [dB]');
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

subplot(2,2,3), plot(p,mag2db(abs(mean_w(1,:)))),hold on,...
    plot(p,mag2db(abs(mean_w(1,:)+std_w(1,:)))), ...
    plot(p,mag2db(abs(mean_w(1,:)-std_w(1,:)))),
grid on, title('Broadband Welch sample mean for M=N/8'),
xlabel('Frequency [rad]'),
ylabel('Amplitude [dB]');
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

subplot(2,2,4), plot(p,mag2db(abs(mean_w(2,:)))),hold on,...
    plot(p,mag2db(abs(mean_w(2,:)+std_w(2,:)))),...
    plot(p,mag2db(abs(mean_w(2,:)-std_w(2,:)))),
grid on, title('Broadband Welch sample mean for M=N/16'),
xlabel('Frequency [rad]'),
ylabel('Amplitude [dB]');
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

figure(2)

subplot(2,2,1), plot(p,mag2db(abs(mean_bt_n(1,:)))),hold on,...
    plot(p,mag2db(abs(mean_bt_n(1,:)+std_bt_n(1,:)))),...
    plot(p,mag2db(abs(mean_bt_n(1,:)-std_bt_n(1,:)))),
grid on, title('Narrowband BT sample mean for M=N/8'),
xlabel('Frequency [rad]'),
ylabel('Amplitude [dB]');
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

subplot(2,2,2), plot(p,mag2db(abs(mean_bt_n(2,:)))),hold on,...
    plot(p,mag2db(abs(mean_bt_n(2,:)+std_bt_n(2,:)))),...
    plot(p,mag2db(abs(mean_bt_n(2,:)-std_bt_n(2,:)))),
grid on, title('Narrowband BT sample mean for M=N/16'),
xlabel('Frequency [rad]'),
ylabel('Amplitude [dB]');
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

subplot(2,2,3), plot(p,mag2db(abs(mean_w_n(1,:)))),hold on,...
    plot(p,mag2db(abs(mean_w_n(1,:)+std_w_n(1,:)))),...
    plot(p,mag2db(abs(mean_w_n(1,:)-std_w_n(1,:)))),
grid on, title('Narrowband Welch sample mean for M=N/8'),
xlabel('Frequency [rad]'),
ylabel('Amplitude [dB]');
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

subplot(2,2,4), plot(p,mag2db(abs(mean_w_n(2,:)))),hold on,...
    plot(p,mag2db(abs(mean_w_n(2,:)+std_w_n(2,:)))),...
    plot(p,mag2db(abs(mean_w_n(2,:)-std_w_n(2,:)))),
grid on, title('Narrowband Welch sample mean for M=N/16'),
xlabel('Frequency [rad]'),
ylabel('Amplitude [dB]');
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

figure(3)

subplot(2,1,1),title('Broadband variances'),hold on, ... 
    plot(p,std_bt(1,:).^2), plot(p,std_bt(2,:).^2), plot(p,std_w(1,:).^2), plot(p,std_w(2,:).^2),
xlabel('Frequency [rad]'),
ylabel('Variance');
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
legend('BT N/8','BT N/16','Welch N/8','Welch N/16')

subplot(2,1,2),title('Narrowband variances'),hold on, ...
    plot(p,std_bt_n(1,:).^2), plot(p,std_bt_n(2,:).^2), plot(p,std_w_n(1,:).^2), plot(p,std_w_n(2,:).^2),
xlabel('Frequency [rad]'),
ylabel('Variance');
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
legend('BT N/8','BT N/16','Welch N/8','Welch N/16')

% true spectrum
figure(4)
plot(log(periodogram(e_filt)))
hold on
plot(log(periodogram(e_filt_n)))
