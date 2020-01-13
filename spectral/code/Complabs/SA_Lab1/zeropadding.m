function zeropadding
global N
N=64;        % Number of samples
phi_1=0;     % phase offset of first sinusoid
phi_2=0;     % phase offset of second sinusoid
A_1=10;      % amplitude of first sinusoid
A_2=5;       % amplitude of second sinusoid
f_1=0.2;     % frequency of first sinusoid
f_2=f_1+1/N; % frequency of second sinusoid
sig2=1;      % noise variance

% Time variable
t=(0:N-1).';

% Generate noise sequence
e=sqrt(sig2)*randn(N,1); % white Gaussion noise with variance sig2 

% Observed Signal: two sinusoids plus white gaussian noise
y=A_1*sin(2*pi*f_1*t+phi_1)+A_2*sin(2*pi*f_2*t+phi_2)+e;

% window
v=ones(N,1);

% GUI code
%*********
global stop
global slide
global quit_button

stop=0;
f_gui=figure('Position',[680 300 800 600]);
plot_axis=axes('position',[0.1 0.4 0.35 0.5],'NextPlot','add');
axis([-pi pi 0 1700])
xlabel('\omega [rad]')                    % axis labels
ylabel('\phi_p')                          
title('Periodogram') 

plot_axis_z=axes('position',[0.55 0.4 0.35 0.5],'NextPlot','add');
axis([2*pi*f_1-pi/10 2*pi*f_2+pi/10 0 1700])
xlabel('\omega [rad]')                          % axis labels
ylabel('\phi_p')                                % 
title('Zoomed Periodogram (zoomable)')          % title

plot_data=axes('position',[0.1 0.16 0.80 0.1],'NextPlot','add');
axis([1 1024 -12 12])                 % set axis range
title('Data')                                   % title 
xlabel('Samples')                               % axis labels
ylabel('Amplitude')  
  
slide=uicontrol(f_gui,'Style','slider',...
		'Min',N,'Max',N*16,'Value',N,'String','banan','Units', ...
		'normalized','Position',[0.1 0.02 0.8 0.05]);

quit_button=uicontrol(f_gui,'style','pushbutton','String','Quit','Callback', ...
	       @quit_call,'Units','normalized','position',[0.90 0.90 0.1 0.1]);

text(1040,10,'Total Length')                    % heading
       
[phi_p_f,w_f]=periodogram_se(y,v,8192*16);
while stop~=1
  %Remove previous plots
  cla(plot_axis)
  cla(plot_axis_z)
  cla(plot_data)
  
  L=round(get(slide,'Value')/2)*2;
  set(slide,'Value',L); %make sure length is even
  [phi_p,w]=periodogram_se(y,v,L); % periodogram
                                     % erase earlier 
  plot(plot_axis,w_f,phi_p_f,'k')                     % plot "theoretic" periodogram                                   % allow superposition plots
  plot(plot_axis,w,phi_p)                             % plot periodogram
  plot(plot_axis,2*pi*[f_1 f_1],[1 1700],'r:')        % mark frequencies
  plot(plot_axis,-2*pi*[f_1 f_1],[1 1700],'r:')       %
  plot(plot_axis,2*pi*[f_2 f_2],[1 1700],'r:')        %  
  plot(plot_axis,-2*pi*[f_2 f_2],[1 1700],'r:')       %
   
  % zoomed periodogram plot
  plot(plot_axis_z,w_f,phi_p_f,'k')                   % plot "theoretical"                                       % allow superposition
  plot(plot_axis_z,w,phi_p)                           % periodogram
  plot(plot_axis_z,2*pi*[f_1 f_1],[1 1700],'r:')      % mark frequencies
  plot(plot_axis_z,-2*pi*[f_1 f_1],[1 1700],'r:')     %
  plot(plot_axis_z,2*pi*[f_2 f_2],[1 1700],'r:')      %
  plot(plot_axis_z,-2*pi*[f_2 f_2],[1 1700],'r:')     %
  legend(plot_axis_z,'Theor.','Zerop.','True');       % only needed once really...

  zoom on    
  
  % plot data used                                  
  plot(plot_data,y(:).','r')                          % Plot nonzero data                                        % allow superposition
  plot(plot_data,zeros(1,L),'k')                      % plot samples used
  plot(plot_data,[1 1],[-4 4],'k')                    %
  plot(plot_data,[1 1]*L,[-4 4],'k')                  %

                           
  
  %display data length used (data+zeropad)
  title(plot_data,['Data: ' num2str(L) ' samples'])  
  
  %wait for new Value given using the slider
  waitfor(slide,'Value')
end
close(f_gui)

%******************************
% Periodogram function
%*******************************
function [phi_p,w]=periodogram_se(y,v,L)
% The windowed periodogram estimate
%
%  [phi_p,w]=periodogram_se(y,v,L)
%
% y - data sequence of length N
% v - Apodization (time window) of data (symetric tapering of length N)
% L - data length + zero padding (N + number of zeros)
%
% phi_p - the periodogram spectral estimate at w_k=2*pi*k/L-pi where k=0:L-1
% w - frequency points wk=2*pi*k/L-pi where k=0:L-1

y=y(:); % make sure data is in a column vector
v=v(:); % make sure window is in a column vector

N=length(y); % data length

if N~=length(v)
  disp('*** Error:  Window length does not match data size ***')
  return
end

phi_p=1/N*abs(fft(y.*v,L)).^2; % periodogram estimate

phi_p=fftshift(phi_p); % give spectrum -pi:pi instead of 0:2pi 

w=((0:L-1)/L-0.5)*2*pi;
return


%**************
% exit gui function
%**************
function quit_call(~,~)
global slide
global stop
global N
stop=1;                          % lit stop flag
set(slide,'Value',mod(get(slide,'Value')-N+1,15*N+1)+N) %Map to new value in given interval N to 16*N
return

