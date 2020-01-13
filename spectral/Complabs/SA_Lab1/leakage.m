function leakage
%Updated 2016, MB

%GUI
%%%%%%%%%%
global stop
global slide_freq
global slide_amp
global window_button
global quit_button
global h1
global h2
global h3
global h4
global y
global N
global f_1
global trigger_obj
global f_gui

N=256;
phi_1=0;     % phase offset of first sinusoid
phi_2=0;     % phase offset of second sinusoid
A_1=1;       % amplitude of first sinusoid
A_2=1;       % amplitude of second sinusoid
f_1=0.2;     % frequency of first sinusoid
sig2=0;      % noise variance
alpha=4;
% Time variable
t=(0:N-1).';

% Generate noise sequence
e=sqrt(sig2)*randn(N,1); % white Gaussion noise with variance sig2 

% Compute Chebychew window
window_func=chebwin(N,60); %Use chebwin with 6dB sidelobe level


stop=0;

f_gui=figure('Position',[500 300 820 550]);

plot_axis_lin=axes('position',[0.10 0.55 0.4 0.4]);
plot_axis_dB=axes('position',[0.58 0.55 0.4 0.4]);


slide_pos=[0.09 0.17 0.4 0.05];
slide_freq=uicontrol(f_gui,'Style','slider',...
		'Min',0,'Max',12,'Value',alpha,'Units', ...
		'normalized','Position',slide_pos,'Callback',@newAlpha);

slide_amp_pos=[0.58 0.17 0.4 0.05];
slide_amp=uicontrol(f_gui,'Style','slider',...
		'Min',-3,'Max',0,'Value',0,'Units', ...
		'normalized','Position',slide_amp_pos,'Callback',@newAmp);
quit_button_pos=[0.90 0 0.10 0.10];
quit_button=uicontrol(f_gui,'style','pushbutton','String','Quit','Callback', ...
	       @quit_call,'Units','normalized','position',quit_button_pos);
window_button_pos=[0.09 0.08 0.18 0.07];
window_button=uicontrol(f_gui,'style','checkbox','value',0,'Callback',...
		   @window_check,'String','60dB Cheb Win','Units','normalized','position',...
			 window_button_pos);
cont_pos=[0.70 0.0 0.10 0.10];
trigger_obj=uicontrol(f_gui,'style','pushbutton','String','Compute spectrum','Units','normalized','position',cont_pos,'Visible','off');
plotButton=uicontrol(f_gui,'style','pushbutton','String','Plot signal and ACS','Units','normalized',...
   'position',[0.75 0.0 0.15 0.10],'Callback',@plotSignal);
% sampleBox=uicontrol(f_gui,'style','edit','String',256,'Units','normalized','position',...
% 			 [0.40 0.08 0.10 0.07],'callback',@sampleCall,'Fontsize',10);
         
  text_axis=axes('position',[0.09 0.25 0.9 0.2],'visible','off');
  
  text(0,1,'Settings:')
%   text(0.30,-0.65,'N =')
  
  text(0,0.7,['2\pi f_0=' num2str(2*pi*f_1)])
  h1=text(0.2,0,['\alpha=' num2str(alpha)]);
  h2=text(0,0.4,['2\pi(f_0+\alpha/N)=' num2str(2*pi*(f_1+alpha/N))]);
  text(0.55,0.7,['a_1=' num2str(A_1)])
  h3=text(0.74,0,['a_2=' num2str(A_2)]);
  h4=text(0.55,0.4,['a_2=' num2str(A_2)]);
  text(0,0,'0')
  text(0.42,0,'12')
  text(0.55,0,'0.001')
  text(0.96,0,'1')
  

 
  
while stop~=1
  alpha=get(slide_freq,'value');         % read frequency from slider
  A_2=10^(get(slide_amp,'value'));       % read amplitude from slider

  y=A_1*sin(2*pi*f_1*t+phi_1)+A_2*sin(2*pi*(f_1+alpha/N)*t+phi_2)+e;

  if get(window_button,'value')
     [phi,w]=periodogram_se(y,window_func,N*8);          % periodogram
  else
     [phi,w]=periodogram_se(y,ones(N,1),N*8);  % periodogram RECT window!
  end
 
  phi=phi/max(phi);                    % normalize
  if alpha==0 && abs(phi_2-pi)<5*eps
      phi=zeros(N*8,1); %To avoid numerical issues when the sinusoids cancel
  end
  
  axes(plot_axis_lin)
  hold off
  plot(w,phi)                          % plot linear scale
  hold on
  plot([1 1]*2*pi*f_1,[0 1],'r')         % half power indicator
  plot([1 1]*2*pi*(f_1+alpha/N),[0 A_2^2],'r')         % half power indicator
  zoom on
  axis([2*pi*(f_1-5/N) 2*pi*(f_1+20/N)  0 1])
  title('Periodogram estimate, linear scale')
  xlabel('\omega')
  ylabel('Normalized \phi_p')
  axes(plot_axis_dB)
  hold off
  plot(w,10*log10(phi+eps))            % plot dB scale
  hold on
  plot([1 1]*2*pi*f_1,[-65 0],'r')         
  plot([1 1]*2*pi*(f_1+alpha/N),[-65 20*log10(A_2)],'r') 
  zoom on
  axis([2*pi*(f_1-5/N) 2*pi*(f_1+20/N)  -65 0])
  title('Periodogram estimate, dB scale')
  xlabel('\omega')
  ylabel('\phi_p [dB]')
  
  set(trigger_obj,'enable','off')
  waitfor(trigger_obj,'enable','on')     
  
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
% v - Apodization (time window) of data (symmetric tapering of length N)
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

function newAmp(~,~) %Two inputs are required for callback functions
    global slide_amp
    global h3
    global h4
    global trigger_obj
    set(h3,'String',['a_2=' num2str(10^get(slide_amp,'value'))])
    set(h4,'String',['a_2=' num2str(10^get(slide_amp,'value'))])
    set(trigger_obj,'enable','on')
return

function newAlpha(~,~)
    global slide_freq
    global h1
    global h2
    global f_1
    global N
    global trigger_obj
    set(h1,'String',['\alpha=' num2str(get(slide_freq,'value'))])
    set(h2,'String',['2\pi(f_0+\alpha/N)=' num2str(2*pi*(f_1+get(slide_freq,'value')/N))])
    set(trigger_obj,'enable','on')
return

function plotSignal(~,~)
    global y
    N=length(y);   
    figure('Position',[600 500 820 360]);
    subplot(1,2,1)
    plot(0:N-1,y)
    xlim([0 N-1])
    title('Signal')
    xlabel('k')
    ylabel('y(k)')
    subplot(1,2,2)
    plot(-N+1:N-1,xcorr(y,'biased'))
    xlim([-N+1 N-1])
    title('ACS estimate (biased)')
    xlabel('k')
    ylabel('$\hat{r}(k)$','interpreter','latex')
return

% function cont_pushed(~,~)
% global cont
% set(cont,'visible','off')
% return

function window_check(~,~)
global trigger_obj
set(trigger_obj,'enable','on')
return

% function sampleCall(src,~)
% str=get(src,'String');
% if isempty(str2num(str))
%     set(src,'string','0');
%     warndlg('Input must be numerical');
% end
% return

%**************
% exit gui function
%**************
function quit_call(~,~)
global stop
global trigger_obj
stop=1;                          % lit stop flag
set(trigger_obj,'enable','on')
return

