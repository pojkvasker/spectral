function resolution
%Updated 2017, MB

%GUI
%%%%%%%%%%
global stop
global slide_freq
global slide_phase
global hamming_button
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
alpha=3;

s=rng;
% Compute Hamming window
ham_win=hamming(N);


stop=0;

f_gui=figure('Position',[500 300 820 550]);

plot_axis_lin=axes('position',[0.10 0.55 0.4 0.4]);
plot_axis_dB=axes('position',[0.58 0.55 0.4 0.4]);


slide_pos=[0.09 0.17 0.4 0.05];
slide_freq=uicontrol(f_gui,'Style','slider',...
		'Min',0,'Max',3,'Value',alpha,'Units', ...
		'normalized','Position',slide_pos,'Callback',@newAlpha);

slide_phase_pos=[0.58 0.17 0.4 0.05];
slide_phase=uicontrol(f_gui,'Style','slider','Min',0,'Max',2*pi,'Value',0,'Units', ...
		'normalized','Position',slide_phase_pos,'Callback',@newPhi);
quit_button_pos=[0.90 0 0.10 0.10];
quit_button=uicontrol(f_gui,'style','pushbutton','String','Quit','Callback', ...
	       @quit_call,'Units','normalized','position',quit_button_pos);
hamming_button_pos=[0.09 0.08 0.18 0.07];
hamming_button=uicontrol(f_gui,'style','checkbox','value',0,'Callback',...
		   @hamming_check,'String','Hamming window','Units','normalized','position',...
			 hamming_button_pos);
cont_pos=[0.70 0.0 0.10 0.10];
trigger_obj=uicontrol(f_gui,'style','pushbutton','String','Compute spectrum','Units','normalized','position',cont_pos,'visible','off');
plotButton=uicontrol(f_gui,'style','pushbutton','String','Plot signal and ACS','Units','normalized',...
   'position',[0.75 0.0 0.15 0.10],'Callback',@plotSignal);
%sampleBox=uicontrol(f_gui,'style','edit','String',256,'Units','normalized','position',...
%			 [0.40 0.08 0.10 0.07],'callback',@sampleCall,'Fontsize',10);

     
Ncontrol = uicontrol('style','edit','units','normalized','position',[0.85 .35 0.10 0.05],...
            'string',num2str(N),'callback',@N_call);

  text_axes=axes('position',[0.09 0.25 0.9 0.2],'visible','off');
  text(0,1,'Settings:')
  %text(0.30,-0.65,'N =')
  text(0.80,0.62,['N ='])
  text(0,0.7,['2\pi f_0=' num2str(2*pi*f_1)])
  h1=text(0.2,0,['\alpha=' num2str(alpha)]);
  h2=text(0,0.4,['2\pi(f_0+\alpha/N)=' num2str(2*pi*(f_1+alpha/N))]);
  text(0.55,0.7,['\phi_1=' num2str(phi_1)])
  h3=text(0.74,0,['\phi_2=' num2str(phi_2)]);
  h4=text(0.55,0.4,['\phi_2=' num2str(phi_2)]);
  text(0,0,'0')
  text(0.42,0,'3')
  text(0.55,0,'0')
  text(0.96,0,'2\pi')

while stop~=1
  alpha=get(slide_freq,'value');         % read frequency from slider
  phi_2=get(slide_phase,'value');        % read phase from slider
  %axes(text_axis)
  %cla

  % Generate noise sequence
  rng(s); %Set to current seed for each noise generation
  e=sqrt(sig2)*randn(N,1); % white Gaussion noise with variance sig2 
  y=A_1*sin(2*pi*f_1*(0:N-1).'+phi_1)+A_2*sin(2*pi*(f_1+alpha/N)*(0:N-1).'+phi_2)+e;

  if get(hamming_button,'value')
     [phi,w]=periodogram_se(y,ham_win,N*8);          % periodogram
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
  plot([1 1]*2*pi*(f_1+alpha/N),[0 1],'r')         % half power indicator
  axis([2*pi*(f_1-4/N) 2*pi*(f_1+4/N)  0 1])
  title('Periodogram estimate, linear scale')
  xlabel('\omega')
  ylabel('Normalized \phi_p')
  axes(plot_axis_dB)
  hold off
  plot(w,10*log10(phi+eps))            % plot dB scale
  hold on
  plot([1 1]*2*pi*f_1,[-50 0],'r')         
  plot([1 1]*2*pi*(f_1+alpha/N),[-50 0],'r') 
  axis([2*pi*(f_1-4/N) 2*pi*(f_1+4/N)  -50 0])
  title('Periodogram estimate, dB scale')
  xlabel('\omega')
  ylabel('\phi_p [dB]')
  
  set(trigger_obj,'enable','off')
  waitfor(trigger_obj,'enable','on')     
  
end
close(f_gui)
end

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
end

function newPhi(~,~) %Two inputs are required for callback functions
    global slide_phase
    global h3
    global h4
    global trigger_obj
    set(h3,'String',['\phi_2=' num2str(get(slide_phase,'value'))])
    set(h4,'String',['\phi_2=' num2str(get(slide_phase,'value'))])
    set(trigger_obj,'enable','on')
end

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
end

function plotSignal(~,~)
    global y
    global N  
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
end

function hamming_check(~,~)
global trigger_obj
set(trigger_obj,'enable','on')
end

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
end

function N_call(src,~)
    global N
    global trigger_obj
    [num, status] = str2num(get(src,'String'));
    if status==1
        N=num;
    else
        set(src,'string','256');
        N=256;
        waitfor(warndlg('Input must be numerical. Setting N = 256.'));
    end
    set(trigger_obj,'enable','on')
end
