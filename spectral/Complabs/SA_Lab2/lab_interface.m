function lab_interface
global Y
global Y_n
global Y_b
global psd_true
global psd_b_true
global psd_n_true
global process_menu
global method_menu
global MA_menu
global AR_menu
global MA_text
global AR_text
global K_menu
global M_menu
global next_button
global quit_button
global stop
global poles_n_true
global zeros_n_true
global poles_true
global zeros_true
global zeros_b_true
global poles_b_true


sig_2=1;                                % driving noise variance
A_b=[1 -1.3817 1.5632 -0.8843 0.4096];  % AR coefficients broadband
B_b=[1 0.3544 0.3508 0.1736 0.2401];    % MA coefficients broadband
A_n=[1 -1.6408 2.2044 -1.4808 0.8145];  % AR coefficients narrow
B_n=[1 1.5857 0.9604];    % MA coefficients narrow

K_A_b=length(A_b)-1;                      % AR order broadband
K_B_b=length(B_b)-1;                      % MA order broadband
K_A_n=length(A_n)-1;                      % AR order narrowband
K_B_n=length(B_n)-1;                      % MA order narrowband
R=10;                                   % number of realizations
N=256;                                  % samples in each realization
E=randn(R,N);                           % driving noise
Y_b=filter(B_b,A_b,E,[],2);               % generate broadband data
Y_n=filter(B_n,A_n,E,[],2);               % generate narrowband data
Y=Y_b;

poles_n_true=roots(A_n);
zeros_n_true=roots(B_n);
poles_b_true=roots(A_b);
zeros_b_true=roots(B_b);

poles_true=poles_b_true;
zeros_true=zeros_b_true;

w=(0:(4*N-1))/(4*N)*2*pi-pi;                  % frequency axis
A_n_w=A_n*exp(1j*(0:K_A_n).'*w);           % A_n(omega)
B_n_w=B_n*exp(1j*(0:K_B_n).'*w);           % B_n(omega)
psd_n_true=B_n_w.*conj(B_n_w)./(A_n_w.*conj(A_n_w))*sig_2; % true PSD narrowband

A_b_w=A_b*exp(1j*(0:K_A_b).'*w);           % A_b(omega)
B_b_w=B_b*exp(1j*(0:K_B_b).'*w);           % B_b(omega)
psd_b_true=B_b_w.*conj(B_b_w)./(A_b_w.*conj(A_b_w))*sig_2; % true PSD narrowband
psd_true=psd_b_true;

gui_fig=figure('Position',[500 200 1200 800]);
plot_axis=axes('position',[0.05 0.55 0.4 0.4]);
lin_axis=axes('position',[0.05 0.07 0.4 0.4]);
pz_axis=axes('position',[0.55 0.07 0.4 0.4]);

text_axis=axes('position',[0.55 1 0.4 0.4],'Visible','off');
text(0,-0.08,'Process')
text(0,-0.33,'Method')
AR_text=text(0,-0.7,'AR order');
MA_text=text(0,-0.95,'MA order','Visible','off');

process_pos=[0.55 0.90 0.15 0.05];
method_pos=[0.55 0.80 0.15 0.05];
M_pos=[0.6 0.75 0.1 0.05];
MA_pos=[0.55 0.55 0.15 0.05];

AR_pos=[0.55 0.65 0.15 0.05];
K_pos=[0.6 0.75 0.1 0.05];
next_pos=[0.825 0.55 0.15 0.15];
quit_pos=[0.9,0.9,0.10,0.10];
process_menu=uicontrol(gui_fig,'style','popup','string',...
'Broadband|Narrowband','value',1,'callback',@process_call,...
'units','normalized','position',process_pos);

method_menu=uicontrol(gui_fig,'style','popup',...
 'string','YW|LSAR|MYW|LSARMA|periodogram','value',1,...
'callback',@method_call,'units','normalized','position', method_pos);

M_menu=uicontrol(gui_fig,'style','popup','string',...
'M=n|M=2n','value',1,'callback',@trigger_call,...
'units','normalized','position', M_pos,'visible','off');

MA_menu=uicontrol(gui_fig,'style','popup','string',...
'm=0|m=2|m=4|m=6|m=8','value',1,'callback',@trigger_call,...
'units','normalized','position', MA_pos,'visible','off');

AR_menu=uicontrol(gui_fig,'style','popup','string',...
'n=4|n=8|n=12|n=16','value',1,'callback',@trigger_call,...
'units','normalized','position',AR_pos);

K_menu=uicontrol(gui_fig,'style','popup','string',...
'K=n|K=2n|K=3n|K=4n','value',2,'callback',@trigger_call,...
'units','normalized','position',K_pos,'visible','off');

next_button=uicontrol(gui_fig,'style','pushbutton','string','Compute estimates',...
    'units','normalized','position',next_pos,'Visible','off');

quit_button=uicontrol(gui_fig,'style','pushbutton','string','Quit',...
'callback',@quit_call,'units','normalized','position',quit_pos);

stop=0;
n=[4 8 12 16];
m=[0 2 4 6 8];
while stop~=1
  axes(plot_axis)
  switch get(method_menu,'value')
   
   case 1
    %yule-walker
    n_=n(get(AR_menu,'value'));
    a=complex(zeros(n_+1,R));
    a_w=complex(zeros(length(w),R));
    sig2=zeros(1,R);
    poles_=complex(zeros(n_,R));
    zeros_=[];
    for k=1:R
      [a(:,k),sig2(:,k)]=yulewalker(Y(k,:),n_);
      a_w(:,k)=(a(:,k).'*exp(1j*(0:n_).'*w)).';
      poles_(:,k)=roots(a(:,k));
    end
    psd=ones(length(w),1)*sig2./(a_w.*conj(a_w));
    
   case 2
    %lsar
    n_=n(get(AR_menu,'value'));
    a=complex(zeros(n_+1,R));
    a_w=complex(zeros(length(w),R));
    sig2=zeros(1,R);
    poles_=complex(zeros(n_,R));
    zeros_=[];
    for k=1:R
      [a(:,k),sig2(:,k)]=lsar(Y(k,:),n_);
      a_w(:,k)=(a(:,k).'*exp(1j*(0:n_).'*w)).';
      poles_(:,k)=roots(a(:,k));
    end
    psd=ones(length(w),1)*sig2./(a_w.*conj(a_w));
   
   case 3
    % modified yule walker
    n_=n(get(AR_menu,'value'));
    m_=m(get(MA_menu,'value'));
    M_=get(M_menu,'value')*n_;
    a=complex(zeros(n_+1,R));
    a_w=complex(zeros(length(w),R));
    gamma=complex(zeros(m_+1,R));
    gamma_w=zeros(length(w),R);
    poles_=complex(zeros(n_,R));
    zeros_=complex(zeros(m_,R));
    for k=1:R
      [a(:,k),gamma(:,k)]=mywarma(Y(k,:),n_,m_,M_);
      a_w(:,k)=(a(:,k).'*exp(1j*(0:n_).'*w)).';
      poles_(:,k)=roots(a(:,k));
      if m_>1
	gamma_w(:,k)=gamma(1,k)+2*real(gamma(2:end,k).'*...
        exp(1j*(1:m_).'*w)).';
	z=roots(gamma(:,k));
	zeros_(:,k)=z(abs(z)<=1);
      else 
	gamma_w(:,k)=ones(length(w),1)*gamma(1,k);
	zeros_=[];
      end
    end
    psd=gamma_w./(a_w.*conj(a_w));
   case 4
   %lsarma
    n_=n(get(AR_menu,'value'));
    m_=m(get(MA_menu,'value')+1);
    M_=get(M_menu,'value')*n_;
    K=[n_ 2*n_ 3*n_ 4*n_];
    K_=K(get(K_menu,'value'));
    a=complex(zeros(n_+1,R));
    a_w=complex(zeros(length(w),R));
    b=complex(zeros(m_+1,R));
    b_w=complex(zeros(length(w),R));
    sig2=zeros(1,R);
    poles_=complex(zeros(n_,R));
    zeros_=complex(zeros(m_,R));
    for k=1:R
      [a(:,k),b(:,k),sig2(:,k)]=lsarma(Y(k,:),n_,m_,K_);
      a_w(:,k)=(a(:,k).'*exp(1j*(0:n_).'*w)).';
      b_w(:,k)=(b(:,k).'*exp(1j*(0:m_).'*w)).';
      zeros_(:,k)=roots(b(:,k));
      poles_(:,k)=roots(a(:,k));
    end
    psd=(b_w.*conj(b_w))./(a_w.*conj(a_w)).*(ones(length(w),1)*sig2);
   case 5
    %periodogram
    psd=zeros(length(w),R);
    for k=1:R
      psd(:,k)=fftshift(periodogramse(Y(k,:),ones(1,N),length(w)));
    end
  end
  figure(gui_fig)  
  p1=plot(w,real(10*log10(psd./max(psd_true))),'r--');
  hold on 
  p2=plot(w,10*log10(psd_true/max(psd_true)),'b');
  legend([p1(1) p2],'Estimated PSDs','True PSD')
  db_min=min(10*log10(psd_true/max(psd_true)+eps));
  axis([-pi pi db_min*1.1 -db_min*0.10])
  xlabel('\omega')
  ylabel('\phi_y [dB]')
  %zoom on
  hold off
 
  axes(lin_axis)
  pl1=plot(w,psd,'r--');
  hold on 
  pl2=plot(w,psd_true,'b');
  legend([pl1(1) pl2],'Estimated PSDs','True PSD')
  lin_max=max([max(psd_true) max(max(psd))]);
  axis([-pi pi -0.05*lin_max lin_max*1.05])
  xlabel('\omega')
  ylabel('\phi_y (linear scale)')
  %zoom on
  hold off
  
  if get(method_menu,'value')==5
      cla(pz_axis,'reset');
      set(pz_axis,'visible','off')
  else
      axes(pz_axis)
      plot(cos(2*pi*(0:512)/512),sin(2*pi*(0:512)/512),'k')
      hold on
      pp_1=plot(real(poles_), imag(poles_),'xr');
      pz_1=plot(real(zeros_), imag(zeros_),'or');
      pp_2=plot(real(poles_true), imag(poles_true),'xb');
      pz_2=plot(real(zeros_true), imag(zeros_true),'ob');
      if isempty(pz_1)
        legend([pp_2 pz_2 pp_1(1)],'True poles','True zeros',...
           'Estimated poles')
      else
        legend([pp_2 pz_2 pp_1(1) pz_1(1)],'True poles','True zeros',...
           'Estimated poles', 'estimated zeros')
      end
      xlabel('Re\{Z\}')
      ylabel('Im\{Z\} [j]')
      hold off
      axis square
      axis equal
      %zoom on
  end
  
  set(next_button,'enable','off')
  waitfor(next_button,'enable','on')
end
close(gui_fig)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function process_call(~,~)
global process_menu
global Y
global Y_n
global Y_b
global psd_true
global psd_b_true
global psd_n_true
global poles_n_true
global zeros_n_true
global poles_true
global zeros_true
global zeros_b_true
global poles_b_true
global next_button
switch get(process_menu,'value')
 case 1
  Y=Y_b;
  psd_true=psd_b_true;
  poles_true=poles_b_true;
  zeros_true=zeros_b_true;
 case 2
  Y=Y_n;
  psd_true=psd_n_true;
  poles_true=poles_n_true;
  zeros_true=zeros_n_true;
end
set(next_button,'enable','on')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function method_call(~,~)
global next_button
global M_menu
global K_menu
global MA_menu
global AR_menu
global MA_text
global AR_text
global method_menu
switch get(method_menu,'value')
 case 1
  set(K_menu,'visible','off')
  set(M_menu,'visible','off')
  set(MA_menu,'visible','off')
  set(MA_text,'visible','off')
  set(AR_menu,'visible','on')
  set(AR_text,'visible','on')
 case 2
  set(K_menu,'visible','off')
  set(M_menu,'visible','off')
  set(MA_menu,'visible','off')
  set(MA_text,'visible','off')
  set(AR_menu,'visible','on')
  set(AR_text,'visible','on')
 case 3
  set(M_menu,'visible','on')
  set(K_menu,'visible','off')
  set(MA_menu,'string','m=0|m=2|m=4|m=6|m=8')
  set(MA_menu,'visible','on')
  set(MA_text,'visible','on')
  set(AR_menu,'visible','on')
  set(AR_text,'visible','on')
 case 4
  set(K_menu,'visible','on')
  set(M_menu,'visible','off')
  set(MA_menu,'visible','on')
  set(MA_menu,'value',min(4,get(MA_menu,'value')))
  set(MA_menu,'string','m=2|m=4|m=6|m=8')
  set(MA_text,'visible','on')
  set(AR_menu,'visible','on')
  set(AR_text,'visible','on')
 case 5
  set(K_menu,'visible','off')
  set(M_menu,'visible','off')
  set(MA_menu,'visible','off')
  set(AR_menu,'visible','off')
  set(MA_text,'visible','off')
  set(AR_text,'visible','off')
end
set(next_button,'enable','on')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trigger_call(~,~)
global next_button
set(next_button,'enable','on')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quit_call(~,~)
global next_button
global stop
stop=1;
set(next_button,'enable','on')
return