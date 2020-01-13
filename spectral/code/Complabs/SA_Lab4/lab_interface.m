function lab_interface
global next_button
global stop
global method_menu
global part_menu
global m
global m_menu
global K_menu
global Km_text
global w_menu
global N
global w_text

N=64;
R=50;
f1=0.2*2*pi;
f2=f1+1/N*2*pi;
phi1=0;
phi2=0;
sig2=1;
e=sqrt(sig2)*randn(R,N);
y=ones(R,1)*(10*sin(f1*(0:N-1)+phi1)+5*sin(f2*(0:N-1)+phi2))+e;
r_k=50*cos(f1*(0:N-1))+12.5*cos(f2*(0:N-1));
r_k(1)=r_k(1)+sig2;

gui_fig=figure('position',[300 200 800 800]);
offset=[0 0.05 0 0];
plot_axes=axes('position',[0.075 0.075 0.9 0.775]);
plot_axes2=axes('position',[0.5625 0.5  0.4125 0.35],'visible','off');
plot_axes3=axes('position',[0.075 0.075 0.4125 0.35],'visible','off');
plot_axes4=axes('position',[0.5625 0.075 0.4125 0.35],'visible','off');

quit_pos=[0.9,0.9,0.1,0.1];
next_pos=[0.875 0.875 0.1 0.1];
Km_pos=[0.5 0.9, 0.15 0.05];
part_pos=[0.15 0.9 0.15 0.05];
method_pos=[0.325 0.9 0.15 0.05];
w_pos=[0.675 0.9 0.15 0.05];

text_axes=axes('position',[0.15 1 0.65 0.1],'visible','off');
text(0,-0.3,'Part')
text(0.27,-0.3,'Method')
Km_text=text(0.54,-0.3,'K=');
w_text=text(0.81,-0.3,'\omega=');

part_menu=uicontrol(gui_fig,'style','popup','string','(a)|(b)+(c)',...
'units','normalized','position',part_pos,'callback',@part_call);
quit_button=uicontrol(gui_fig,'style','pushbutton','String','Quit',...
'units','normalized','position',quit_pos,'callback',@quit_call);

next_button=uicontrol(gui_fig,'style','pushbutton','String','Next',...
'units','normalized','position',next_pos,'visible','off');

K_menu=uicontrol(gui_fig,'style','popup','String','1|4','units',...
'normalized','position',Km_pos,'callback',@trigger_call);

method_menu=uicontrol(gui_fig,'style','popup','String',...
'Slepian RFB|CAPON','units','normalized','position',method_pos,...
'callback',@method_call);

w_menu=uicontrol(gui_fig,'style','popup','string','0|0.1*2 \pi',...
'units','normalized','position',w_pos,'callback',@trigger_call);

m_menu=uicontrol(gui_fig,'style','popup','string',...
'N/4|N/2-1','units','normalized','position',Km_pos,'callback',...
@trigger_call,'visible','off');

K=[1 4];
m=[N/4 N/2-1];
stop=0;
L_=2^12;
while stop~=1
  set(plot_axes,'position',[0.075 0.075 0.9 0.775])
  axes(plot_axes)
  cla
  axes(plot_axes2)
  cla
  axes(plot_axes3)
  cla
  axes(plot_axes4)
  cla
  set(plot_axes2,'visible','off')
  set(plot_axes3,'visible','off')
  set(plot_axes4,'visible','off')
  switch get(part_menu,'value')
   case 1
    switch get(method_menu,'value')
     case 1
      switch get(K_menu,'value')
       case 1
        w_=(get(w_menu,'value')-1)*2*pi*0.1;
        h=slepian(N,1,1).*exp(1j*w_*(0:N-1).');
        axes(plot_axes)
        H=fftshift(abs(fft(h,L_)));
        H=H/max(max(H));
        HdB=20*log10(eps+H);
        plot((0:(L_-1))/L_*2*pi-pi,HdB)
        axis([-pi pi -60 0 ])
        hold on
        plot([f1 f1; -f1 -f1;f2 f2;-f2 -f2].',[-60 0;-60 0;-60 0;-60 0].','r*-')
        plot([w_ w_],[-60 0],'k*-')
        zoom on
        hold off
        xlabel('\omega')
        ylabel('dB')
        title('First Slepian')
       case 2
        w_=(get(w_menu,'value')-1)*2*pi*0.1;
        h=slepian(N,4,4).*repmat(exp(1j*w_*(0:N-1).'),1,4);
        H=fftshift(abs(fft(h,L_,1)));
        H=H/max(max(H));
        HdB=20*log10(eps+H);
        set(plot_axes,'position',[0.075 0.5 0.4125 0.35])
        set(plot_axes2,'visible','on')
        set(plot_axes3,'visible','on')
        set(plot_axes4,'visible','on')
        axes(plot_axes)
        plot((0:(L_-1))/L_*2*pi-pi,HdB(:,3))
        axis([-pi pi -60 0])
        hold on
        plot([f1 f1; -f1 -f1;f2 f2;-f2 -f2 ].',[-60 0;-60 0;-60 0;-60 0].','r*-')
        plot([w_ w_],[-60 0],'k*-')
        zoom on
        hold off
        xlabel('\omega')
        ylabel('dB')
        title('First Slepian')
        axes(plot_axes2)
        plot((0:(L_-1))/L_*2*pi-pi,HdB(:,4))
        axis([-pi pi -60 0])
        hold on
        plot([f1 f1; -f1 -f1;f2 f2;-f2 -f2 ].',[-60 0;-60 0;-60 0;-60 0].','r*-')
        plot([w_ w_],[-60 0],'k*-')
        zoom on
        hold off
        xlabel('\omega')
        ylabel('dB')
        title('Second Slepian')
        axes(plot_axes3)
        plot((0:(L_-1))/L_*2*pi-pi,HdB(:,1))
        axis([-pi pi -60 0])
        hold on
        plot([f1 f1; -f1 -f1;f2 f2;-f2 -f2 ].',[-60 0;-60 0;-60 0;-60 0].','r*-')
        plot([w_ w_],[-60 0],'k*-')
        zoom on
        hold off
        xlabel('\omega')
        ylabel('dB')
        title('Third Slepian')
        axes(plot_axes4)
        plot((0:(L_-1))/L_*2*pi-pi,HdB(:,2))
        axis([-pi pi -60 0])
        hold on
        plot([f1 f1; -f1 -f1;f2 f2;-f2 -f2 ].',[-60 0;-60 0;-60 0;-60 0].','r*-')
        plot([w_ w_],[-60 0],'k*-')
        zoom on
        hold off
        xlabel('\omega')
        ylabel('dB')
        title('Fourth Slepian')
      end
     case 2
      m=[N/4 N/2-1];
      m_=m(get(m_menu,'value'));
      w_=(get(w_menu,'value')-1)*2*pi*0.1;
      a=exp(-1j*(0:m_).'*w_);
      R_inv=inv(toeplitz(r_k(1:m_+1)));
      h=R_inv*a/(a'*R_inv*a);
      H=fftshift(fft(conj(h),L_));
      H=H/max(abs(H));
      axes(plot_axes)
      plot(2*pi*(0:L_-1)/L_-pi,20*log10(abs(H)+eps))
      axis([-pi pi -60 0])
      hold on
      plot([f1 f1; -f1 -f1;f2 f2;-f2 -f2 ].',[-60 0;-60 0;-60 0;-60 0].','r*-')
      plot([w_ w_],[-60 0],'k*-')
      zoom on
      hold off
      xlabel('\omega')
      ylabel('dB')
    end
   case 2
    phi=zeros(R,L_);
    switch get(method_menu,'value')
     case 1
	  K_=K(get(K_menu,'value'));
      for k=1:R
	   phi(k,:)=rfb(y(k,:),K_,L_).';
      end
      phi=fftshift(phi,2);
     case 2
      m_=m(get(m_menu,'value'));
      for k=1:R
	   phi(k,:)=capon(y(k,:),m_,L_).';
      end
      phi=fftshift(phi,2);
     case 3
      m_=m(get(m_menu,'value'));
      a=complex(zeros(m_+1,R));
      ex=exp(-1j*(0:m_).'*(2*pi*(0:L_-1)/L_-pi));
      sig2_=zeros(R,1);
      for k=1:R
	   [a(:,k),sig2_(k)]=lsar(y(k,:),m_);	
      end
      disp([num2str(size(sig2_)) ' ' num2str(R) ' ' num2str(m_)])
      Aw=a.'*ex;
      phi=repmat(sig2_,1,L_)./(Aw.*conj(Aw));
    end 
    axes(plot_axes)
    phi=phi/max(max(phi));
    plot(2*pi*(0:L_-1)/L_-pi,10*log10(eps+phi),'b')
    axis([-pi pi -60 0])
    hold on
    plot([f1 f1; -f1 -f1;f2 f2;-f2 -f2 ].',[-60 0;-60 0;-60 0;-60 0].','r*-')
    hold off
    xlabel('\omega')
    ylabel('dB')
    title('Estimates from 50 Monte-Carlo realizations')
  end
  set(next_button,'enable','off')
  waitfor(next_button,'enable','on')
end

close(gui_fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trigger_call(~,~)
global next_button
set(next_button,'enable','on')
return

function part_call(~,~)
global next_button
global part_menu
global method_menu
global m
global N
global w_menu
global w_text
global m_menu
global K_menu
global Km_text
switch get(part_menu,'value')
 case 1
  set(w_menu,'visible','on')
  set(w_text,'visible','on')
  set(method_menu,'string','Slepian RFB|CAPON')
  set(m_menu,'string','N/4|N/2-1')
  if get(m_menu,'value')>2
    set(m_menu,'value',1)
  end
  if get(method_menu,'value')>2
    set(method_menu,'value',1)
    set(K_menu,'visible','on')
    set(m_menu,'visible','off')
    set(Km_text,'String','K=')
  end
  m=[N/4 N/2-1];
 case 2
  set(w_menu,'visible','off')
  set(w_text,'visible','off')
  set(method_menu,'string','Slepian RFB|CAPON|LSAR')
  set(m_menu,'string','N/4|N/2-1|8|16|30')
  m=[N/4 N/2-1 8 16 30];
end
set(next_button,'enable','on')
return

function method_call(~,~)
global next_button
global method_menu
global K_menu
global m_menu
global Km_text
switch get(method_menu,'value')
 case 1
  set(K_menu,'visible','on')
  set(m_menu,'visible','off')
  set(Km_text,'String','K=')
 case 2
  set(m_menu,'visible','on')
  set(K_menu,'visible','off')
  set(Km_text,'String','m=')
 case 3
  set(m_menu,'visible','on')
  set(K_menu,'visible','off')
  set(Km_text,'String','m=n=')
end
set(next_button,'enable','on')
return

function quit_call(~,~)
global stop
global next_button
stop=1;
set(next_button,'enable','on')
return