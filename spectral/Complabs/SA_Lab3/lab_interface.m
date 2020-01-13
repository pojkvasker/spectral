function lab_interface
global next_button
global stop
global method_menu
global ML_menu
global m_menu
global ML_text
global sig2_text
global method_text
global input_text
global m_text
global gui_fig

N=64;
R=50;
n=4;
sig_2=1;
f1=0.24*pi;
f2=0.26*pi;
w_true=[-f1 f1 -f2 f2];
k=0:24;
ACS_true=50*cos(f1*k)+12.5*cos(f2*k);
e=randn(R,N);

gui_fig=figure('Position',[500 400 820 440]);
plot_axes=axes('position',[0.1 0.12 0.7 0.75]);
text_axes=axes('position',[0.1 0.12 0.7 0.75],'visible','off');
method_text=text(1.05,1,'Method:');
input_text=text(1.05,1.13,'Input:');
ML_text=text(1.05,0.86,'M=L=');
m_text=text(1.05,0.86,'m=','Visible','off');
sig2_text=text(1.05,0.735,'\sigma^2=');

ML_pos=[0.825 0.7 0.15 0.05];
m_pos=[0.825 0.7 0.15 0.05];
method_pos=[0.825 0.8 0.15 0.05];
next_pos=[0.875 0.025 0.1 0.1];
quit_pos=[0.875 0.025 0.1 0.1];
input_pos=[0.825 0.9 0.15 0.05];
sig2_pos=[0.825 0.6 0.15 0.05];
method_menu=uicontrol(gui_fig,'style','popup', 'string',...
'HOYW|MinNorm|MUSIC|ESPRIT','callback',@method_call,'value',1,...
'units','normalized','position',method_pos);

input_menu=uicontrol(gui_fig,'style','popup', 'string',...
'True ACS |y(t) ','callback',@trigger_call,'value',1,...
'units','normalized','position',input_pos);

next_button=uicontrol(gui_fig,'style','pushbutton', 'string',...
'Next','units','normalized','position',next_pos,'Visible','off');

quit_button=uicontrol(gui_fig,'style','pushbutton', 'string',...
'Quit','callback',@quit_call,...
'units','normalized','position',quit_pos);

ML_menu=uicontrol(gui_fig,'style','popup', 'string',...
'4|8|12','callback',@trigger_call,'value',1,...
'units','normalized','position',ML_pos);

m_menu=uicontrol(gui_fig,'style','popup', 'string',...
'5|8|12','callback',@trigger_call,'value',1,...
'units','normalized','position',m_pos,'visible','off');

sig2_menu=uicontrol(gui_fig,'style','popup', 'string',...
'1|0.125|0.0125|0','callback',@trigger_call,'value',4,...
'units','normalized','position',sig2_pos);

stop=0;
m=[5 8 12];
ML=[4 8 12];
sig2=[1 0.125 0.0125 0];
while stop~=1
  switch get(input_menu,'value')
   case 1
    switch get(method_menu,'value')
     case 1
      ML_=ML(get(ML_menu,'value'));
      sig2_=sig2(get(sig2_menu,'value'));
      ACS=ACS_true;
      ACS(1)=ACS(1)+sig2_;
      w=hoyw_r(ACS(2:end),n,ML_,ML_);
     case 2
      m_=m(get(ML_menu,'value'));
      sig2_=sig2(get(sig2_menu,'value'));
      ACS=ACS_true;
      ACS(1)=ACS(1)+sig2_;
      w=minnorm_r(ACS(1:m_),n);
     case 3
      m_=m(get(ML_menu,'value'));
      sig2_=sig2(get(sig2_menu,'value'));
      ACS=ACS_true;
      ACS(1)=ACS(1)+sig2_;
      w=music_r(ACS(1:m_),n);
     case 4
      m_=m(get(ML_menu,'value'));
      sig2_=sig2(get(sig2_menu,'value'));
      ACS=ACS_true;
      ACS(1)=ACS(1)+sig2_;
      w=esprit_r(ACS(1:m_),n);
    end
    axes(plot_axes)
    p1=plot([w(:).';w(:).'],[zeros(1,n);ones(1,n)],'k','LineWidth',3);
    hold on
    p2=plot([w_true(:).';w_true(:).'],[zeros(1,n);ones(1,n)],'*r-');
    set(plot_axes,'ytick',[])
    hold off
    xlabel('\omega')
    title('Estimates using true ACS sequences')
    zoom on
    legend([p1(1) p2(1)],'Estimated frequencies','True Frequencies')
    axis([-pi pi 0 1])
   case 2
    t=0:N-1;
    sig2_=sig2(get(sig2_menu,'value'));
    y=10*sin(repmat(f1*t,R,1)+repmat(2*pi*rand(R,1),1,N))+...
      5*sin(repmat(f2*t,R,1)+repmat(2*pi*rand(R,1),1,N))+...
    sqrt(sig2_)*e;
    w=zeros(n,R);
    switch get(method_menu,'value')
     case 1      
      ML_=ML(get(ML_menu,'value'));
      for k=1:R
	   w(:,k)=hoyw(y(k,:),n,ML_,ML_);
      end    
     case 2
      m_=m(get(m_menu,'value'));
      for k=1:R
       w(:,k)=minnorm(y(k,:),n,m_);
      end    
     case 3
      m_=m(get(m_menu,'value'));
      for k=1:R
	   w(:,k)=music(y(k,:),n,m_);
      end
     case 4
      m_=m(get(m_menu,'value'));
      for k=1:R
	   w(:,k)=esprit(y(k,:),n,m_);
      end
     
    end
    for k=find(w>pi)
      w(k)=w(k)-2*pi;
    end
    for k=find(w<-pi)
      w(k)=w(k)+2*pi;
    end
    axes(plot_axes)
    hist(w(:),-pi:0.01*pi/3:pi);
    p1=findobj(gca,'Type','patch');
    set(p1,'FaceColor','k','EdgeColor','k')
    ylabel('number of estimates in frequency bin')
    xlabel('\omega')
    hold on
    y_max=get(plot_axes,'ylim');
    y_max=y_max(2);
    p2=plot([w_true(:).';w_true(:).'],[zeros(1,n);y_max*ones(1,n)],'*r-');
    title('histogram')
    %legend(p1,'True frequencies')
    legend('Estimated frequencies','True frequencies')
    axis([-pi pi 0 y_max])
    
    zoom on
    hold off
  end
  set(next_button,'enable','off')
  waitfor(next_button,'enable','on')
end
close(gui_fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function method_call(~,~)
global method_menu
global ML_menu
global m_menu
global ML_text
global m_text
global next_button

switch get(method_menu,'value')
 case 1
  set(m_menu,'visible','off')
  set(m_text,'Visible','off')
  set(ML_menu,'visible','on')
  set(ML_text,'Visible','on')
 case 2
  set(ML_menu,'visible','off')
  set(ML_text,'Visible','off')
  set(m_menu,'visible','on')
  set(m_text,'Visible','on')
 case 3
  set(ML_menu,'visible','off')
  set(ML_text,'Visible','off')
  set(m_menu,'visible','on')
  set(m_text,'Visible','on')
 case 4
  set(ML_menu,'visible','off')
  set(ML_text,'Visible','off')
  set(m_menu,'visible','on')
  set(m_text,'Visible','on')
end
set(next_button,'enable','on')
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trigger_call(~,~)
global next_button
set(next_button,'enable','on')
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quit_call(~,~)
global stop
global next_button
stop=1;
set(next_button,'enable','on')
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
