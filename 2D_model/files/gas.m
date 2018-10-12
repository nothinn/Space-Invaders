function gas
clear;
clear global;
global rad lx ly m n;
global ha hf kkx kky xd yd
global x y
global fx % number of fixed
global n vx vy dt tmot hs
global hsl hslm xh yh th % hystory
global stpf
global nesf
global newf
newf=false; % new flag, if want one more time

nesf=false; % not enought space flag
                    while 1 % to have possible repeat all if newf=true
% inital properties window;
while 1 % repeat until random ball positions will be generated
    hfip=initial_props;
    hnds=guidata(hfip); % read handles structure

    stpf=false; % stop flag
    hslm=5; % maxiaml lenght
    xh=[];
    yh=[];
    th=[];
    hsl=0;
    fx=0; 
    n=str2num(get(hnds.n,'string'));	% number of balls
    dt=0.1;             %time step
    tmot=0;              %curent time
    lx=str2num(get(hnds.lx,'string')); ly=str2num(get(hnds.ly,'string'));      % box size
    ihr=get(hnds.hr,'value'); % if high randomization
    if ihr
        rndf=0.8; % randomiztion factor
    else
        rndf=0.4; % randomiztion factor
    end
    
    % radii
    rrr=str2num(get(hnds.rad,'string'));
    if get(hnds.rr,'value');
        rad=rrr*(1+rndf*rand(1,n));	% bals radii
    else
        rad=rrr*ones(1,n);	% bals radii
    end
    
    % masses
    mmm=str2num(get(hnds.m,'string'));
    if get(hnds.rm,'value');
        m=mmm*(1+rndf*rand(1,n));	% bals radii
    else
        m=mmm*ones(1,n);	% bals radii
    end
    
    if get(hnds.rpm,'value');
        % radii proportianal to mass
        rad=rrr*m/mean(m);
    end

    if get(hnds.rc,'value');
        col=[0.8*rand(n,1),0.8*rand(n,1),0.8*rand(n,1)]; % balls colors
    else
        col=[0.1*ones(n,1),0.8*ones(n,1),0.1*ones(n,1)];
    end
    
    close(hnds.figure1);
    
    out=set_random(30);  % random balls position
    if nesf % not enought space for balls
        nes_warn;
    else
        break;
    end
end
x=out(1,:);	      % coordinaties
y=out(2,:);
vx=out(3,:);   % initial velocities
vy=out(4,:);
%h=line(x,y); % 
hs=zeros(1,n); % handles
pr=12; % plot circle resolution
xd=zeros(pr,n); % ball line coordinaties
yd=zeros(pr,n);
al=0:2*pi/pr:2*pi-2*pi/pr; % angles for draw ball line
hf=figure;
set(hf,'NumberTitle','off','Name','2d gas');
ha=axes;
title('hint: move balls with mouse, click it');
set(ha,'Nextplot','add');
for nc=1:n
    xd(:,nc)=rad(nc)*cos(al);
    yd(:,nc)=rad(nc)*sin(al);
    hs(nc)=plot_ball(x(nc),y(nc),nc,col(nc,:));
    set(hs(nc),'UserData',nc); % put number to userdata
    set(hs(nc),'UserData',nc);
    set(hs(nc),'ButtonDownFcn',['bdf(' num2str(hs(nc),'%20.20f') ')']);
end
%set(h,'LineStyle','none','color','b','Marker','o','MarkerSize',2.4*rad(1),'EraseMode','xor');
% axis([0, lx, 0, ly]);
% axis('square');
%pause;
xlim([0 lx]);
ylim([0 ly]);
axis equal;
xlim([0 lx]);
ylim([0 ly]);
set(ha,'Units','normalized');
set(ha,'position',[0.05 0.03 0.65 0.9]);

uicontrol('parent',hf,'units','normalized','Style', 'pushbutton', 'String', 'stop',...
    'Position', [0.75 0.9 0.2 0.05], 'Callback', 'stp');

uicontrol('parent',hf,'units','normalized','Style', 'pushbutton', 'String', 'random velocities',...
    'Position', [0.75 0.8 0.2 0.05], 'Callback', 'rnd_v');

uicontrol('parent',hf,'units','normalized','Style', 'pushbutton', 'String', 'accelerate',...
    'Position', [0.75 0.7 0.2 0.05], 'Callback', 'acc');

uicontrol('parent',hf,'units','normalized','Style', 'pushbutton', 'String', 'decelerate',...
    'Position', [0.75 0.6 0.2 0.05], 'Callback', 'dec');

uicontrol('parent',hf,'units','normalized','Style', 'pushbutton', 'String', 'all to right',...
    'Position', [0.75 0.5 0.2 0.05], 'Callback', 'ar');

uicontrol('parent',hf,'units','normalized','Style', 'pushbutton', 'String', 'all up',...
    'Position', [0.75 0.4 0.2 0.05], 'Callback', 'au');

uicontrol('parent',hf,'units','normalized','Style', 'pushbutton', 'String', 'another parameters...',...
    'Position', [0.75 0.27 0.2 0.05], 'Callback', 'm_new');

uicontrol('parent',hf,'units','normalized','Style', 'pushbutton', 'String', 'finish',...
    'Position', [0.75 0.1 0.2 0.05], 'Callback', 'fin');

resize; % to calulate kk itiially
set(hf,'ResizeFcn','resize');

dt=0;
tic;
while 1
    if stpf
        break;
    end
    if newf
        break;
    end
   tmot=tmot+dt;         
    balls;	% makes moving at time step dt
    for nc=1:n
        set(hs(nc),'XData', x(nc)+xd(:,nc),'YData',y(nc)+yd(:,nc));
    end
    drawnow;
    dt=toc-tmot;
    if fx~=0
        % hystory
        if hsl<=hslm
            xh=[xh x(fx)];
            yh=[yh y(fx)];
            th=[th toc];
            hsl=hsl+1;
        else
            xh=[xh(2:end) x(fx)];
            yh=[yh(2:end) y(fx)];
            th=[th(2:end) toc];
            dth=th(end)-th(1);
            vx(fx)=(xh(end)-xh(1))/dth;
            vy(fx)=(yh(end)-yh(1))/dth;
        end
    end
 end;

drawnow;
close(hf);
drawnow;

if newf
    newf=false;
else
    newf=false;
    break;
end

                    end