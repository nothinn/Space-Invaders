function bdf(hfl)
% Button Down Function
% ha - axes handler
% hfl - fill handle (circle)

global hf x y fx vx vy
% 
% xys

if strcmpi(get(hf,'SelectionType'),'normal') % only if left click

    no=get(hfl,'Userdata');
    fx=no;
    vx(no)=0;
    vy(no)=0;

    % initial pointer location in pixels:
    pl=get(0,'PointerLocation');
    xp=pl(1);
    yp=pl(2);
    
    % initial position of ball
    x1=x(no);
    y1=y(no);


    % what to do when move:
    wbmf=['move(' num2str(hfl,'%20.20f') ',' num2str(xp,'%20.20f') ',' num2str(yp,'%20.20f') ',' num2str(x1,'%20.20f') ',' num2str(y1,'%20.20f') ');'];
    set(hf,'WindowButtonMotionFcn',wbmf);

    % what to do when release button:
    %wbuf=['set(' num2str(hf,'%20.20f') ',''WindowButtonMotionFcn'','''');']; % not do moving any more
    wbuf='buf';
    set(hf,'WindowButtonUpFcn',wbuf);
end