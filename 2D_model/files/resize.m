function resize
global ha kkx kky
drawnow;
% calculate kk, kk is scale coefficient between pixels and units of axes:
set(ha,'units','pixels');
ap=get(ha,'position');
set(ha,'units','normalized');

xlm=get(ha,'Xlim');
xun=xlm(2)-xlm(1);
ylm=get(ha,'Ylim');
yun=ylm(2)-ylm(1);

aba=get(ha,'PlotBoxAspectRatio');

if (aba(1)/aba(2))>(ap(3)/ap(4))
    kkx=xun/ap(3);
    ap4s=aba(2)*ap(3)/aba(1);
    kky=yun/ap4s;
else
    ap3s=aba(1)*ap(4)/aba(2);
    kkx=xun/ap3s;
    kky=yun/ap(4);
end