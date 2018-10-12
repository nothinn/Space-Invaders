function h=plot_ball(x,y,nc,c)
global ha xd yd;
h=fill(x+xd(:,nc),y+yd(:,nc),c,'parent',ha);
