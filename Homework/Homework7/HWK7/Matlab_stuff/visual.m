clear;close all;clc

nx = 40;
ny = 40;
nr = nx;
ns = ny;
hr = 2/(nr-1);
hs = 2/(ns-1);

r = zeros(nr,1);
for ir = 1:nr
    r(ir) = -1+(ir-1)*hr;
end
s = zeros(ns,1);
for is = 1:ns
    s(is) = -1+(is-1)*hs;
end
%xcoord = @(r,s) -2+(2+r+0.2*sin(5*s))*cos(s);
%ycoord = @(r,s) (2+r+0.3*sin(0.5*s))*sin(s);

xcoord = @(r,s) 4*r+2*sin(r)*0.05*sin(r);
ycoord = @(r,s) sin(3*r)+2*sin(s);

%xcoord = @(r,s) r-3*s;
%ycoord = @(r,s) s+r;

x = zeros(nr,ns);
y = zeros(nr,ns);
for i = 1:nr
    for j = 1:ns
        x(i,j) = xcoord(r(i),s(j));
        y(i,j) = ycoord(r(i),s(j));
    end
end


plot(x,y,'k',x',y','k')
hold on
axis equal;%([-4.1 4.1 -3.1 3.1]
plot(xcoord(1,-1),ycoord(1,-1),'*r')
plot(xcoord(1,1),ycoord(1,1),'*b')
plot(xcoord(-1,1),ycoord(-1,1),'*k')
plot(xcoord(-1,-1),ycoord(-1,-1),'*m')
plot(xcoord(0,0),ycoord(0,0),'xr')
xlabel('x','fontsize',16)
ylabel('y','fontsize',16)
title('Test Mesh #2','fontsize',18)
%figure
%[R,S]=meshgrid(r,s);
%plot(R,S,'k',R',S','k')
%hold on
%axis equal
%plot(1,-1,'*r')
%lot(1,1,'*b')
%plot(-1,1,'*k')
%plot(-1,-1,'*m')
%plot(0,0,'xr')
%xlabel('r','fontsize',16)
%ylabel('s','fontsize',16)
