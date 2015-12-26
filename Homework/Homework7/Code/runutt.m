clear
clc

nrl = 11;
nsl = nrl;
hr = 2/(nrl-1);
hs = 2/(nsl-1);
r = zeros(nrl+2,1);

for i = 1:nrl+2
    r(i) = -1+(i-2)*hr;
end
s = zeros(nsl+2,1);
for i = 1:nsl+2
    s(i) = -1+(i-2)*hs;
end
x = s+r;
y = r-s;
rx = 0.5*ones(nrl+2,nsl+2);
ry = 0.5*ones(nrl+2,nsl+2);
sx = 0.5*ones(nrl+2,nsl+2);
sy = -0.5*ones(nrl+2,nsl+2);
u = zeros(nrl+2,nsl+2);
a = ones(nrl+2,nsl+2);
J = rx.*sy-ry.*sx;
for i = 1:nrl+2
    for j = 1:nsl+2
        u(i,j) = x(i)*x(i)+y(j)*y(j);%*y(j);%*x(i)+1*y(j)*y(j);
    end
end
mesh(x,y,u)

u_tt = uttsolver(u,rx,ry,sx,sy,a,hr,hs,J,nrl,nsl)

uttsolver(u,rx,ry,sx,sy,a,hr,hs,J,nrl,nsl)