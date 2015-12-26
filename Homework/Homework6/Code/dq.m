function void = dq

format long 
format compact
% Choose number of dscretization points.
nr = 20;
ns = 1*nr;
% This results in a certain grid spacing.
hr = 2/nr;
hs = 2/ns;
% Set up grids in r and s
r = -1 + hr*(0:nr)';
s = -1 + hs*(0:ns)';
% Store them in arrays for conveninece 
R = zeros(nr+1,ns+1);
S = zeros(nr+1,ns+1);
for i = 1:nr+1
    for j = 1:ns+1
        R(i,j) = r(i);
        S(i,j) = s(j);
    end
end
% Choose some mapping. 
X = (R+2).*cos(S*pi/2);
Y = (R+2).*sin(S*pi/2);
% Plot the "phyical domain"
figure(1)
plot(X,Y,'k',X',Y','k')
axis equal

% Construct a made up integrand
RAD = sqrt(X.^2+Y.^2);
U = -(RAD-1).*(RAD-3);
figure(2)
surf(X,Y,U)
axis equal


% Construct the metric
XR = der(X,hr);
YR = der(Y,hr);
XS = der(X',hs)';
YS = der(Y',hs)';
JAC = XR.*YS-XS.*YR;

% Now integrate  
Integrand = U.*JAC;

% Version 1, integrate first in r to get a function of s alone.  
% N.B. this is the Riemann formula not Trapezoidal.
I_1Ds = zeros(ns+1,1);
for j = 1:ns+1
    I_1ds(j) = hr*sum(Integrand(1:nr,j));
end
% Then integrate in s
I1 = hs*sum(sum(I_1ds(1:ns)));

% Version 1, integrate first in r to get a function of s alone.  
% N.B. this is the Riemann formula not Trapezoidal.
I_1Ds = zeros(nr+1,1);
for i = 1:nr+1
    I_1ds(i) = hs*sum(Integrand(i,1:ns));
end
% Then integrate in s
I2 = hr*sum(I_1ds(1:nr));

% Exact value
Ie = 8/3*pi;

% Display results
disp(sprintf('Exact value:     %0.15g  \nVersion 1 value: %0.15g \nVersion 2 value: %0.15g',Ie,I1,I2))
disp(sprintf('\n(h_max,Error) : (%0.15e,%0.15e)',max(hr,hs),abs(I1-Ie)))
disp(sum(RAD))

function dw_dz = der(w,hz);
% Differentiates in the direction of the first index. 
% Apply with transpose to get differentiation wrt second index

dw_dz = zeros(size(w));

dw_dz(2:end-1,:) = (w(3:end,:)-w(1:end-2,:))/(2*hz);
dw_dz(1,:)   = (-w(3,:) + 4*w(2,:)-3*w(1,:))/(2*hz);
dw_dz(end,:) = -(-w(end-2,:) + 4*w(end-1,:)-3*w(end,:))/(2*hz);



