%Practice Routine: Discretization in Space
function u_tt = uttsolver(u,rx,ry,sx,sy,a,hr,hs,J,nrl,nsl)
%input: Q is the product of J, partial, and a...
%output: second derivative of u w.r.t. time


%Compute Q: all eight of them
%note one Q has the dimensions (nrl+2,nsl+2)
Q1 = zeros(nrl+2,nsl+2);
Q2 = zeros(nrl+2,nsl+2);
Q3 = zeros(nrl+2,nsl+2);
Q4 = zeros(nrl+2,nsl+2);
Q5 = zeros(nrl+2,nsl+2);
Q6 = zeros(nrl+2,nsl+2);
Q7 = zeros(nrl+2,nsl+2);
Q8 = zeros(nrl+2,nsl+2);

for m = 1:nrl+2
    for n = 1:nsl+2
        Q1(m,n) = J(m,n)*rx(m,n)*a(m,n)*rx(m,n);
        Q2(m,n) = J(m,n)*rx(m,n)*a(m,n)*sx(m,n);
        Q3(m,n) = J(m,n)*ry(m,n)*a(m,n)*ry(m,n);
        Q4(m,n) = J(m,n)*ry(m,n)*a(m,n)*sy(m,n);
        Q5(m,n) = J(m,n)*sx(m,n)*a(m,n)*rx(m,n);
        Q6(m,n) = J(m,n)*sx(m,n)*a(m,n)*sx(m,n);
       Q7(m,n) = J(m,n)*sy(m,n)*a(m,n)*ry(m,n);
        Q8(m,n) = J(m,n)*sy(m,n)*a(m,n)*sy(m,n);
    end
end
u
Ju_tt = zeros(nrl,nsl);
for m = 2:nrl+1
    for n = 2:nsl+1
        Q = Q1;
        Ju_tt(m-1,n-1) = (1/(2*hr^2))*(u(m+1,n)*(Q(m+1,n)+Q(m,n))-u(m,n)*(Q(m+1,n)+2*Q(m,n)+Q(m-1,n))+u(m-1,n)*(Q(m,n)+Q(m-1,n)));
        Q = Q2;
        Ju_tt(m-1,n-1) = Ju_tt(m-1,n-1) + (1/(4*hr*hs))*(Q(m,n+1)*(u(m+1,n+1)-u(m+1,n-1))-Q(m,n-1)*(u(m+1,n-1)-u(m-1,n-1)));
        Q = Q3;
        Ju_tt(m-1,n-1) = Ju_tt(m-1,n-1) + (1/(2*hr^2))*(u(m+1,n)*(Q(m+1,n)+Q(m,n))-u(m,n)*(Q(m+1,n)+2*Q(m,n)+Q(m-1,n))+u(m-1,n)*(Q(m,n)+Q(m-1,n)));
        Q = Q4;
        Ju_tt(m-1,n-1) = Ju_tt(m-1,n-1) + (1/(4*hr*hs))*(Q(m,n+1)*(u(m+1,n+1)-u(m+1,n-1))-Q(m,n-1)*(u(m+1,n-1)-u(m-1,n-1)));
        Q = Q5;
        Ju_tt(m-1,n-1) = Ju_tt(m-1,n-1) + (1/(4*hr*hs))*(Q(m+1,n)*(u(m+1,n+1)-u(m+1,n-1))-Q(m-1,n)*(u(m-1,n+1)-u(m-1,n-1)));
        Q = Q6;
        Ju_tt(m-1,n-1) = Ju_tt(m-1,n-1) + (1/(2*hs^2))*(u(m,n+1)*(Q(m,n+1)+Q(m,n))-u(m,n)*(Q(m,n+1)+2*Q(m,n)+Q(m,n-1))+u(m,n-1)*(Q(m,n)+Q(m,n-1)));
        Q = Q7;
        Ju_tt(m-1,n-1) = Ju_tt(m-1,n-1) + (1/(4*hr*hs))*(Q(m+1,n)*(u(m+1,n+1)-u(m+1,n-1))-Q(m-1,n)*(u(m-1,n+1)-u(m-1,n-1)));
        Q = Q8;
        Ju_tt(m-1,n-1) = Ju_tt(m-1,n-1) +...
            (1/(2*hs^2))*(u(m,n+1)*(Q(m,n+1)+Q(m,n))...
            -u(m,n)*(Q(m,n+1)+2*Q(m,n)+Q(m,n-1))...
            +u(m,n-1)*(Q(m,n)+Q(m,n-1)));
    end
end

u_tt = Ju_tt./J(2:nrl+1,2:nsl+1);
