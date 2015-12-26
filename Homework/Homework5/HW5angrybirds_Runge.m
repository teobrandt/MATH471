%MATH471: Homework 5 Angry Birds
close all
clear
clc

global gamma1 C gamma2 kappa rho delta N

%Initial location


N = 10;
gamma1 = 20;
gamma2 = 30;
kappa = 1;
rho = 1;
delta = 0.001;

t0 = 0;
tf = 5;
h = 0.01;
t = t0:h:tf;

B = zeros(2*numel(t),N);
for n = 1:N
    B(1,n) = 3*(rand(1)-0.5);
    B(2,n) = 3*(rand(1)-0.5);
end
C = c_t(t);
for t_ = 2:numel(t)
    index = t_*2-1;
    old = index - 2;
    k1 = gamma1*(C(:,t_ - 1) - B(old:(old+1),1));
    k2 = B(old:(old+1),1) + 0.5*k1*h;
    k3 = B(old:(old+1),1) + 0.5*k2*h;
    k4 = B(old:(old+1),1) + k3*h;
    B(index:(index+1),1) = B(old:(old+1),1) + (1/6)*(k1+2*k2+2*k3+k4)*h;
    
    B_last = B(old:(old+1),:);
    Bbar = B_bar(B_last);
    for n = 2:N
        k1 = gamma2*(B(old:(old+1),1)-B(old:(old+1),n)) + kappa*(Bbar - B(old:(old+1),n)) + repel(B(old:(old+1),n),B(old:(old+1),2:N));
        k2 = B(old:(old+1),n) + 0.5*k1*h;
        k3 = B(old:(old+1),n) + 0.5*k2*h;
        k4 = B(old:(old+1),n) + k3*h;
        B(index:(index+1),n) = B(old:(old+1),n) + (1/6)*(k1+2*k2+2*k3+k4)*h;
    end
end
count = 1;
for t_ = 1:numel(t)
    plot(C(1,t_),C(2,t_),'xr','markersize',10)
    axis([-2.5 2.5 -2.5 2.5])
    hold on
    for n = 1:N
        index = t_*2-1;
        plot(B(index,n),B(index+1,n),'.','markersize',(N+1-n))
    end
    
    pause(0.001)
    M(count) = getframe(figure(1));
    hold off
    count = count + 1;
end

vidObj = VideoWriter('angrybirds.avi');
open(vidObj);
writeVideo(vidObj,M);
close(vidObj);