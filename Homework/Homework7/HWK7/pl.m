clear;close all;clc
a = dir('u0*.txt');
x = load('x0000000.txt');
y = load('y0000000.txt');
fig = figure;
aviobj = avifile('exhaustsim.avi','compression','None');
for i = 1:5:length(a)
    u = load(a(i).name);
    clf
    contour(x,y,u,linspace(-1.5,1.5,30))
    hold on   
    plot(x(:,1),y(:,1),'k',...
         x(:,end),y(:,end),'k',...
         x(1,:),y(1,:),'k',...
         x(end,:),y(end,:),'k','linewidth',2)
    axis equal
    colorbar
    drawnow
    aviobj = addframe(aviobj, getframe(fig));
end

aviobj = close(aviobj);
