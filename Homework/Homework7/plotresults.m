
ufiles = dir('u*.txt');
vfiles = dir('v*.txt');

load x.txt
load a.txt

for i = 1:length(ufiles)
    u = load(ufiles(i).name);    
    v = load(vfiles(i).name);
    
    plot(x,a,'b--',x,u,'k',x,v,'r')
    axis([-1 2 -4.1 4.1])
    drawnow
end

max(abs(u))
max(abs(v+cos(x)))