%MATH471: Convergence with MMS
clear;close all;clc
%Error: test 1
n = [5 10 20 40 80 160];
n = n.^2;
error = [
    0.66631260214370269;
    3.40731237408417043E-002;
    8.92598114289569304E-003;
    1.85550587520385635E-003;
    4.18963690041942972E-004;
    9.97811317795268784E-005];
loglog(n,error,'-*')
set(gca,'fontsize',16)
grid on
title('Error Test','fontsize',22)
xlabel('log(number of grid points)','fontsize',20)
ylabel('log(avg. error)','fontsize',20)
hold on
error = [
    3.5498381642409003
    0.48355369389251002;
    0.10241675873930009;
    3.70942280545531744E-002;
    2.06178374928129501E-002;
    1.47894097887570734E-002];
loglog(n,error,'-*')
legend('test 1: simple domain','test 2: exhaust domain')