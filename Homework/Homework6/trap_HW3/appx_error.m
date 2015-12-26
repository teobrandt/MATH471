clear
close all
clc

gqdata = dlmread('gquad.txt');
gq_iteration = gqdata(:,1);
gq_f_x = gqdata(:,2);
gq_Appx_abs_error = gqdata(:,3);

gq2data = dlmread('gquad2.txt');
gq2_iteration = gq2data(:,1);
gq2_f_x = gq2data(:,2);
gq2_Appx_abs_error = gq2data(:,3);

tadata = dlmread('trappx.txt');
ta_iteration = tadata(:,1);
ta_f_x = tadata(:,2);
ta_Appx_abs_error = tadata(:,3);

ta2data = dlmread('trappx2.txt');
ta2_iteration = ta2data(:,1);
ta2_f_x = ta2data(:,2);
ta2_Appx_abs_error = ta2data(:,3);

C1 = 1.71;
C2 = 2;
alpha1 = 3;
alpha2 = -6;

loglog(gq_iteration, gq_Appx_abs_error, gq2_iteration, gq2_Appx_abs_error, ta_iteration, ta_Appx_abs_error, ta2_iteration, ta2_Appx_abs_error, ta2_iteration, C1.^(-alpha1.*ta2_iteration), ta2_iteration, C2.^(alpha2.*ta2_iteration));
axis([0, 1325, 10E-16, 10E1])
grid on;
title('Error of Trapezoid Rule and Gauss Quadrature'); xlabel('n'); ylabel('error'); legend('Gauss: Pi', 'Gauss: Pi^2', 'Trapezoid: Pi', 'Trapezoid: Pi^2', 'C_{1} = 1.71, a_{1} = 3', 'C_{2} = 2.0, a_{2} = -6');

print('ErrorPlot', '-dpng')

