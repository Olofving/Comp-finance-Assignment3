%% Assignment 2 - Finite difference methods to price European Options
%
% Written by Peili Guo (Peili.Guo.7645@student.uu.se) and 
% Sijia Wang (Sijia.Wang.7090@student.uu.se)
% This report is for Computational Finance: Pricing and Valuation
% Assignment 2 Group 11. 
%
% In this project we implemented both explicit and implicit euler's method
% in Matlab to price a European call option. We studied the convergency,
% stability, computation cost and how the solution varies with \gamma. 
% 
clear all;
close all;

K = 15;
r = 0.1;
sigma = 0.25;
T = 0.5;
gamma = 1;
Smax = 4 * K;
tn = 6; %num of points in the time grid
sn = 61; % number of points in the space (price grid)
trange = [11 51 101 151 201];
srange = [61 121 181 241 301];
count = 0;
    % Euler explicit 6
%for i = 1: length(trange)
 
[dt1,ds1,se1,ve1] = fdemplicit(K,r,sigma,T,trange(1),sn);

[dt2,ds2,se2,ve2] = fdemplicit(K,r,sigma,T,trange(2),sn);

[dt3,ds3,se3,ve3] = fdemplicit(K,r,sigma,T,trange(3),sn);
tic
[dt4,ds4,se4,ve4] = fdemplicit(K,r,sigma,T,trange(4),sn);
toc

% Euler implicit
[dt1,ds1,si1,vi1] = fdimplicit(K,r,sigma,T,gamma,trange(1),sn);
[dt2,ds2,si2,vi2] = fdimplicit(K,r,sigma,T,gamma,trange(2),sn);
[dt3,ds3,si3,vi3] = fdimplicit(K,r,sigma,T,gamma,trange(3),sn);
tic
[dt4,ds4,si4,vi4] = fdimplicit(K,r,sigma,T,gamma,trange(4),sn);
toc
%exact solution
for i=1:sn
    exact1(i) = bsexact(sigma, r, K, T, si1(i));
    exact2(i) = bsexact(sigma, r, K, T, si2(i));
    exact3(i) = bsexact(sigma, r, K, T, si3(i));
    exact4(i) = bsexact(sigma, r, K, T, si4(i));
end
% compute the largest error
err_e1 = abs(ve1 - exact1);
err_i1 = abs(vi1 - exact1);
err_e2 = abs(ve2 - exact2);
err_i2 = abs(vi2 - exact2);
err_e3 = abs(ve3 - exact3);
err_i3 = abs(vi3 - exact3);
err_e4 = abs(ve4 - exact4);
err_i4 = abs(vi4 - exact4);
%end

% plot 1 showing the stability of the results using explicit and implicit
% method
f1 = figure('position', [0,0,1000,450]);
subplot(1,2,1);
plot(se1,ve1,'r*-')
hold on 
plot(se3,ve3,'k+-')
legend("dt =0.05", "dt = 0.005")
xlabel('Price');
ylabel('Value of Option');
title({'results using different dt size with Euler Explicit FDM';' '});

subplot(1,2,2);
plot(si1,vi1,'r*-')
hold on 
plot(si3,vi3,'k+-')
legend("dt =0.05", "dt = 0.005")
xlabel('Price');
ylabel('Value of Option');
title({'results using different dt size with Euler Implicit FDM',' '});


f2 = figure('position', [0, 0, 700, 500]);
plot(si1,err_i1,'g')
hold on 
plot(si2,err_i2,'b')
plot(si3,err_i3,'r')
plot(si4,err_i4,'k')
legend("dt = 0.05", "dt = 0.01", 'dt = 0.005', 'dt = 0.0025')
xlabel('Price');
ylabel('Absolute error');
title('Absolute error using different time step sizes (ds = 1)', 'FontSize', 18);
%% Accuracy using implicit method 
count = 0 
for tt = trange
    exactt=[];
    count = count+1
    [dtt,dst,sit,vi] = fdimplicit(K,r,sigma,T,gamma,tt,301);
    for i=1:length(sit)
        exactt(i) = bsexact(sigma, r, K, T, sit(i));
    end
    dtplot(count) = dtt;
    errtt(count) = max(abs(exactt - vi));
end
count = 0;
for ss = srange
    exactt=[];
    count = count+1
    [dtt,dst,sit,vi] = fdimplicit(K,r,sigma,T,gamma,201,ss);
    for i=1:length(sit)
        exactt(i) = bsexact(sigma, r, K, T, sit(i));
    end
    dsplot(count) = dst;
    errss(count) = max(abs(exactt - vi));
end

f4 = figure('position', [0,0,1000,450]);
subplot(1,2,1);
loglog(dtplot,errtt,'r*-')
xlabel('dt');
ylabel('max error');
title({'max error vs dt (ds = 0.2)';' '});

subplot(1,2,2);
loglog(dsplot,errss,'b*-')
xlabel('ds');
ylabel('max error');
title({'max error vs ds (dt = 0.0025)';' '});



%% experiment with different gamma 
count = 0;
for g = 0.5:0.1:1
    count = count + 1;
    [dt,ds,sg,vjg(count,:)] = fdimplicit(K,r,sigma,T,g,101,61)
end
f3 = figure('position', [0, 0, 700, 500]);
plot(sg,vjg)
xlim([10 20])
legend("gamma = 0.5", "gamma = 0.6", 'gamma = 0.7', 'gamma = 0.8','gamma = 0.9','gamma = 1','Location','northwest')
xlabel('Price');
ylabel('Value of Option');
title('V(s) changes when \gamma varies');


%% plots
% plot with different size of t, s.
% tictoc on time.
%plot(s,last_v,'r*')
% plot(s,vi,'gd')
% 
% hold on
% plot(s,exact,'b')
% plot(s,ve,'k*-')
% 
%  figure
%  plot(s,abs(vi'-exact),'r')
%  hold on
%  plot(s,abs(ve-exact),'b')
%  
%  %% task to be done 
%  % plot with different size of t, s.
 % tictoc on time.
 % to plot V(sj) with different \gamma for sj = 15. 
