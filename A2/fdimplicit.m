function [dt,ds,s,vj] = fdimplicit(K,r,sigma,T,g,tn,sn)

%% function to compute the call option price use implicit euler 
% parameter
% K = 15; strike price 
% r = 0.1; risk free interest rate
% sigma = 0.25; volatility
% T = 0.5; option time 
% g = gamma, and set gamma = 1 to inspect ;
% tn: number of time points on discretization, and dt = (T - T0)/(tn-1)
% sn: number of price points on discretizatin, and ds = (S - S0)/(sn-1)

Smax = 4 * K; %the price range to compute
t = linspace(0,T,tn); % time vector
s = linspace(0,Smax,sn); % price vector
dt = t(2)-t(1); 
ds = s(2)-s(1); 
vj = max(s-K,0); % boundary condition, at time = T, we know the value of option
vj = vj(:); %force a column vector
sig2 = sigma*sigma;
%construct the tridiagonal matrix
a = zeros(1,sn);
b = zeros(1,sn);
c = zeros(1,sn);
for i = 1:1:sn %iterate at each grid point
    a2(i) = 0.5*dt*(r*i - sig2*i^(2*g)*ds^(2*g-2));
    b2(i) = 1 + dt*(sig2*i^(2*g)*ds^(2*g-2) + r);
    c2(i) = -0.5*dt*(r*i + sig2*i^(2*g)*ds^(2*g-2));
    a(i) = 0.5*dt*(r*s(i)/ds - sig2*(s(i)^(2*g))/(ds*ds));
    b(i) = 1 + dt*(sig2*(s(i)^(2*g)/(ds*ds)) + r);
    c(i) = -0.5*dt*(r*s(i)/ds + sig2*(s(i)^(2*g)/(ds*ds)));
end
D = diag(a(3:end-1),-1) + diag(b(2:end-1)) + diag(c(2:end-2),1);
 
for n = tn:-1:2
    v = zeros(sn,1);
    %v(1) = 0; already initialized 
    v(sn,1) = Smax-K*exp(-r*(T-t(n-1)));    
    f = vj(2:end-1); % take vj(2:end-1) exclude the boundary points 
    f(end) = f(end) - c(end-1)*v(end); % for the v(t(i-1), s(n-1))
    v(2:end-1) = D\f;
    vj = v;
end
vj = vj';
end