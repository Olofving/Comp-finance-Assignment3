function [dt,ds,s,last_v] = fdemplicit(K,r,sigma,T,tn,sn)

%% function to compute the call option price use emplicit euler 
% parameter
% K = 15; strike price 
% r = 0.1; risk free interest rate
% sigma = 0.25; volatility
% T = 0.5; option time 
% tn: number of time points on discretization, and dt = (T - T0)/(tn-1)
% sn: number of price points on discretizatin, and ds = (S - S0)/(sn-1)

Smax = 4 * K; %the price range to compute
t = linspace(0,T,tn); % time vector
s = linspace(0,Smax,sn); % price vector
dt = t(2)-t(1); 
ds = s(2)-s(1); 
last_v = max(s-K,0); % boundary condition, at time = T, we know the value of option
last_v = last_v(:); %force a column vector

for n=tn:-1:2
    v(1) = 0;
    v(sn) = Smax-K*exp(-r*(T-t(n-1)));
    for j=2:sn-1
        v(j)=last_v(j)+r*s(j)*dt/(2*ds)*(last_v(j+1)-last_v(j-1))...
             +sigma^2*0.5*(s(j))^2*(dt/(ds^2))*(last_v(j+1)-...
             2*last_v(j)+last_v(j-1))-dt*r*last_v(j);
    end
    last_v=v;
end
end