%-------------------------------------------%
%---Assignment 1 of Computational Finance---%
%-------------------------------------------%

clear all;
close all;

%% setting parameters

s0 = 14; 
K = 15;
r = 0.1;
sigma = 0.25;
T = 0.5;
gamma = 1;
num_t=3000; %number of time steps
M=2000; %number of sample paths


dt=T/num_t;
range = 0:dt:T;

%% MC simulation
for j = 1:M
    S_p(j,1) = s0;
    V(j,1) = max(S_p(j,1) - K, 0);
    
    S_m(j,1) = s0;
    VP(j,1) = max(S_m(j,1) - K, 0);
    

    for i = 2:num_t+1
        
        rand_num=randn();
        S_p(j,i) = S_p(j,i-1) + r*S_p(j,i-1)*dt + sigma*(S_p(j,i-1)^gamma)*rand_num*sqrt(dt);
        S_m(j,i)=S_m(j,i-1) + r*S_m(j,i-1)*dt + sigma*(S_m(j,i-1)^gamma)*(-rand_num)*sqrt(dt);
        V(j,i) = max(S_p(j,i)-K, 0);
        VP(j,i)=(max(S_p(j,i)-K, 0)+max(S_m(j,i)-K, 0))/2;
    end
    
end

EqV = mean(V);
V0 = exp(-r*T) * EqV;

EqVP = mean(VP);
V0P = exp(-r*T) * EqVP;

%% calculate exact solution
for i = 1:num_t+1
    
    bigT = range(i);%T-range(i); 
    
    teachers(i)=bsexact(sigma, r, K, bigT, s0);
    
end

%% plot result
plot(V0,'r')
hold on
plot(V0P,'b')
hold on
plot(teachers)