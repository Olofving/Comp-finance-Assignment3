% Computational Finance
% Finite Difference Method

% Group 4
clc 
clear

%% Parameter Values
M = 100;
N = 10000;

K = 15; 
r = 0.1;
sig = 0.25;
T = 0.5;
gamma = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];

s_max = 4*K;
s = linspace(0,s_max,M);

%% Exact - bsexact 
% Analytic solution for frice function
bsexa = bsexact(sig, r, K, T, s);

%% FDM
M = 100;
N = 10000;

v = zeros(M,N);
V = zeros(M,N);
V_imp = zeros(M,N);


t = linspace(0,T,N);
dt = T/N;
s = linspace(0,s_max,M);
ds = s_max/M;

v(:,N) = max(s-K,0);
V(:,N) = max(s-K,0);
V_imp(:,N) = max(s-K,0);

A = zeros(M,M);

%% Explicit- and Implicit compared to the analytical solution
[B, B_imp] = matrix_assembly(M,A,r,s,ds,dt,sig);
%[B_gamma,B_imp_gamma] = matrix_assembly_gamma(M,A_gamma,r,s,ds,dt,sig,gamma);

v = iterative(M,N,v,s_max,K,T,t,r,s,ds,dt,sig);
V = exp_euler(N,V,B,s_max,K,T,t,r);
V_imp = imp_euler(N,V_imp,B_imp,s_max,K,T,t,r);

figure
plot(s,V_imp(:,1)) 
hold on
plot(s,V(:,1))
plot(s,bsexa)
title('Explicit-, Implicit and analytical solution')
ylabel('Option Price')
xlabel('Stock Price')
legend('V_{explicit}','V_{implicit}','bsexact')
%%
% The implicit and the explicit solution seem to match the analytical
% solution exactly. 
%% Stability
M_stab = 100;
N_stab = 100;
t_stab = linspace(0,T,N_stab);
dt_stab = T/N_stab;
s_stab = linspace(0,s_max,M_stab);
ds_stab = s_max/M_stab;

V_stab = zeros(M_stab,N_stab);
V_stab(:,N_stab) = max(s_stab-K,0);
A_stab = zeros(M_stab,M_stab);
[B, B_imp] = matrix_assembly(M_stab,A_stab,r,s_stab,ds_stab,dt_stab,sig);
V_stab = exp_euler(N_stab,V_stab,B,s_max,K,T,t_stab,r);
figure
plot(s_stab,V_stab(:,1))
title('Unstable Explicit Solution')
ylabel('Option Price')
xlabel('Stock Price')
%%
% The explicit euler method is conditionally stable. The stability
% condition for this method is such that the time-step devided by the
% space-step squared has to be bigger than one. In the above unstable
% simulation of the option price, the time- and space steps are equal in
% number. 
%% Error plotting
figure
semilogy(s,abs(V(:,1)'-bsexa))
hold on
semilogy(s,abs(V_imp(:,1)'-bsexa))
title('Absolute Error')
ylabel('Absolute Error')
xlabel('Stock Price')
legend('V_{explicit} - bsexact','V_{implicit} - bsexact')
ylim([1e-5 1e-0])
%%
% The error of the explicit and implicit solution seem to match. The
% Simulation is run with 100 spacial-points and 10000 time-points.
%% Error as functions of space points
M_var = [10 50 100 500 1000 1500];
N_s = 5000;

for i = 1:length(M_var)
    t_s = linspace(0,T,N_s);
    dt = T/N_s;
    m = M_var(i);
    s_s = linspace(0,s_max,m);
    ds = s_max/m;
    
    bsexa = bsexact(sig, r, K, T, s_s);
    
    V_imp = zeros(m,N_s);
    V_imp(:,N_s) = max(s_s-K,0);
    A = zeros(m,m);
    [B, B_imp] = matrix_assembly(m,A,r,s_s,ds,dt,sig);
    V_imp = imp_euler(N_s,V_imp,B_imp,s_max,K,T,t,r);
    err_s(i) = max(abs(V_imp(:,1)'-bsexa));  
end
figure
loglog(M_var,err_s)
title('Error as functions of space points')
ylabel('Absolute Error')
xlabel('Space Points')
%% Error as functions of time points
M_s = 500;
N_var = [10 100 500 1000 1500];

for i = 1:length(N_var)
    m = M_s;
    n_var = N_var(i);
    t_t = linspace(0,T,n_var);
    dt = T/n_var;
    s_t = linspace(0,s_max,m);
    ds = s_max/m;
    
    bsexa = bsexact(sig, r, K, T, s_t);
    
    V_imp = zeros(m,n_var);
    V_imp(:,n_var) = max(s_t-K,0);
    A = zeros(m,m);
    [B, B_imp] = matrix_assembly(m,A,r,s_t,ds,dt,sig);
    V_imp = imp_euler(n_var,V_imp,B_imp,s_max,K,T,t_t,r);
    err_t(i) = max(abs(V_imp(:,1)'-bsexa));  
end
figure
loglog(N_var,err_t)
title('Error as functions of time points')
ylabel('Absolute Error')
xlabel('Time Points')

%% Option Price as a function of gamma
M_gamma = 1000;
N_gamma = 1000;
t_gamma = linspace(0,T,N_gamma);
dt_gamma = T/N_gamma;
s_gamma = linspace(0,s_max,M_gamma);
ds_gamma = s_max/M_gamma;
A_gamma = zeros(M_gamma,M_gamma);
V_imp_gamma = zeros(M_gamma,N_gamma);
V_imp_gamma(:,N) = max(s_gamma-K,0);
figure
for i = 1:length(gamma)
    g = gamma(i);
    [B_gamma,B_imp_gamma] = matrix_assembly_gamma(M_gamma,A_gamma,r,s_gamma,ds_gamma,dt_gamma,sig,g);
    V_imp_gamma = imp_euler(1000,V_imp_gamma,B_imp_gamma,s_max,K,T,t_gamma,r);
    plot(s_gamma,V_imp_gamma(:,1))
    hold on
end
title('Option Price as a function of gamma')
ylabel('Option Price')
xlabel('Stock price')
legend('0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1')
%%
% The gamma is summposed to act to decrease the volatility as the price of
% the stock increases. The decreased stock volatility for more valuble
% stock leads to a more stable stock and the growth more constant. This in
% turn leads to a more predictable and 'possibly' a higher stock increase.
% This causes the shape of the option price in the figure above. 


%% Explicit Euler
function V = exp_euler(N,V,B,s_max,K,T,t,r)
    for i = 0:N-2
        n = N-i;
        V(:,n-1) = B*V(:,n);
        V(1,n-1) = 0;
        V(end,n-1) = s_max - K*exp(-r*(T-t(n-1)));
        %disp(sum(sum(v-V)))
    end
end

function v = iterative(M,N,v,s_max,K,T,t,r,s,ds,dt,sig)
    for i = 0:N-2
        n = N-i;
        v(1,n-1) = 0;
        v(M,n-1) = s_max - K*exp(-r*(T-t(n-1)));
        for j = 2:M-1
            v(j,n-1) = v(j,n) + r*s(j)*(0.5*dt/ds)*(v(j+1,n)-v(j-1,n)) + ...
                0.5*sig^2*s(j)^2*(dt/ds^2)*(v(j+1,n) - 2*v(j,n) + v(j-1,n)) - dt*r*v(j,n);
        end
    end
end
%% Implicit
function V_imp = imp_euler(N_imp,V_imp,B_imp,s_max,K,T,t,r)
    for i = 0:N_imp-2
        n_imp = N_imp-i;
        V_imp(:,n_imp-1) = B_imp*V_imp(:,n_imp);
        V_imp(1,n_imp-1) = 0;
        V_imp(end,n_imp-1) = s_max - K*exp(-r*(T-t(n_imp-1))); 
    end
end
%% Matrix Form
function [B, B_imp] = matrix_assembly(M,A,r,s,ds,dt,sig)
    for j = 2:M-1
        A(j,j-1) = -0.5*r*s(j)/ds + 0.5*sig^2*s(j)^2/ds^2;
        A(j,j) = -sig^2*s(j)^2/ds^2 - r;
        A(j,j+1) = 0.5*r*s(j)/ds + 0.5*sig^2*s(j)^2/ds^2;
    end
    B_imp = inv(eye(M) - dt*A);
    B = (eye(M) + dt*A);
end 

function [B_gamma,B_imp_gamma] = matrix_assembly_gamma(M,A_gamma,r,s,ds,dt,sig,gamma)
    for j = 2:M-1
        A_gamma(j,j-1) = -0.5*r*s(j)/ds + 0.5*sig^2*s(j)^(2*gamma)/ds^2;
        A_gamma(j,j) = -sig^2*s(j)^(2*gamma)/ds^2 - r;
        A_gamma(j,j+1) = 0.5*r*s(j)/ds + 0.5*sig^2*s(j)^(2*gamma)/ds^2;
    end
    B_imp_gamma = inv(eye(M,M) - dt*A_gamma);
    B_gamma = (eye(M,M) + dt*A_gamma);
end