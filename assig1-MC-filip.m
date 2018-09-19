% Computational Finance
% Monte Carlo simulation of European call option

% Group 4
clc 
clear

%% Parameter Values
s = 14;
K = 15; 
r = 0.1;
sig = 0.25;
T = 0.5;
gamma = 1;%0.8;

%% Exact - bsexact 
bsexa = bsexact(sig, r, K, T, s);

%% Antithetic variables
samples = 10:100:1010;              % Varying no. samples
iterations = 1000;                  % Fix itereations
no_sims = 20;                       % No. simulations to take average over
dt = (T/iterations);
E = zeros(length(samples),2);
for p = 1:no_sims
    for i = 1:length(samples)
        Niter = samples(i);
        for k=1:Niter               % Sample trajectories
            S0 = s;
            Sp0 = s;
            Sn0 = s;
            for t=1:iterations      % Ito time discretisation
                q = randn;
                Sp1 = Sp0 + r*Sp0*dt + sig*(Sp0^gamma)*sqrt(dt)*q;
                Sn1 = Sn0 + r*Sn0*dt + sig*(Sn0^gamma)*sqrt(dt)*(-q);
                S1 = S0 + r*S0*dt + sig*(S0^gamma)*sqrt(dt)*q;

                Sp0 = Sp1;
                Sn0 = Sn1;
                S0 = S1;
            end
            SpT(k) = Sp1;
            SnT(k) = Sn1;
            ST(k) = S1;

            VpT = max(Sp1-K,0);
            VnT = max(Sn1-K,0);
            VaT(k) = (VpT+VnT)/2;
            VT(k) = max(S1-K,0);
        end
        E_Va = mean(VaT);
        E_V = mean(VT);
        Err(i,1) = abs(bsexa - E_V);
        Err(i,2) = abs(bsexa - E_Va);
    end
    E(:,1) = E(:,1) + Err(:,1);
    E(:,2) = E(:,2) + Err(:,2);
end
avgErr = E./no_sims;                  % Average error
figure
plot(samples,avgErr(:,1),samples,avgErr(:,2))
legend('std. var','antithetic var')
title('Error of "normal"- vs antithetic variables over 20 simulations')
xlabel('samples')
ylabel('Average Error')

%% Plotting the antithetic vs. normal simulation
% The average error when using antithetic variables tends to be lower than
% when not using them. This follows the theory as this reduces the variance
% of the pricing function by taking the average of two values with
% covariance -1.

%% Sample error as function of sample plaths
samp = [1 10 20 30 50 70 100];   % Vary samples
iterations = 100000;                % Fix iterations
dt = (T/iterations);
for i = 1:length(samp)
    for k = 1:samp(i)
        S0 = s;
        for t=1:iterations          % Ito time discretisation 
            S1 = S0 + r*S0*dt + sig*(S0^gamma)*sqrt(dt)*randn;
            S0 = S1;
        end
        ST(k) = S1;
        VT(k) = max(S1-K,0);        % Value function
    end
    V = exp(-r*T)*mean(VT);         % Value of option
    err(i) = abs(bsexa-V);          % Error compared to analytical sol.
end
k = polyfit(log(samp),log(err),1);
ff = samp.^k(1)*exp(k(2));
figure
loglog(samp,err,samp,ff)
title('Sample error as function of sample plaths')
xlabel('Sample paths')
ylabel('Sample Error')

%% Plotting Sampling error
% The Error of the sampleing error as a function of the number of sample
% paths is rather unclear. The error should decrease with the number of
% sample paths. This pattern is visible in some simulations but not others. 

%% Discretization error as a function of the time step. 
sample = 1000000;                    % Fix samples
iterations = 1:1:15;                  % varying no. iterations
for i = 1:length(iterations)
    dt = (T/iterations(i));
    for k = 1:sample
        S0 = s;
        for t=1:iterations(i)         % Ito time discretisation
            S1 = S0 + r*S0*dt + sig*(S0^gamma)*sqrt(dt)*randn;
            S0 = S1;
        end
        ST(k) = S1;
        VT(k) = max(S1-K,0);          % Value function 
    end
    E_V = mean(VT);
    V = exp(-r*T)*E_V;                % Val. of option
    err(i) = abs(bsexa-V);            % Error compared to analytical sol.
end
k = polyfit(log(T./iterations),log(err),1);
gg = (T./iterations).^k(1)*exp(k(2));
figure
loglog(T./iterations,err,T./iterations,gg)
title('Discretization error as a function of the time step')
xlabel('time step')
ylabel('Disc. Error')

%% Plotting Disc. error
% The error of the discretization as a function of the step size decreases
% with a decrease in step size. The slope of the log-log-plot is 1. This
% correspons to a decrease of order 'step size'.  

%% V as function of gamma
iterations = 1000;       % Time steps
samp = 2000;
dt = (T/iterations);
gam = 0.5:0.01:1;   % Gamma vector [0.5 , 1]
V_gam = zeros(length(gam),1);
i = 1;
for gamma = gam
    for k=1:samp                   % Number of trajectories
        Sp0 = s;
        Sn0 = s;
        for t=1:iterations               % Time Steps
            q = randn;
            Sp1 = Sp0 + r*Sp0*dt + sig*(Sp0^gamma)*sqrt(dt)*q;
            Sn1 = Sn0 + r*Sn0*dt + sig*(Sn0^gamma)*sqrt(dt)*-q;

            Sp0 = Sp1;
            Sn0 = Sn1;
        end
        SpT(k) = Sp1;
        SnT(k) = Sn1;

        VpT = max(Sp1-K,0);
        VnT = max(Sn1-K,0);
        VaT(k) = (VpT+VnT)/2;
    end
    E_Va = mean(VaT);
    V_gam(i) = exp(-r*T)*E_Va;        % Value of option
    i = i + 1;
end
% Plotting
figure
plot(gam,V_gam)
title('V as a function of gamma')
xlabel('Gamma')
ylabel('V')

%% Plotting V as function of gamma
% The price of a European Call Option as a function of gamma decreases with
% a decrease of gamma. With the decrease of gamma the volatility decreases.
% This results in a more stable stock and therefore a lower proce of the
% call option as the predicted stock price and the strike price are very
% close. 