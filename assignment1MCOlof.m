clc; clear; close all;

%% Parameters & vectors
S0 = 14;
K = 15;
r = 0.1;
sigma = 0.25;
T = 0.5;
gamma = 1;
M = [10, 100, 1000, 10000, 100000];
timestep = [0.01, 0.001, 0.0001];
sizeofM = length(M);
sizeoft = length(timestep);
t0 = 0;
Vexact = bsexact(sigma,r,K,T,S0);

errorM = zeros(1,sizeofM); %sample error
errorManti = zeros(1,sizeofM);
errorT = zeros(1,sizeoft); %time error

%% Monte carlo and Euler method

%Varying number of simulations
for j = 1:sizeofM
    V = zeros(1,M(j));
    Vanti = zeros(1,M(j));
    deltat = timestep(sizeoft);
    for i = 1:M(j)
        %Euler method
        S = S0;
        Splus = S0;
        Sminus = S0;
        for t = t0:deltat:T
           S = S + r*S*deltat + sigma*(S^gamma)*sqrt(deltat)*randn;
           Splus = Splus + r*Splus*deltat + sigma*(Splus^gamma)*sqrt(deltat)*randn;
           Sminus = Sminus + r*Sminus*deltat - sigma*(Sminus^gamma)*sqrt(deltat)*randn;
        end
        V(i) = max(S-K,0);
        Vanti(i) = (max(Splus-K,0) + max(Sminus-K,0))/2;
    end

    EQ = mean(V);
    V0 = exp(-r*T)*EQ;
    V0anti = exp(-r*T)*mean(Vanti);

    errorM(j) = abs(V0 - Vexact);
    errorManti(j) = abs(V0anti - Vexact);
end

%Varying number of time steps 
for j = 1:sizeoft
    V = zeros(1,M(sizeofM));
    deltat = timestep(j);
    
    for i = 1:M(sizeofM)
        %Euler method
        S = S0;
        for t = t0:deltat:T
           S = S + r*S*deltat + sigma*(S^gamma)*sqrt(deltat)*randn;
        end
        V(i) = max(S-K,0);
    end

    EQ = mean(V);
    V0 = exp(-r*T)*EQ;

    errorT(j) = abs(V0 - Vexact);
end

%% Plots
figure()
plot(M,errorM,M,errorManti)
loglog(M,errorM)
hold on
scatter(M,errorM)
loglog(M,errorManti)
scatter(M,errorManti)
legend('errorM','errorManti')
title(['Sample error, dt = ' num2str(timestep(sizeoft))])
figure()
loglog(timestep,errorT)
title(['discretization error, number of simulations = ' num2str(M(sizeofM))])

