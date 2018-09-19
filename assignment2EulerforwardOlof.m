clc; clear; close all;

%% Parameters and vectors

K = 15;
r = 0.1;
sigma = 0.25;
T = 0.5;
smax = 4*K;
gamma = 1;
timesteps = [50000]; % Number of steps in time
spacesteps = [10, 1000]; % Number of steps in space


%% Boundary & final conditions, simulations and vectors
for l = 1:length(timesteps)
    t = linspace(0,T,timesteps(l));
    s = linspace(0,smax,spacesteps(end));
    V = zeros(spacesteps(end),timesteps(l)); %V(space,time)
    Vexact = zeros(spacesteps(end),1);
    for i = 1:spacesteps(end)
        Vexact(i) = bsexact(sigma,r,K,T,s(i));
    end

    %Boundary condition
    V(1,:) = 0;
    V(spacesteps(end),:) = smax - K*exp(-r*(T-t));

    %Final condition
    V(:,timesteps(l)) = max(s-K,0);

    %Euler forward using loops, Explicit method
    for n = timesteps(l):-1:2
        deltat = t(n)-t(n-1); 
        for j = 2:(spacesteps(end) - 1)
            deltas = s(j) - s(j-1);
            V(j,n-1) = V(j,n) + r*s(j)*deltat/(2*deltas)*(V(j+1,n)-V(j-1,n)) ...
            + (sigma^2/2)*(s(j)^(2*gamma))*deltat/(deltas^2)*(V(j+1,n) - 2*V(j,n) + V(j-1,n)) ...
            - deltat*r*V(j,n);
        end
    end
    
    errorMax(l) = max(abs(V(:,1)-Vexact));
end
%% Plots
%Comparing Numerical and analytical
figure()
plot(s,V(:,1),s,Vexact)
xlim([55 60])
legend('Numerical','Analytical')
xlabel('Price')
ylabel('Value')
%Error propagation plot
figure()
plot(s,abs(V(:,1)-Vexact))
xlabel('price S')
ylabel('Absolute error')
%size of error
figure()
loglog(timesteps,errorMax)
hold on
scatter(timesteps,errorMax)
grid on
xlabel('Number of time discretizations')
ylabel('Maximum absolute error')

