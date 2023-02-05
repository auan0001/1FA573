%% code for Home exercise no. 1. 
clear
% add functions from subdir
addpath ./functions
set(0,'defaulttextinterpreter','latex')
% analysis of two coupled 1st order differential equations 
 
%% code for Home exercise no. 1. 
% analysis of two coupled 1st order differential equations 
 
% setting initial conditions
y_0 = 1; % start-displacement 
p_0 = 1; % end-displacement 
a = 0; % start-value in time 
b = 1; % end-value in time 
N = 1000; % number of points 
h = (b-a)/N; % step-size in time 
t = linspace(a, b, N); % time-values to be controlled 

% constants for the analytic sltn
c1 = y_0; 
c2 = p_0/(2*pi); 

% numerically determining the solution values 

displacement = zeros(1, N); 
momentum = zeros(1, N); 

displacement(1) = y_0; 
momentum(1) = p_0; 

for i = 2:N
    % general idea: each solution uses the other's solution
    displacement(i) = RK4(@y_dt,momentum(i-1), displacement(i-1), t(i), h); 
    momentum(i) = RK4(@p_dt, displacement(i-1), momentum(i-1), t(i), h); 
end 


%% PLOTTING THE SOLUTIONS 
% Gap size between plotted data points
gap = N/50;

% Displacement
figure
subplot(2,3,1)
plot(t, analytic(c1, c2, t)) % analytic, for comp. 
hold on 
grid on
plot(t(1:gap:end), displacement(1:gap:end), 'rs') % numerical solution INSANE! 
title('Displacement')
leg1 = legend('$y$', '$y_n$', 'location', 'southwest');
set(leg1, 'Interpreter','latex');
xlabel('$t \in [0,1]$')
ylabel('$y$')

% Momentum
subplot(2,3,4)
plot(t, analyticp(c1, c2, t)) % analytic, for comp.
hold on
grid on
plot(t(1:gap:end), momentum(1:gap:end), 'rs') % numerical solution INSANE!
title('Momentum')
leg2 = legend('$p$', '$p_n$', 'location', 'northwest');
set(leg2, 'Interpreter','latex');
xlabel('$t \in [0,1]$')
ylabel('$p$')

% Phase plane
subplot(2,3,[2,3,5,6])
plot(analytic(c1, c2, t), analyticp(c1, c2, t)) % analytic, for comp.
hold on
grid on
plot(displacement(1:gap:end), momentum(1:gap:end), 'rs') % numerical solution INSANE!
title('Phase plane')
leg2 = legend('$(y,p)$', '$(y_n,p_n)$', 'location', 'northwest');
set(leg2, 'Interpreter','latex');
ylabel('Displacement $y$')
xlabel('$y$')
ylabel('$p$')

sgtitle(strcat('Analytic and RK4 numerical solutions for $N=', num2str(N), '$'))
hold off
