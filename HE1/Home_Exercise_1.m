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
N = 100; % number of points 
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

% plotting the sltns
hold on 
plot(t, analytic(c1, c2, t)) % analytic, for comp. 
%plot(t, analyticp(c1, c2, t)) % analytic, for comp.

plot(t, displacement, 'rs') % numerical solution INSANE! 

hold off
