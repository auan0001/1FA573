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
%% DEFINING NUMERIC ODE SOLVERS 

% RK2-function which returns the function value at the next time-incident 
function res = RK2(func, y_calc, y_cur, t, h)

    k = h*func(t, y_calc); 
    res = y_cur + h*func(t + 0.5*h, y_calc + 0.5*k); 

end 


% RK2-function which returns the function value at the next time-incident 
function res = RK4(func, y_calc, y_cur, t, h)

    k1 = h*func(t, y_calc); 
    k2 = h*func(t+h*(0.5), y_calc + k1/2);
    k3 = h*func(t+h*(0.5), y_calc + k2/2);
    k4 = h*func(t+h, y_calc + k3); 

    res = y_cur + (1/6)*(k1 + 2*k2 + + 2*k3 + k4);

end 

%% 
%defining our diferential equations and analytic solution 

% ODE for displacement 
function y_prim = y_dt(t, p)
    
    y_prim = p; 
     
end 

% ODE for momentum 
function p_prim = p_dt(t, y)

    p_prim = -(4*pi^2)*y;  

end 

% analytic solution for displacement 
function ana = analytic(c1, c2, t)

   ana = c1*cos(2*pi*t) + c2*sin(2*pi*t); 

end 

% analytic solution for momentum 
function anap = analyticp(c1, c2, t)

   anap = c1*2*pi*cos(2*pi*t) - c2*2*pi*sin(2*pi*t); 

end 