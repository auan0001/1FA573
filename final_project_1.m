% Code for final project, part 1. 

r_max = 1000; % completely arb. value 
r = linspace(0, r_max);
numerical = zeros(length(r), 1); 

V = 2; % PLACEHOLDER - CHECK!! 
E = 1; % PLACEHOLDER - CHECK!! 

for i = 1:length(numerical)
    
    % at each point in the iteration the value of b is changed 
    numerical(i) = num_theta(r(i), V , E, r_max); 
    
end 

% plotting the analytic vs the numeric sltns
% ANSWER IS INSANE 
hold on 

plot(r, analytic_1(r))
plot(r, numerical, 'rs-')

hold off

%% Functions for determining the numeric sltn

% numeric deflection function

function t = num_theta(b, V , E, r_max)

tol = 0.0001; 
N = 1000; 
r_min = bisec(0, r_max, b, V, E, @r_min, tol);

first = bode(@integrand_1, b, r_max, b, E, V, N); 
second = bode(@integrand_2, r_min, r_max, b, E, V, N); 

t = (2*b) * (first + second);

end 

% numeric integrand solver
function integ = bode(f, l, u, b, E, V, N)

% l - lower bound 
% u - upper bound 
% f - function handle under consideration 
% N - number of integration steps, N must be a multiple of 4

h = (u-l)/N; % step-length

integral = 0; % prealocation 

for n = 1:4:N-3
    
    %the step is so fkn huge that the term has to be split in two  
    
    part_1 = 7*f(l+(n-1)*h, b, V , E) + 32*f(l+(h*n), b, V , E); 
    part_2 = 12*f(l+h*(n+1), b, V , E) + 32*f(l+h*(n+2), b, V , E) + 7*f(l+h*(n+3), b, V , E); 
    integral = integral + (((2*h)/45) * (part_1 + part_2)); 
 
end

integ = integral; %returning the final value 

end 

% first integrand in expression I8
% unused variables maintained for generality 
function y = integrand_1(r, b, V , E)

y = (1/(r^.2)) * sqrt(1 - ((b^2)/(r^2))); 

end 

% second integrand in expression I8
function y = integrand_2(r, b, V , E)

y = (1/(r^.2)) * sqrt(1 - ((b^2)/(r^2)) - (V/E)); 

end 

% implementing bisec. method for finding the minimal of r_min
% implementation of bisection method for finding roots of a function f 
function root = bisec(l, u, b, V, E, f, tol)

% l - lower bound 
% u - upper bound 
% f - function handle under consideration 

while abs(l - u) > tol 
    c = l + ((u-l)/2); % finding the midpoint 
    
    if sign(f(l, b, V, E)) == sign(f(c, b, V, E))
        l = c; 
    else
        u = c; 
    end 
end 

root = l + ((u-l)/2); 
    
end 

% definition of the r_min-function
function r = r_min(r, b, V, E)

r = 1 - ((b^2)/(r^2)) - (V/E);

end 

%% Analytic solutions - for reference 

function theta = analytic_1(b)

r_max = 1000; % pre-given value of r_max 

theta = pi - asin(b./r_max); 

end 

function theta = analytic_2(b, E, V0)

theta = 2 * (asin(b/(r_max*sqty(1 - (V0/E)))) - asin(b/r_max));

end 