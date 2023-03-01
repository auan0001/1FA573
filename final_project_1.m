% Code for final project, part 1. 

r_max = 1000; % completely arb. value 
r = linspace(0, r_max, 300); 
numerical = zeros(length(r), 1); 

V = 1; % PLACEHOLDER - CHECK!! 
E =-1; % PLACEHOLDER - CHECK!! 

for i = 1:length(numerical)
    
    % at each point in the iteration the value of b is changed 
    numerical(i) = num_theta_1(r(i), V , E, r_max); 
    
end 

% plotting the analytic vs the numeric sltns
hold on 

plot(r, analytic_1(r, E, V, r_max))
plot(r, numerical, 'rs-')

hold off

%% Functions for determining the numeric sltn

% numeric deflection function
% VALID FOR E<U0
function t = num_theta_1(b, V , E, r_max)

tol = 0.0001; % tolerance when finding r_min value 
N = 1000;  % total number of integration steps ;
r_min = bisec(0, r_max, b, V, E, @r_min, tol);

lower1 = 0.01; % we take a value close to 0, but not 0. 
upper1 = sqrt(r_max - b);

lower2 = sqrt(r_max - r_min); 
upper2 = sqrt(r_max - b);

first = bode(@integrand_3, lower1, upper1, b, E, V, N, r_min);
second = bode(@integrand_4, lower2, upper2, b, E, V, N, r_min);

t = first - second;

end 

% numeric deflection function
% VALID FOR E>U0
% DOES NOT WORK FOR SOME REASON
function t = num_theta_2(b, V , E, r_max)

tol = 0.0001; % tolerance when finding r_min value 
N = 1000;  % total number of integration steps ;
r_min = bisec(0, r_max, b, V, E, @r_min, tol);

lower1 = sqrt(b-r_min); % we take a value close to 0, but not 0. 
upper1 = sqrt(r_max - r_min);

lower2 = 0.01; 
upper2 = sqrt(r_max - r_min);

first = bode(@integrand_5, lower1, upper1, b, E, V, N, r_min);
second = bode(@integrand_6, lower2, upper2, b, E, V, N, r_min);

t = first - second;

end 


% numeric integrand solver
function integ = bode(f, l, u, b, E, V, N, r_min)

% l - lower bound 
% u - upper bound 
% f - function handle under consideration 
% N - number of integration steps, N must be a multiple of 4

h = (u-l)/N; % step-length

integral = 0; % prealocation 

for n = 1:4:N-3
    
    %the step is so fkn huge that the term has to be split in two  
    
    part_1 = 7*f(l+(n-1)*h, b, V , E, r_min) + 32*f(l+(h*n), b, V , E, r_min); 
    part_2 = 12*f(l+h*(n+1), b, V , E, r_min) + 32*f(l+h*(n+2), b, V , E, r_min) + 7*f(l+h*(n+3), b, V , E, r_min); 
    integral = integral + (((2*h)/45) * (part_1 + part_2)); 
 
end

integ = integral; %returning the final value 

end 

%% ORIGINAL INTEGRAND EXPRESSIONS - DO NOT USE 
% first integrand in expression I8
% unused variables maintained for generality 
function y = integrand_1(r, b, V , E, r_min)

y = 1/((r^.2) * sqrt(1 - ((b^2)/(r^2)))); 

end 

% second integrand in expression I8
function y = integrand_2(r, b, V , E, r_min)

y = 1/((r^.2) * sqrt(1 - ((b^2)/(r^2)) - (V/E)));

end

%% TRANSFORMED INTEGRANDS FOR NUMERICAL SLTNS 

% VALID FOR E<U0

% first integrand after conversion 
% unused variables maintained for generality 
function y = integrand_3(p, b, V , E, r_min)

y = (4*p*b)/((p^2 + b)^2 * sqrt(1 - b^2/(((p^2 + b)^2)))); 

end 

% second integrand after conversion
function y = integrand_4(p, b, V , E, r_min)

y = (4*p*b)/((p^2 + b)^2 * sqrt(1 - b^2/(((p^2 + b)^2)) - V/E)); 

end 

%% TRANSFORMED INTEGRANDS FOR NUMERICAL SLTNS 

% VALID FOR E>U0

% first integrand after conversion 
% unused variables maintained for generality 
function y = integrand_5(p, b, V , E, r_min)

y = (4*p*b)/((p^2 + r_min)^2 * sqrt(1 - b^2/((p^2 + r_min)^2))); 

end 

% second integrand after conversion
function y = integrand_6(p, b, V , E, r_min)

y = (4*p*b)/((p^2 + r_min)^2 * sqrt(1 - b^2/((p^2 + r_min)^2) - V/E)); 

end 

%%
% implementing bisec. method for finding the minimum of r_min
function root = bisec(l, u, b, V, E, f, tol)

% l - lower bound 
% u - upper bound 
% b - impact parameter 
% V - potential 
% E - total energy, constant 
% f - function handle under consideration 
% tol - tolerance of calculation

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
% THIS GIVES INFINITE INPUTS 
function rad = r_min(r, b, V, E)

rad = 1 - ((b^2)/(r^2)) - (V/E);

end 

%% Analytic solutions - for reference 

% VALID FOR E>U0
function theta = analytic_1(b, E, V0, r_max)

r_max = 1000; % pre-given value of r_max 

theta = pi - asin(b./r_max); 

end 


% VALID FOR E<U0
function theta = analytic_2(b, E, V0, r_max)

theta = 2 * (asin(b./(r_max*sqrt(1 - (V0/E)))) - asin(b./r_max));

end 