% FINAL PROJECT 2. 

% Code for final project, part 2. 

% definition of variables 
V0 = 10; 
a = 1; 
E = 0.1 *V0; % should be taken in intervals of 0.1V0 and 100V0 
r_max = 3*a; % R-max should be set to 3*a
b_values = linspace(0, r_max, 300); 
numerical = zeros(length(b_values), 1); 

for i = 1:length(numerical)
    
    % at each point in the iteration the value of b is changed 
    numerical(i) = num_theta(b_values(i), E, r_max); 
    
end 

% plotting the analytic vs the numeric sltns
hold on 
plot(b_values, numerical, 'rs-')
grid on
xlabel('Impact parameter')
ylabel('Deflection angle in radians')
hold off

%% DEFINITION OF FUNCTIONS 

% numeric deflection function
function t = num_theta(b, E, r_max)

r_min = bisec(b, r_max, b, E, @r_min);


first = bode(@integrand_1, 0, sqrt(r_max-b), b, E, b, b); 
second = bode(@integrand_2, 0, sqrt(r_max-r_min), b, E, r_min);  

t = 2*b*(first - second);

end 

% transformed first integrand in expression I8
function y = integrand_1(p, b, E, r_min)

if p == 0 
    y = 0; 
else 
    y = 2*p/((p^2 + b)^2 * sqrt(1 - b^2/(p^2+b)^2)); 
end 

end 

% transformed second integrand in expression I8
function y = integrand_2(p, b, E, r_min)

if p == 0 
    y = 0; 
else 
    y = 2*p/((p^2 + r_min)^2 * sqrt(1 - b^2/(p^2+r_min)^2 - V(p^2 + r_min)/E)); 
end 

end

% numeric integrand solver
function integ = bode(f, l, u, b, E, r_min, bnd)

% f - function handle under consideration 
% l - lower bound 
% u - upper bound 
% b - impact parameter 
% E - total energy of particle 
% r_min - minimum r-value 

N = 10000;  % total number of integration steps ;
h = (u-l)/N; % step-length

integral = 0; % prealocation 

for n = 1:4:N-3
    
    %the step is so fkn huge that the term has to be split in three 
    part_1 = 7*f(l+(n-1)*h, b, E, r_min) + 32*f(l+(h*n), b, E, r_min); 
    part_2 = 12*f(l+h*(n+1), b, E, r_min) + 32*f(l+h*(n+2), b, E, r_min);  
    part_3 = 7*f(l+h*(n+3), b, E, r_min);
    
    integral = integral + (((2*h)/45) * (part_1 + part_2 + part_3)); 
 
end

integ = integral; %returning the final value 

end 

% implementing bisec. method for finding the minimum of r_min
function root = bisec(l, u, b, E, f)

% l - lower bound 
% u - upper bound 
% b - impact parameter 
% V - potential 
% E - total energy, constant 
% f - function handle under consideration 

tol = 0.001; % tolerance of calculation

while abs(l - u) > tol 
    c = l + ((u-l)/2); % finding the midpoint 
    
    if sign(f(l, b, E)) == sign(f(c, b, E))
        l = c; 
    else
        u = c; 
    end 
end 

root = l + ((u-l)/2);
    
end 

% definition of the r_min-function
function rad = r_min(r, b, E)

rad = 1 - b^2/r^2 - V(r)/E;

end 

% definition of the potential-function
function pot = V(r)

% magic constants. choose at your peril 
V0 = 10; 
a = 1; 

pot = 4*V0 * ((a/r)^12 - (a/r)^6); 

end 