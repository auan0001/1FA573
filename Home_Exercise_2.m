% Code for home exercise number 2. 
% problem focuses on implementation of numeric integration for solution of
% Green's equation 

% Defining parameter values 
r_max = 30; 
a = 4; % given constant 
N = 1000; % number of steps 
r_values = linspace(0, r_max/3, 1000); 
y_numeric = zeros(1, length(r_values)); % prealocation

% numeric solutions. This is not efficient, but does not need to be. 
for i = 1:length(r_values)
    
    y_numeric(i) = numeric_phi(r_values(i), a, r_max); 
    
end 

% plotting the solution 
hold on

plot(r_values, analytic(a,r_values))  
plot(r_values, y_numeric)
legend('analytic solution','numeric solution')

hold off

%% DEFINITION OF LOCAL FUNCTIONS 

% defining a numeric solution to the entire problem 

function y = numeric_phi(r, a, r_max)

int_steps  = 1000; % number of integration steps 

% we split the sltn in two for readability 
first = phi_great(a, r)*bode(@less_integ, 0, r, int_steps, a);  
second = phi_less(a, r)*bode(@great_integ, r, r_max, int_steps, a); 

y = first + second;  

end 

% numerical integration using Bode's method 
% Tests can be found in lecture_1_2.m

function integ = bode(f, a, b, N, a_val)

% f - function handle, to be integrated 
% a , b - upper and lower bounds 
% N, number of steps 
% a_val - given a-constant 

% assuming N is a multiple of 4

h = (b-a)/N; % step-length

integral = 0; % prealocation 

for n = 1:4:N-3 % CHECK THIS
  
    part_1 = 7*f(a_val, a+(n-1)*h) + 32*f(a_val, a+(h*n)); 
    part_2 = 12*f(a_val, a+h*(n+1)) + 32*f(a_val, a+h*(n+2)) + 7*f(a_val, a+h*(n+3)); 
    integral = integral + (((2*h)/45) * (part_1 + part_2)); 
 
end

integ = integral; %returning the final value 

end 


% defining integrand for phi_less term 
function y = less_integ(a, r)

source = -((r*exp(-r))/2); 
y = phi_less(a, r) * source; 

end 

% defining integrand for phi_great term 
function y = great_integ(a, r)

source = -((r*exp(-r))/2); 
y = phi_great(a, r) * source; 

end 


% first slt one to the homogeneous differential eqn, given 
function y = phi_less(a, r)

    y = (1/sqrt(2*a))*(exp(a*r) - exp(-a*r)); 

end 

% second slt one to the homogeneous differential eqn, given 
function y = phi_great(a, r)

    y = -(1/sqrt(2*a))*exp(-a*r); 

end 


% defining the given analytic sltn
function y = analytic(a, r)

    y = ((1/(1-a^2))^2)*(exp(-a*r) - exp(-r).*(1+0.5*((1-a^2).*r))); 

end 