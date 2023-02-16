%% CODE FOR HOME-EXERCISE 3. 
% evaluating integrals using Gauss-Legendre quadrature 

% defining integrands to be solved. 

N = 4:4:24; 

error_legendre = zeros(1, length(N)); 
error_simps = zeros(1, length(N)); 
error_bode = zeros(1, length(N)); 

for i = 1:length(N)
    
    error_legendre(i) = abs((pi/2) - Gauss_Legendre(@integrand,N(i))); 
    error_simps(i) = abs((pi/2) - simps(@integrand,-1, 1, N(i))); 
    error_bode(i) = abs((pi/2) - bode(@integrand,-1, 1, N(i))); 
    
end 

hold on 

plot(N, error_legendre, 'rs-')
plot(N, error_simps, 'bs-')

hold off 


% defining the integrand function 
function y = integrand(x)

y = sqrt(1 - x.^2); 

end 

%% IMPLEMENTING "NORMAL" INTEGRAL SOLVERS SIMPS AND BODE 

function integ = simps(func, a, b, N)

% a = lower end, b = upper end, N = number of integration steps

h = (b-a)/N; 

integral = 0; %prealocation 

for n = 1:2:N-1
    
    integral = integral+(func(a+(n-1)*h)+4*func(a+n*h)+func(a+(n+1)*h));
    
end 

integ = (h/3)*integral; 

end 

function integ = bode(f, a, b, N)

% assuming N is a multiple of 4

h = (b-a)/N; % step-length

integral = 0; % prealocation 

for n = 1:4:N-3
    
    %the step is so fkn huge that the term has to be split in two  
    
    part_1 = 7*f(a+(n-1)*h) + 32*f(a+(h*n)); 
    part_2 = 12*f(a+h*(n+1)) + 32*f(a+h*(n+2)) + 7*f(a+h*(n+3)); 
    integral = integral + (((2*h)/45) * (part_1 + part_2)); 
 
end

integ = integral; %returning the final value 

end 
%% 
% implementing a complete Gauss-Legendre integral solver using locally
% defined functions ASSUMES INTEGRAL FROM -1 TO 1 
function integ = Gauss_Legendre(func, l)

tol = 0.0001; % tolerance for root-finding 

x_values = legroots(l, tol); % finding the roots of the legendre function 

h_values = zeros(1, l); 

for i = 1:l
    
    h_values(i) = weights(x_values(i), l); 
    
end 

integ = h_values*func(x_values)'; % numerical integration 

end  

%% implementation of necessary functions 

% a function for finding function weights 
function w = weights(root, l)

w = 2/((1-root^2)*(legder(root, l))^2);

end 

% a function for finding all the roots of a legendre polynomial of order l 
function roots = legroots(l, tol)

% l - order of legendre polynomial 
% tol - tolerance of sltn 

% this is the point of failure of the method. increase step to increase the
% order of acceptable l, but this increases process time. 
step = (1-(-1))/(20*l); % a very rough step-value 

a = -1; % setting the first a - value 
b = a; 

roots = zeros(1,l); % each legendre-polynomial has l roots 

for i = 1:l
    
    % while-loop for finding appropriate step-values 
    while sign(legendre(a,l)) == sign(legendre(b,l)) 
        
        b = b + step; 
    end 
    
    roots(i) = bisec(a, b, l, @legendre, tol); 
    
    a = b; 
end 

end 

% implementation of bisection method for finding roots of a function f 
function root = bisec(a, b, l, f, tol)

% a - lower bound 
% b - upper bound 
% l - degree of legendre polynomial 
% f - function handle under consideration 

if sign(f(a, l)) == sign(f(b, l))
    error('BOTH POINTS HAVE THE SAME SIGN')
end 

while abs(b - a) > tol 
    c = a + ((b-a)/2); % finding the midpoint 
    
    if sign(f(a, l)) == sign(f(c, l))
        a = c; 
    else
        b = c; 
    end 
end 

root = a + ((b-a)/2); 
    
end 

% function for finding the derivative of legendre polynomial l at point x
function y = legder(x,l)

% given formula for determining legendre polynomial derivative
y = (-l*x*legendre(x,l)+l*legendre(x,l-1))/(1-x^2);

end 

% recursive function for finding the value of legendre polynomial l at point x
function y = legendre(x, l)


if l == 0 
    y = 1;
    
elseif l == 1
    y = x; 
    
else 
    y = ((2*l-1)*x*legendre(x,l-1)-(l-1)*legendre(x,l-2))/l;
   
end 

end 