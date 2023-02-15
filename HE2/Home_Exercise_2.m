% Code for home exercise number 2. 
% problem focuses on implementation of numeric integration for solution of
% Green's equation 
addpath ./functions
set(0, 'defaulttextinterpreter', 'latex')
hAxes.TickLabelInterpreter = 'latex';

% Defining parameter values 
r_max = 30; 
a = 4; % given constant 
N = 1000; % number of steps 

N = logspace(2,4,3)

for i = 1:length(N)

  r_values = linspace(0, r_max, N(i)); 
  y_numeric = zeros(1, length(r_values)); % prealocation

  % numeric solutions. This is not efficient, but does not need to be. 
  for j = 1:length(r_values)
    y_numeric(j) = numeric_phi(r_values(j), a, r_max, N(i)); 
  end 
end

% plotting the solution 

gap = N/50;

abs_error = abs(analytic(a,r_values) - y_numeric);

figure;
plot(r_values, abs_error)
grid on
hold off
title('Absolute error')
xlabel('$r$')
ylabel('$|\phi-\phi_n|$')

figure;
plot(r_values, analytic(a,r_values))  
hold on
grid on
plot(r_values(1:gap:end), y_numeric(1:gap:end), 'rs')
leg = legend({'Analytical $\phi(r)$', 'Numerical $\phi_n(r)$'}, 'location', 'northeast');
set(leg, 'Interpreter','latex')
title('Numerical solution')
xlabel('$r$')
ylabel('$\phi$')
hold off
