%% code for Home exercise no. 1. 
clear
% add functions from subdir
addpath ./functions
set(0,'defaulttextinterpreter','latex')
% analysis of two coupled 1st order differential equations 
 
% setting initial conditions
y_0 = 1; % start-displacement 
p_0 = 1; % end-displacement 
a = 0; % start-value in time 
b = 1; % end-value in time 

% constants for the analytic sltn
c1 = y_0; 
c2 = p_0/(2*pi); 

% numerically determining the solution values 

% logspace error sweep
N = logspace(2,5,4);

figure
for steps = 1:length(N)

  grid = N(steps);
  h(steps) = (b-a)/grid; % step-size in time 
  t = linspace(a, b, grid); % time-values to be controlled 
  displacement = zeros(1, grid); 
  momentum = zeros(1, grid); 

  displacement(1) = y_0; 
  momentum(1) = p_0; 

  % hold on
  for i = 2:grid

    % general idea: each solution uses the other's solution
    displacement(i) = RK4(@y_dt,momentum(i-1), displacement(i-1), t(i), h(steps)); 
    momentum(i) = RK4(@p_dt, displacement(i-1), momentum(i-1), t(i), h(steps));
  end 
  % Absolute errors
  err = abs(displacement-analytic(c1,c2,t));
  err_p = abs(momentum-analyticp(c1,c2,t));

  % Plot either mom or displ.

  % displ plot
  % hold on
  % grid on
  % plot(t, err, '-.') % numerical solution INSANE! 
  % title('Absolute displacement error in the time domain')
  % xlabel('$t \in [0,1]$')
  % ylabel('$\tilde{e}_y = |y_n - \tilde{y}_n|$')
  % end displ plot

  % mom plot
  hold on
  grid on
  plot(t,err_p, '-.') % numerical solution INSANE! 
  title('Absolute momentum error in the time domain')
  xlabel('$t \in [0,1]$')
  ylabel('$\tilde{e}_p = |p_n - \tilde{p}_n|$')
  % end mom plot

  % Total distance error
  total_dist_err(steps) = err*err';  
  total_dist_err_p(steps) = err_p*err_p';  
end
% legend
hl = legend(strcat('$h=', string(num2cell(h))+'$'), 'location', 'northwest')
set(hl, 'Interpreter','latex')
hold off

%% TOTAL DISTANCE ERROR

figure
loglog(total_dist_err, (b-a)./N, '--v')
hold on
loglog(total_dist_err_p, (b-a)./N, '--v')
grid on
title('Total distance error vs step size')
xlabel('$h$')
ylabel('$\tilde{e} = e^T_n e_n$')
leg = legend({'$\tilde{e}_y$', '$\tilde{e}_p$'}, 'location', 'northwest')
set(leg, 'Interpreter','latex')
hold off
