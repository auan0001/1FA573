clear
set(0, 'defaulttextinterpreter', 'latex')
hAxes.TickLabelInterpreter = 'latex';

% Parameters
N_b = 400; % _b for bisection
a = -1;
b = 1;
tol = 1e-5;

% Soln to integral
analytic = pi/2;

% 
N = 4:4:28;

% Error loop
for i = 1:length(N)
  h = (b-a)/N(i);
  h_b = (b-a)/N_b;
  x_b = linspace(a, b, N_b);

  % Legendre polynomial
  coeff = legendre_coeff(N(i));
  P_n = @(x) polyval(coeff, x);

  integral_gl = gauss_legendre(@integrand, P_n, N(i), x_b, N_b, h_b, tol);
  err_gl(i) = abs_err(integral_gl, analytic);

  integral_bode = bode(@integrand, a, b, h);
  err_bode(i) = abs_err(integral_bode, analytic);

  integral_simp = simpson(@integrand, a, b, h);
  err_simp(i) = abs_err(integral_simp, analytic);
end

% Error loop to get mean computational time
for i = 1:length(N)
  for j = 1:length(N)
    h = (b-a)/N(i);
    h_b = (b-a)/N_b;
    x_b = linspace(a, b, N_b);

    % Legendre polynomial
    coeff = legendre_coeff(N(i));
    P_n = @(x) polyval(coeff, x);

    tic;
    integral_gl = gauss_legendre(@integrand, P_n, N(i), x_b, N_b, h_b, tol);
    t_gl(i,j) = toc;

    tic;
    integral_bode = bode(@integrand, a, b, h);
    t_bode(i,j) = toc;

    tic;
    integral_simp = simpson(@integrand, a, b, h);
    t_simp(i,j) = toc;
  end
end

handle(1) = figure;
semilogy(N, err_gl, '.-')
hold on
grid on
semilogy(N, err_bode, '.-')
semilogy(N, err_simp, '.-')
err_leg = legend({'G-L', 'Bode', 'Simpson'}, 'location', 'northeast');
set(err_leg, 'Interpreter','latex')
title('Absolute error')
xlabel('N')
ylabel('$L^1$ Error')
hold off


handle(2) = figure;
% Times from mean of matrix rows
semilogy(N, mean(t_gl, 2), '.-')
hold on
grid on
semilogy(N, mean(t_bode, 2), '.-')
semilogy(N, mean(t_simp, 2), '.-')
t_leg = legend({'G-L', 'Bode', 'Simpson'}, 'location', 'northwest');
set(t_leg, 'Interpreter','latex')
title('Mean computation times')
xlabel('N')
ylabel('Time [s]')
hold off

handle(3) = figure;
% Times from mean of matrix rows
loglog(mean(t_gl, 2), err_gl, '.-')
hold on
grid on
loglog(mean(t_bode, 2), err_bode, '.-')
loglog(mean(t_simp, 2), err_simp, '.-')
t_leg = legend({'G-L', 'Bode', 'Simpson'}, 'location', 'northeast');
set(t_leg, 'Interpreter','latex')
title('Efficiency')
xlabel('Time [s]')
ylabel('$L^1$ Error')
hold off

% Save the figures
savefig(handle, 'HE3.fig')

%% FUNCTIONS

function err = abs_err(numeric, analytic)
  % Absolte error
  err = abs(analytic - numeric);
end

function x_n = bisect_roots(f, x, N, h, tol)
  % Bisection method, for multiple roots, given a function handle

  % Compute alternating signs
  x_sign = sign(f(x));
  x_alt_sign = x_sign(1:N-1).*x_sign(2:N);
  % Find alternating signs
  x_alt_sign_idx = x(x_alt_sign < 0);
  % Bisect the root. Might need to adjust step size
  % depending on the discretization level
  for i = 1:length(x_alt_sign_idx)
    x_n(i) = bisect(x_alt_sign_idx(i), x_alt_sign_idx(i)+2*h, f, tol);
  end
end

% defining the integrand function 
function y = integrand(x)
  y = sqrt(1 - x.^2); 
end 

function I = gauss_legendre(integrand, P_n, l, x_bisect, N_bisect, h, tol)
  % Find roots, generate weights and compute integral
  x_n = bisect_roots(P_n, x_bisect, N_bisect, h, tol);
  I = sum(weights(x_n, l).*integrand(x_n));
end

function I = bode(f,a,b,h)
  x = a:4*h:b-3*h;
  I = 2*h/45*sum((7*f(x) + 32*f(x+h) + 12*f(x+2*h)+ 32*f(x+3*h) + 7*f(x+4*h)));
end

function I = simpson(f,a,b,h)
  x = a:h:b;
  I = h/3*(f(x(1)) + 2*sum(f(x(3:2:end-2))) + 4*sum(f(x(2:2:end))) + f(x(end)));
end

function w = weights(x_n, order)
  % Compute weights
  coeff = legendre_coeff(order+1);
  w = 2*(1-x_n.^2)./((order+1)^2.*(polyval(coeff, x_n).^2));
end 

function root = bisect(a, b, f, tol)
  % Implementation of bisection method for finding roots of a function f 

  if sign(f(a)) == sign(f(b))
    error('BOTH POINTS HAVE THE SAME SIGN')
  end 

  while abs(b - a) > tol 
    c = a + ((b-a)/2); % finding the midpoint 

    if sign(f(a)) == sign(f(c))
      a = c; 
    else
      b = c; 
    end 

  end 

  root = a + ((b-a)/2); 
end 

function coeff = legendre_coeff(l)
  % Compute a vector of coefficients for P_n of order n
  if l == 0
    coeff = 1;
  elseif l == 1
    coeff = [1 0];
  else
    coeff = ((2*l-1)*[legendre_coeff(l-1),0] - (l-1)*[0,0,legendre_coeff(l-2)])/l;
  end
end
