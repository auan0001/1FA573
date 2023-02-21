% Parameters
clear

N_bisect = 400;
N = 24;
a = -1;
b = 1;
tol = 1e-5;
h = (b-a)/N;

h_bisect = (b-a)/N_bisect;
x_bisect = linspace(a,b,N_bisect);

% Legendre polynomial
coeff = legendre_coeff(N);
P_n = @(x) polyval(coeff, x);

% Integrate
tic
integral_gl = gausslegendre(@integrand, P_n, x_bisect, N_bisect, N, h_bisect, tol);
toc
disp("G-L = " + integral_gl)

tic
integral_bode = bode(@integrand, a, b, h);
toc
disp("Bode = " + integral_bode)

tic
integral_simpson = simpson(@integrand, a, b, h);
toc
disp("Simpson = " + integral_simpson)

function x_n = bisect_roots(f, x, N, h, tol)
  % Bisection method, for multiple roots, given a function handle

  % Compute alternating signs
  x_sign = sign(f(x));
  x_alt_sign = x_sign(1:N-1).*x_sign(2:N);
  % Find alternating signs
  x_alt_sign_idx = x(find(x_alt_sign < 0));
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

%% IMPLEMENTING "NORMAL" INTEGRAL SOLVERS SIMPS AND BODE 

function I = gausslegendre(integrand, P_n, x_bisect, N_bisect, l, h, tol)
  % Find roots, generate weights and compute integral
  x_n = bisect_roots(P_n, x_bisect, N_bisect, h, tol);
  I = sum(weights(x_n, l).*integrand(x_n));

end

function I = bode(f,a,b,h)
  x = a:4*h:b-3*h;
  I = sum(2*h/45*(7*f(x) + 32*f(x+h) + 12*f(x+2*h)+ 32*f(x+3*h) + 7*f(x+4*h)));
end

function I = simpson(f,a,b,h)
  x = a:h:b;
  I = h/3*(f(x(1))+2*sum(f(x(3:2:end-2)))+4*sum(f(x(2:2:end)))+f(x(end)));
end

%% implementation of necessary functions 

function w = weights(x_n, order)
  % Compute weights
  coeff = legendre_coeff(order+1);
  w = 2*(1-x_n.^2)./((order+1)^2.*(polyval(coeff, x_n).^2));

end 


function root = bisect(a, b, f, tol)
  % Implementation of bisection method for finding roots of a function f 

  % a - lower bound 
  % b - upper bound 
  % f - function handle under consideration 

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
