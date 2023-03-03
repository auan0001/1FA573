clear

% TeX for plots
set(0, 'defaulttextinterpreter', 'latex')
hAxes.TickLabelInterpreter = 'latex';

% Spatial domain
N = 100;
rmax = 10;
bmin = 1e-1;
rmax_numerical = rmax-1e-6;

% Array containing N impact b-values
impact = linspace(bmin,rmax_numerical,N);

% Constant potential and energy
V = 1;
E = 3*V;

% Integrands
FNI1 = @(r,b) 1/r^2 / sqrt(1-(b/r)^2);
FNI2 = @(r,b) 1/r^2 / sqrt(1-(b/r)^2-(V/E));

% Square root term in the denominator of FNI2
denominator = @(r,b) 1-(b/r).^2-(V/E);

% Analytical solution
FN = @(b) 2*asin(b/(rmax*sqrt(1-V/E))) - 2*asin(b/rmax);

% Integrate for each rmin
for i = 1:N
  rmin = search(rmax, 1e-12, denominator, 0.2, impact(i));
  I1(i) = rectangle_rule(FNI1,impact(i),sqrt(rmax - impact(i)),impact(i),N);
  I2(i) = rectangle_rule(FNI2,rmin, sqrt(rmax - rmin),impact(i),N);
end

% Numerical solution
theta_n = 2*impact.*(I1-I2);
theta = real(FN(impact));

% Plots
subplot(1,2,1)
plot(impact, theta)
grid on
hold on
plot(impact,theta_n,'.')
hold off
title('Scattering angle')
xlabel('$b$')
ylabel('$\theta$')
legend('Analytical', 'Numerical', 'Location', 'Best')
subplot(1,2,2)
semilogy(impact,abs(theta_n-theta),'*-')
title('Absolute error')
xlabel('$b$')
ylabel('$L^1(\theta)$')
hold off
grid on

% rectangle rule for function f on domain
% transformed from [r1,r2] -> [0,(r2-r1)^1/2] 
function int = rectangle_rule(f,r1,r2,b,N)
  h = r2/N;
  s = 0; % s as in sum
  for iu = 1:N
    u = h*(iu-0.5);
    r = u^2 + r1;
    s = s + u*f(r,b);
    int = 2*h*s;
  end
end

% Shooting search method used to find rmin
function x0 = search(x0, tol, f, dx, b)
  while (dx > tol)
    x0 = x0 - dx;
    if(f(x0,b) < 0)
      x0 = x0 + dx;
      dx = dx/2;
    end
  end
end
