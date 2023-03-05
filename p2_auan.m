% V : r dependent
% a : 1 (should not matter? for the LJ-pot)
% b : impact, varies as in Part 1 on [bmin, rmax]
% rmax : 3*a (rmax = 3 if a = 1)
% E = 0.1V to 100V (use same spacing as for b)
% Remark: No analytical solutions exists for the deflection
% angle using a radial L-J potential.

% 1. L-J POTENTIAL
% 1.1 Find rmin from radial dependent V-term in denom square root term
% 1.2 Redefine integrand 2 with radial dep V-term
% 1.3 Integrate all triples {E_i, b_i, rmin_i}

% 2 CROSS-SECTION FROM L-J POTENTIAL
% 2.1 Use diff(theta_n) to get a vector of size N-1.
% 2.2 
clear

% TeX for plots
set(0, 'defaulttextinterpreter', 'latex')
hAxes.TickLabelInterpreter = 'latex';

% Parameters
N = 100;
Vmin = 0.1;
Vmax = 100;
rmax = 3;

bmin = 1e-2;
rmax_numerical = rmax-1e-4;

% Array containing N impact b-values
impact = linspace(bmin,rmax_numerical,N);
E = linspace(Vmin,Vmax,N);

% Integrate for each {E_i, b_i (impact), rmin_i}
for i = 1:N
  rmin(i) = search(rmax, 1e-12, @denominator, @FNV, E(i), 0.2, impact(i));
  I1(i) = rectangle_rule( @FNI1, impact(i), sqrt(rmax - impact(i)), @FNV, E(i), impact(i), N);
  I2(i) = rectangle_rule( @FNI2, rmin(i), sqrt(rmax - rmin(i)), @FNV, E(i), impact(i), N);
end
theta_n = 2*impact.*(I1-I2);

% Indexing helper since finite difference results in N->N-1 points
idx = 1:length(theta_n)-1;

% Find indexes just before theta > 0
idx0 = 1:(find(theta_n < 0, 1)-1);

% Derivatives
dtheta_db = diff(theta_n);
db_dtheta = 1./diff(theta_n);

% Cross-section
cross_sec = impact(idx0)./sin(theta_n(idx0)).*abs(db_dtheta(idx0));

% Plot skip param
sk = 1:N/50:N;

% Plots
subplot(2,2,1)
plot(impact(sk), theta_n(sk), '.-')
grid on
hold on
hold off
title('Scattering angle for Lennard-Jones potential')
xlabel('$E \in [0.1V_0, 100V_0]$')
ylabel('$\theta$')
legend('Numerical', 'Location', 'Best')

subplot(2,2,2)
plot(theta_n(idx0), cross_sec, '.-')
hold on
grid on
hold off
xlabel('$\theta > 0$')
ylabel('$\frac{d\sigma}{d\Omega}=\frac{b}{\sin(\theta)}|\frac{db}{d\theta}|$')
title('Scattering cross section')

subplot(2,2,3)
plot(theta_n(idx0), dtheta_db(idx0), '.-')
hold on
grid on
hold off
xlabel('$\theta > 0$')
ylabel('$\frac{d\theta}{db}$')
title('Computed derivative')

subplot(2,2,4)
plot(theta_n(idx0), abs(db_dtheta(idx0)), '.-')
hold on
grid on
hold off
xlabel('$\theta > 0$')
ylabel('$|\frac{db}{d\theta}|$')
title('Computed absolute inverse derivative')

% rectangle rule for function f on domain
% transformed from [r1,r2] -> [0,(r2-r1)^1/2] 
function int = rectangle_rule(f,r1,r2,V,E,b,N)
  h = r2/N;
  s = 0; % s as in sum
  for iu = 1:N
    u = h*(iu-0.5);
    r = u^2 + r1;
    s = s + u*f(r,b,V,E);
    int = 2*h*s;
  end
end

% Shooting search method used to find rmin
function x0 = search(x0, tol, f, V, E, dx, b)
  while (dx > tol)
    x0 = x0 - dx;
    if(f(x0,b,V,E) < 0)
      x0 = x0 + dx;
      dx = dx/2;
    end
  end
end

function f = FNV(r)
  f = 8*((r).^(-12)-(r).^(-6));
end

function f = FNI1(r,b,~,~)
  f = 1/r^2 / sqrt(1-(b/r)^2);
end

function f = FNI2(r,b,V,E)
  f = 1/r^2 / sqrt(1-(b/r)^2-(V(r)/E));
end

function f = denominator(r,b,V,E)
  f = 1-(b/r)^2-(V(r)/E);
end

% Analytical
function f = FN(r,b,V,E,rmax)
  f = 2*asin(b./(rmax*sqrt(1-V(r)./E))) - 2*asin(b./rmax);
end
