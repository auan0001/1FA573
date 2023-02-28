clear

E = 2;
V = 1;
rmax = 20;
N = 1000;

% To find rmin by bisection
% denom = @(r,b) 1-r.^2./b.^2 - V/E;

% Anon funcs
int1 = @(r,b) 2*b*(1./(r.^2)) .* sqrt(1 - ((b^2)./(r.^2))); 
int2 = @(r,b) 2*b*(1./(r.^2)) .* sqrt(1 - ((b^2)./(r.^2)) - (V/E)); 

% E>V
% analytical = @(b) 2*asin(b/(rmax*sqrt(1-V/E)))- 2*asin(b/rmax) 
% analytical = @(b) b/(rmax*sqrt(1-V/E)) 

% Equally spaced b wrt to rmax
b = linspace(0.1, rmax, 20);

% rmin found analytically (?)
rmin = b*sqrt(1-E/V);
%
for i = 1:length(rmin)
  plot(b(i), bode(int1, b(i), rmax, N, b(i)), '.')
  hold on
end

function I = bode(f,r1,r2,N,b)
  h = (r2-r1)/N;
  x = r1:4*h:r2-3*h;
  I = 2*h/45*sum((7*f(x,b) + 32*f(x+h,b) + 12*f(x+2*h,b)+ 32*f(x+3*h,b) + 7*f(x+4*h,b)));
end
