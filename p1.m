clear

%% WARNING! Lots of magic numbers tuned not to get overflow

E = 2;
V = 1;
rmax = 10;

% b vector
bi = linspace(0.0001, rmax, 20);

% Anon integrands. V and E kept constant
int1 = @(p,b,rmin) (4.*p.*b)./((p.^2 + rmin).^2 .* sqrt(1 - b.^2/((p.^2 + rmin).^2))); 
int2 = @(p,b,rmin) (4.*p.*b)./((p.^2 + rmin).^2 .* sqrt(1 - b.^2./((p.^2 + rmin).^2) - V/E)); 

% E>V
analytical = @(b) 2*asin(b/(rmax*sqrt(1-V/E))) - 2*asin(b/rmax) 

% Bisect rmin roots depending on b
for i = 1:length(bi)
  rmin(i) = bisec(0, rmax, bi(i), V, E, @r_min, 1e-5);
end

% Filter out roots using a tolerance
rmin = rmin(rmin<rmax-1e-4)

for i = 1:length(rmin)
  I1(i) = adaptive_simpsons(int1, sqrt(bi(i)-rmin(i)), sqrt(rmax-rmin(i)), bi(i)+1e-3, rmin(i), 1e-10,0,32);
  I2(i) = adaptive_simpsons(int2, 1e-8, sqrt(rmax-rmin(i)), bi(i)+1e-3, rmin(i), 1e-10,0,32);
end
subplot(2,1,1)
plot(bi(1:length(rmin)), analytical(bi(1:length(rmin))), '-')
grid on
hold on
plot(bi(1:length(rmin)), (I1-I2), 'rs')
subplot(2,1,2)
semilogy(bi(1:length(rmin)), abs(analytical(bi(1:length(rmin)))-(I1-I2)), '.-')
grid on
hold off

% Parameters par1, par2 are b and rmin
function I = adaptive_simpsons(f,a,b,par1,par2,tol,depth,maxdepth)

  m = (a + b) / 2;
  fa = f(a,par1,par2);
  fb = f(b,par1,par2);
  fm = f(m,par1,par2);

  I1 = (b-a) / 6 * (fa + 4*fm + fb);

  I2a = (m-a) / 6 * (fa + 4*f((a+m)/2,par1,par2) + fm);
  I2b = (b-m) / 6 * (fm + 4*f((m+b)/2,par1,par2) + fb);
  I2 = I2a + I2b;

  err = abs(I2 - I1);

  if err <= tol || depth >= maxdepth
    I = I2;
    return;
  else
    Ia = adaptive_simpsons(f,a,m,par1,par2,tol/2,depth+1,maxdepth);
    Ib = adaptive_simpsons(f,m,b,par1,par2,tol/2,depth+1,maxdepth);
    I = Ia + Ib;
  end
end

% implementing bisec. method for finding the minimum of r_min
function root = bisec(l, u, b, V, E, f, tol)

% l - lower bound 
% u - upper bound 
% b - impact parameter 
% V - potential 
% E - total energy, constant 
% f - function handle under consideration 
% tol - tolerance of calculation

while abs(l - u) > tol 
    c = l + ((u-l)/2); % finding the midpoint 
    
    if sign(f(l, b, V, E)) == sign(f(c, b, V, E))
        l = c; 
    else
        u = c; 
    end 
end 

root = l + ((u-l)/2);
    
end 

% definition of the r_min-function
function rad = r_min(r, b, V, E)

rad = 1 - ((b^2)/(r^2)) - (V/E);

end 
