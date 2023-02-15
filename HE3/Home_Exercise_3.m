a = -1;
b = 1;
N = 100;

h = (b-a)/N;

tol = 10e-8;
order = 3;

x = linspace(a,b,N);

Pn = P(x,3);
plot(x,Pn);

function P_n = P(x,l)

  if (l == 0)
    P_n = 1;
  end

  if (l == 1)
    P_n = x;
  else
    P_n = ((2*l)*x.*P(x,l)-(l).*P(x,l-1))./l
  end

end
