import numpy as np
import matplotlib.pyplot as plt

def P(x,l):
    if (l == 0):
        return 1
    if (l == 1):
        return x
    else:
        return ((2*l-1)*x*P(x,l-1)-(l-1)*P(x,l-2))/l

def dP(x,l):
    if (l == 0):
        return 1
    if (l == 1):
        return x
    else:
        return (-l*x*P(x,l)+l*P(x,l-1))/(1-x*x)

def bisection(f,a,b,tol):
    if f(a)*f(b) > 0:
        return None
    c = (a+b)/2
    if np.abs(f(c)) < tol:
        return c
    elif f(c)*f(a) > 0:
        return bisection(f,c,b,tol)
    elif f(c)*f(b) > 0:
        return bisection(f,a,c,tol)

a = -1
b = 1
N = 1000

h = (b-a)/N

tol = 10e-8
order = 3

x = np.linspace(a,b,N)
f = lambda x: P(x,order)
'''
Pseudocode MATLAB:
    f - legendre function

    x_sign = sign(f(x))
    x_alt_sign = x_sign(1:N-1).*x_sign(2:N)
    x_alt_sign_idx = find(x_alt_sign < 0)
    for x in x_alt_sign_idx
       bisection(f,x,x+h,tol) 
'''


x_sign = np.sign(f(x))
x_alt_sign = x_sign[0:N-1]*x_sign[1:N]
x_alt_sign_idx = x[np.where(x_alt_sign < 0)]
x_sz = x_alt_sign_idx.size
roots = np.zeros(x_sz)
for i in range (0,x_sz):
    roots[i] = bisection(f,x_alt_sign_idx[i],x_alt_sign_idx[i]+h,tol)
print(roots)

plt.plot(x,f(x))
plt.plot(roots,f(roots),'o')
plt.xlabel('x')
plt.ylabel('P(x,l=' + str(order) + ')')
plt.grid()
plt.show()
