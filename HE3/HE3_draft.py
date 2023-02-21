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

# Integration domain and grid size
a = -1
b = 1
N = 300
x = np.linspace(a,b,N)

# Step size
h = (b-a)/N

# Bisection tolerance
tol = 10e-5

# Legendre polynomial order
order = 14

# Legendre polynomials, wrapper functions
P_n = lambda x: P(x,order)
P_n_prime = lambda x: dP(x,order)
f = lambda x: np.sqrt(1-x*x)

# Analytical answer to integral
I_analytical = np.pi/2

'''
Pseudocode for finding roots in MATLAB:
    f - legendre function

    x_sign = sign(f(x))
    x_alt_sign = x_sign(1:N-1).*x_sign(2:N)
    x_alt_sign_idx = find(x_alt_sign < 0)
    for x in x_alt_sign_idx
       bisection(f,x,x+h,tol) 
'''


# Finding roots of Legendre polynomial
x_sign = np.sign(P_n(x))
x_alt_sign = x_sign[0:N-1]*x_sign[1:N]
x_alt_sign_idx = x[np.where(x_alt_sign < 0)]
x_sz = x_alt_sign_idx.size
x_n = np.zeros(x_sz)
for i in range (0,x_sz):
    x_n[i] = bisection(P_n,x_alt_sign_idx[i],x_alt_sign_idx[i]+h,tol)
print(x_n.size)

# Weights
w_n = 2/((1-x_n*x_n)*P_n_prime(x_n)*P_n_prime(x_n))

# Calculating the integral, printing absolute error
I_numerical = np.sum(w_n*f(x_n))
print(I_numerical)
print(abs(I_numerical-I_analytical))

# Plots

# plt.plot(x,f(x))
# plt.plot(x_n,f(x_n),'o')
# plt.xlabel('x')
# plt.ylabel('P(x,l=' + str(order) + ')')
# plt.grid()
# plt.show()
