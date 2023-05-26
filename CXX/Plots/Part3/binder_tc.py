import matplotlib.pyplot as plt
import numpy as np
plt.rc('text',usetex = True)
plt.rc('pgf', texsystem='pdflatex')  # or luatex, xelatex...
plt.rc('xtick', labelsize=18) 
plt.rc('ytick', labelsize=18)
fig, ax = plt.subplots(1,2)


# Lattice dims with base 2
low_dim = 3
high_dim = 5
part = 3
col = {'temp' : 0,
       'order': 1,
       'chi': 2,
       'cb': 3,
       'u': 4}

# Path to stored data
path = '/home/auan/1FA573/CXX/Data/Part' + str(part)

for i in range(low_dim,high_dim+1):
    N = 2**i
    x = np.genfromtxt(path + '/ising_' + str(N) + '_' + str(part) + '.dat',
                      skip_header=1,
                      skip_footer=0,
                      names=True,
                      dtype=float,
                      delimiter=' ',
                      usecols=(col['temp'],)).astype(float)

    y = np.genfromtxt('/home/auan/1FA573/CXX/Data/Part3/ising_' + str(N) + '_' + str(part) + '.dat',
                       skip_header=1,
                       skip_footer=0,
                       names=True,
                       dtype=float,
                       delimiter=' ',
                       usecols=(col['u'],)).astype(float)

    ax[0].plot(x, y, label = r'L='+str(N), marker='.', linestyle='--', linewidth='1', markersize='3')

ax[0].vlines(x=2.2691853,ymin=min(y),ymax=max(y), label = r'$T_c$', color='grey', linestyle='--')
ax[0].set_title(r'4th Order Cumulant (Binder parameter)', fontsize = '18')
ax[0].set_xlabel(r'$\frac{k_BT}{J}$', fontsize = '18')
ax[0].set_ylabel(r'$U_L$', fontsize = '18')
ax[0].grid(linestyle='dashed')
plt.legend(loc='best', fontsize=18) 

xz = np.genfromtxt('/home/auan/1FA573/CXX/Data/Part3/binder_tc.txt',
                     skip_header=0,
                     skip_footer=0,
                     names=True,
                     dtype=float,
                     delimiter=' ',
                     usecols=(3,)).astype(float)

y1 = np.genfromtxt('/home/auan/1FA573/CXX/Data/Part3/binder_tc.txt',
                     skip_header=0,
                     skip_footer=0,
                     names=True,
                     dtype=float,
                     delimiter=' ',
                     usecols=(0,)).astype(float)

y2 = np.genfromtxt('/home/auan/1FA573/CXX/Data/Part3/binder_tc.txt',
                     skip_header=0,
                     skip_footer=0,
                     names=True,
                     dtype=float,
                     delimiter=' ',
                     usecols=(1,)).astype(float)

y3 = np.genfromtxt('/home/auan/1FA573/CXX/Data/Part3/binder_tc.txt',
                     skip_header=0,
                     skip_footer=0,
                     names=True,
                     dtype=float,
                     delimiter=' ',
                     usecols=(2,)).astype(float)

ax[1].set_title(r'Phase transition around $T_c$', fontsize = '18')
ax[1].plot(xz, y1, label = r'L=8', marker='.', linestyle='--', linewidth='1', markersize='3')
ax[1].plot(xz, y2, label = r'L=16', marker='.', linestyle='--', linewidth='1', markersize='3')
ax[1].plot(xz, y3, label = r'L=32', marker='.', linestyle='--', linewidth='1', markersize='3')
ax[1].vlines(x=2.2691853,ymin=min(y3),ymax=max(y3), label = r'$T_c$', color='grey', linestyle='--')
ax[1].set_xlabel(r'$\frac{k_BT}{J}$', fontsize = '18')
ax[1].set_ylabel(r'$U_L$', fontsize = '18')
plt.tight_layout()
plt.grid(linestyle='dashed')
plt.legend(loc='best', fontsize=18) 
fig.suptitle(r'Approximating phase transition for $J=1$, $B=0$ over $200000$ sweeps', fontsize='22')
lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox = True, shadow = True)
plt.show()
