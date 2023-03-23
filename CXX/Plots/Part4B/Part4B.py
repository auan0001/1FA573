import matplotlib.pyplot as plt
import numpy as np
# plt.style.use('seaborn-colorblind')
params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{amsmath}', r'\usepackage{amsfonts}', r'\usepackage{amssymb}']}
plt.rcParams.update(params)
plt.rc('text',usetex = True)
plt.rc('pgf', texsystem='pdflatex')  # or luatex, xelatex...
plt.rc('xtick', labelsize=18) 
plt.rc('ytick', labelsize=18)

# Lattice dims with base 2
start = 1
stop = 10
N = 32
part1 = '3'
part2 = '4B'
col = {'temp' : 0,
       'order': 1,
       'chi': 2,
       'cb': 3,
       'u': 4}

# Path to stored data
path1 = '/home/auan/1FA573/CXX/Data/Part' + str(part1)
path2 = '/home/auan/1FA573/CXX/Data/Part' + str(part2)

fig, ax = plt.subplots(1,3)
sym = col['order']
x = np.genfromtxt(path1 + '/ising_' + str(N) + '_' + str(part1) + '.dat',
                  skip_header=1,
                  skip_footer=0,
                  names=True,
                  dtype=float,
                  delimiter=' ',
                  usecols=(col['temp'],)).astype(float)

ym0 = np.genfromtxt(path1 + '/ising_' + str(N) + '_' + str(part1) + '.dat',
                  skip_header=1,
                  skip_footer=0,
                  names=True,
                  dtype=float,
                  delimiter=' ',
                  usecols=(sym,)).astype(float)
ax[0].plot(x, ym0, label = r'$L_H='+str(N)+'$', marker='.', linestyle='--', linewidth='1', markersize='3')

yh0 = np.genfromtxt(path2 + '/ising_' + str(N) + '_' + str(part2) + '.dat',
                  skip_header=1,
                  skip_footer=0,
                  names=True,
                  dtype=float,
                  delimiter=' ',
                  usecols=(sym,)).astype(float)
ax[0].plot(x, yh0, label = r'$L_H='+str(N)+'$', marker='.', linestyle='--', linewidth='1', markersize='3')

ax[0].grid()
ax[0].set_title(r'Order parameter', fontsize = '18')
ax[0].set_xlabel(r'$k_BT$', fontsize = '18')
ax[0].set_ylabel(r'$|M|$', fontsize = '18')


sym = col['chi']
x = np.genfromtxt(path1 + '/ising_' + str(N) + '_' + str(part1) + '.dat',
                  skip_header=1,
                  skip_footer=0,
                  names=True,
                  dtype=float,
                  delimiter=' ',
                  usecols=(col['temp'],)).astype(float)

ym1 = np.genfromtxt(path1 + '/ising_' + str(N) + '_' + str(part1) + '.dat',
                  skip_header=1,
                  skip_footer=0,
                  names=True,
                  dtype=float,
                  delimiter=' ',
                  usecols=(sym,)).astype(float)
ax[1].plot(x, ym1, label = r'$L_H='+str(N)+'$', marker='.', linestyle='--', linewidth='1', markersize='3')

yh1 = np.genfromtxt(path2 + '/ising_' + str(N) + '_' + str(part2) + '.dat',
                  skip_header=1,
                  skip_footer=0,
                  names=True,
                  dtype=float,
                  delimiter=' ',
                  usecols=(sym,)).astype(float)
ax[1].plot(x, ym1, label = r'$L_H='+str(N)+'$', marker='.', linestyle='--', linewidth='1', markersize='3')

ax[1].grid()
ax[1].set_title(r'Susceptibility', fontsize = '18')
ax[1].set_xlabel(r'$k_BT$', fontsize = '18')
ax[1].set_ylabel(r'$\chi$', fontsize = '18')

sym = col['cb']
x = np.genfromtxt(path1 + '/ising_' + str(N) + '_' + str(part1) + '.dat',
                  skip_header=1,
                  skip_footer=0,
                  names=True,
                  dtype=float,
                  delimiter=' ',
                  usecols=(col['temp'],)).astype(float)

ym2 = np.genfromtxt(path1 + '/ising_' + str(N) + '_' + str(part1) + '.dat',
                  skip_header=1,
                  skip_footer=0,
                  names=True,
                  dtype=float,
                  delimiter=' ',
                  usecols=(sym,)).astype(float)
ax[2].plot(x, ym2, label = r'$L_M='+str(N)+'$', marker='.', linestyle='--', linewidth='1', markersize='3')

yh2 = np.genfromtxt(path2 + '/ising_' + str(N) + '_' + str(part2) + '.dat',
                  skip_header=1,
                  skip_footer=0,
                  names=True,
                  dtype=float,
                  delimiter=' ',
                  usecols=(sym,)).astype(float)
ax[2].plot(x, yh2, label = r'$L_H='+str(N)+'$', marker='.', linestyle='--', linewidth='1', markersize='3')

ax[2].grid()
ax[2].set_title(r'Specific heat', fontsize = '18')
ax[2].set_xlabel(r'$k_BT$', fontsize = '18')
ax[2].set_ylabel(r'$C_B$', fontsize = '18')

fig.suptitle(r'Comparison of methods for $L=' + str(N) + '$, $J=0$ over $200000$ sweeps', fontsize='22')
lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox = True, shadow = True)
fig.tight_layout()
plt.show()
