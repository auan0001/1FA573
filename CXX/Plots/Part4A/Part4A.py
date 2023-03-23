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
N = 16
part1 = '4A'
col = {'temp' : 0,
       'order': 1,
       'chi': 2,
       'cb': 3,
       'u': 4}

# Path to stored data
path1 = '/home/auan/1FA573/CXX/Data/Part' + str(part1)

fig, ax = plt.subplots(1,3)
for i in range(start,stop+1,2):
    sym = col['order']
    x = np.genfromtxt(path1 + '/ising_' + str(N) + '_' + str(i) + '_' + str(part1) + '.dat',
                      skip_header=1,
                      skip_footer=0,
                      names=True,
                      dtype=float,
                      delimiter=' ',
                      usecols=(col['temp'],)).astype(float)

    y = np.genfromtxt(path1 + '/ising_' + str(N) + '_' + str(i) + '_' + str(part1) + '.dat',
                       skip_header=1,
                       skip_footer=0,
                       names=True,
                       dtype=float,
                       delimiter=' ',
                       usecols=(sym,)).astype(float)
    ax[0].plot(x, y, label = r'$B='+str(i)+'$', marker='.', linestyle='--', linewidth='1', markersize='3')
ax[0].grid()
ax[0].set_title(r'Order parameter', fontsize = '18')
ax[0].set_xlabel(r'$k_BT$', fontsize = '18')
ax[0].set_ylabel(r'$|M|$', fontsize = '18')


for i in range(start,stop+1,2):
    sym = col['chi']
    x = np.genfromtxt(path1 + '/ising_' + str(N) + '_' + str(i) + '_' + str(part1) + '.dat',
                      skip_header=1,
                      skip_footer=0,
                      names=True,
                      dtype=float,
                      delimiter=' ',
                      usecols=(col['temp'],)).astype(float)

    y = np.genfromtxt(path1 + '/ising_' + str(N) + '_' + str(i) + '_' + str(part1) + '.dat',
                       skip_header=1,
                       skip_footer=0,
                       names=True,
                       dtype=float,
                       delimiter=' ',
                       usecols=(sym,)).astype(float)
    ax[1].plot(x, y, label = r'$B='+str(i)+'$', marker='.', linestyle='--', linewidth='1', markersize='3')
ax[1].grid()
ax[1].set_title(r'Susceptibility', fontsize = '18')
ax[1].set_xlabel(r'$k_BT$', fontsize = '18')
ax[1].set_ylabel(r'$\chi$', fontsize = '18')

for i in range(start,stop+1,2):
    sym = col['cb']
    x = np.genfromtxt(path1 + '/ising_' + str(N) + '_' + str(i) + '_' + str(part1) + '.dat',
                      skip_header=1,
                      skip_footer=0,
                      names=True,
                      dtype=float,
                      delimiter=' ',
                      usecols=(col['temp'],)).astype(float)

    y = np.genfromtxt(path1 + '/ising_' + str(N) + '_' + str(i) + '_' + str(part1) + '.dat',
                       skip_header=1,
                       skip_footer=0,
                       names=True,
                       dtype=float,
                       delimiter=' ',
                       usecols=(sym,)).astype(float)
    ax[2].plot(x, y, label = r'$B='+str(i)+'$', marker='.', linestyle='--', linewidth='1', markersize='3')
ax[2].grid()
ax[2].set_title(r'Specific heat', fontsize = '18')
ax[2].set_xlabel(r'$k_BT$', fontsize = '18')
ax[2].set_ylabel(r'$C_B$', fontsize = '18')

fig.suptitle(r'Varying $B$ with $L=' + str(N) + '$, $J=-1$ over $200000$ sweeps', fontsize='22')
lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox = True, shadow = True)
fig.tight_layout()
plt.show()
