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
low_dim = 3
high_dim = 3
part1 = '4B'
part2 = '3'
col = {'temp' : 0,
       'order': 1,
       'chi': 2,
       'cb': 3,
       'u': 4}

# Path to stored data
path1 = '/home/auan/1FA573/CXX/Data/Part' + str(part1)
path2 = '/home/auan/1FA573/CXX/Data/Part' + str(part2)

for i in range(low_dim,high_dim+1):
    N = 2**i
    sym = col['cb']
    x = np.genfromtxt(path1 + '/ising_' + str(N) + '_' + str(part1) + '.dat',
                      skip_header=1,
                      skip_footer=0,
                      names=True,
                      dtype=float,
                      delimiter=' ',
                      usecols=(col['temp'],)).astype(float)

    yg = np.genfromtxt(path1 + '/ising_' + str(N) + '_' + str(part1) + '.dat',
                       skip_header=1,
                       skip_footer=0,
                       names=True,
                       dtype=float,
                       delimiter=' ',
                       usecols=(sym,)).astype(float)

    # ym = np.genfromtxt(path2 + '/ising_' + str(N) + '_' + str(part2) + '.dat',
    #                    skip_header=1,
    #                    skip_footer=0,
    #                    names=True,
    #                    dtype=float,
    #                    delimiter=' ',
    #                    usecols=(sym,)).astype(float)

    plt.plot(x, yg, label = r'$L_H$='+str(N), marker='.', linestyle='--', linewidth='1')
    # plt.plot(x, ym, label = r'$L_M$='+str(N), marker='.', linestyle='--', linewidth='1')
plt.title(r'Order parameter', fontsize = '22')
plt.xlabel(r'$\frac{k_BT}{J}$', fontsize = '22')
plt.ylabel(r'$|M|$', fontsize = '22')
plt.grid(linestyle='dashed')
plt.legend(loc='best', fontsize=18) 
plt.tight_layout()
plt.show()
