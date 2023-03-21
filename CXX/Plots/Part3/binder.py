import matplotlib.pyplot as plt
import numpy as np
# plt.style.use('seaborn-colorblind')
params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{amsmath}', r'\usepackage{amsfonts}', r'\usepackage{amssymb}']}
plt.rcParams.update(params)
plt.rc('text',usetex = True)
# plt.rc('font',size = 15)
plt.rc('pgf', texsystem='pdflatex')  # or luatex, xelatex...
plt.rc('xtick', labelsize=18) 
plt.rc('ytick', labelsize=18)

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

    plt.plot(x, y, label = r'L='+str(N), marker='.', linestyle='--', linewidth='1')

plt.title(r'4th Order Cumulant (Binder parameter)', fontsize = '22')
plt.xlabel(r'$\frac{k_BT}{J}$', fontsize = '22')
plt.ylabel(r'$U_L$', fontsize = '22')
plt.grid(linestyle='dashed')
plt.legend(loc='best', fontsize=18) 
plt.tight_layout()
plt.show()
