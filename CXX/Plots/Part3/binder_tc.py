import matplotlib.pyplot as plt
import numpy as np
params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{amsmath}', r'\usepackage{amsfonts}', r'\usepackage{amssymb}']}
plt.rcParams.update(params)
plt.rc('text',usetex = True)
plt.rc('pgf', texsystem='pdflatex')  # or luatex, xelatex...
plt.rc('xtick', labelsize=18) 
plt.rc('ytick', labelsize=18)

xx = np.genfromtxt('/home/auan/1FA573/CXX/Data/Part3/binder_tc.txt',
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

plt.title(r'Phase transition at $T_c$', fontsize = '22')
plt.plot(xx, y1, label = r'L=8', marker='.', linestyle='--', linewidth='1')
plt.plot(xx, y2, label = r'L=16', marker='.', linestyle='--', linewidth='1')
plt.plot(xx, y3, label = r'L=32', marker='.', linestyle='--', linewidth='1')
plt.vlines(x=2.2691853,ymin=min(y3),ymax=max(y3), label = r'$T_c$', color='grey', linestyle='--')
plt.xlabel(r'$\frac{k_BT}{J}$', fontsize = '22')
plt.ylabel(r'$U_L$', fontsize = '22')
plt.tight_layout()
plt.grid(linestyle='dashed')
plt.legend(loc='best', fontsize=18) 
plt.show()
