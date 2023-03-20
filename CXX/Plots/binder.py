import matplotlib.pyplot as plt
import numpy as np
# plt.gca().set_aspect('equal')
# plt.style.use('seaborn-colorblind')
params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{amsmath}', r'\usepackage{amsfonts}', r'\usepackage{amssymb}']}
plt.rcParams.update(params)
plt.rc('text',usetex = True)
plt.rc('font',size = 13)
plt.rc('pgf', texsystem='pdflatex')  # or luatex, xelatex...

xx = np.genfromtxt('/home/auan/1FA573/CXX/Data/binder.txt',
                     skip_header=0,
                     skip_footer=0,
                     names=True,
                     dtype=float,
                     delimiter=' ',
                     usecols=(3,)).astype(float)

y1 = np.genfromtxt('/home/auan/1FA573/CXX/Data/binder.txt',
                     skip_header=0,
                     skip_footer=0,
                     names=True,
                     dtype=float,
                     delimiter=' ',
                     usecols=(0,)).astype(float)

y2 = np.genfromtxt('/home/auan/1FA573/CXX/Data/binder.txt',
                     skip_header=0,
                     skip_footer=0,
                     names=True,
                     dtype=float,
                     delimiter=' ',
                     usecols=(1,)).astype(float)

y3 = np.genfromtxt('/home/auan/1FA573/CXX/Data/binder.txt',
                     skip_header=0,
                     skip_footer=0,
                     names=True,
                     dtype=float,
                     delimiter=' ',
                     usecols=(2,)).astype(float)

fig, ax = plt.subplots()
plt.grid(linestyle='dashed')
plt.title(r'4th Order Cumulant (Binder parameter)', fontsize = '22')
plt.plot(xx, y1, label = r'L=8', marker='x', linestyle='--', linewidth='1')
plt.plot(xx, y2, label = r'L=16', marker='x', linestyle='--', linewidth='1')
plt.plot(xx, y3, label = r'L=32', marker='x', linestyle='--', linewidth='1')
plt.vlines(x=2.2691853,ymin=min(y3),ymax=max(y3), label = r'$T_c$', color='grey', linestyle='--')
plt.xlabel(r'$\frac{k_BT}{J}$', fontsize = '24')
plt.ylabel(r'$U_L$', fontsize = '24')
ratio = 1.0
xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)
plt.grid(linestyle='dashed')
plt.legend(loc='best', fontsize=15) 
plt.show()
