import matplotlib.pyplot as plt
import numpy as np
# plt.gca().set_aspect('equal')
# plt.style.use('seaborn-colorblind')
params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{amsmath}', r'\usepackage{amsfonts}', r'\usepackage{amssymb}']}
plt.rcParams.update(params)
plt.rc('text',usetex = True)
plt.rc('font',size = 13)
plt.rc('pgf', texsystem='pdflatex')  # or luatex, xelatex...

data = np.loadtxt('/home/auan/1FA573/CXX/Data/susc3.txt').astype(float)

fig, ax = plt.subplots()
plt.grid(linestyle='dashed')
plt.title(r'Susceptibility', fontsize = '22')
plt.plot(data[3,:], data[0,:], label = r'L=8', marker='x', linestyle='--', linewidth='1')
plt.plot(data[3,:], data[1,:], label = r'L=16', marker='x', linestyle='--', linewidth='1')
plt.plot(data[3,:], data[2,:], label = r'L=32', marker='x', linestyle='--', linewidth='1')
plt.xlabel(r'$\frac{k_BT}{J}$', fontsize = '24')
plt.ylabel(r'$\chi$', fontsize = '24')
ratio = 1.0
xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)
plt.grid(linestyle='dashed')
plt.legend(loc='best', fontsize=15) 
plt.show()
