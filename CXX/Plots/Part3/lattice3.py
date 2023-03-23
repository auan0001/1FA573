import matplotlib.pyplot as plt
import numpy as np
# plt.style.use('seaborn-colorblind')
params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{amsmath}', r'\usepackage{amsfonts}', r'\usepackage{amssymb}']}
plt.rcParams.update(params)
plt.rc('text',usetex = True)
plt.rc('font',size = 15)
plt.rc('pgf', texsystem='pdflatex')  # or luatex, xelatex...
plt.rc('xtick', labelsize=18) 
plt.rc('ytick', labelsize=18)

part = 3
N = 32
col = {'temp' : 0,
       'order': 1,
       'chi': 2,
       'cb': 3,
       'u': 4}

# Path to stored data
path = '/home/auan/1FA573/CXX/Data/Part' + str(part)
fig, ax = plt.subplots(2,2)

plt.axis("off")
lat1 = np.loadtxt(path + '/S_32_0_0' + '.dat', delimiter=' ').astype(int).reshape(N,N)
lat2 = np.loadtxt(path + '/S_32_37_0' + '.dat', delimiter=' ').astype(int).reshape(N,N)
lat3 = np.loadtxt(path + '/S_32_50_0' + '.dat', delimiter=' ').astype(int).reshape(N,N)
lat4 = np.loadtxt(path + '/S_32_64_0' + '.dat', delimiter=' ').astype(int).reshape(N,N)

s1 = ax[0,0].imshow(lat1, cmap='tab20', interpolation='nearest')
ax[0,0].set_title("$k_BT=5$")

ax[0,1].imshow(lat2, cmap='tab20', interpolation='nearest')
ax[0,1].set_title("$k_BT=T_c$")

ax[1,0].imshow(lat3, cmap='tab20', interpolation='nearest')
ax[1,0].set_title("$k_BT=2$")

ax[1,1].imshow(lat4, cmap='tab20', interpolation='nearest')
ax[1,1].set_title("$k_BT=1$")

fig.suptitle("$L=32$ lattice with $J=1$, $B=0$")
[axi.set_axis_off() for axi in ax.ravel()]
plt.show()
