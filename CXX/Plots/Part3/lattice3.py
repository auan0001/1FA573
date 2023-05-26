import matplotlib.pyplot as plt
import numpy as np
plt.rc('text',usetex = True)
plt.rc('font',size = 10)
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
fig, ax = plt.subplots(1,4)

plt.axis("off")
lat1 = np.loadtxt(path + '/S_32_0_0' + '.dat', delimiter=' ').astype(int).reshape(N,N)
lat2 = np.loadtxt(path + '/S_32_37_0' + '.dat', delimiter=' ').astype(int).reshape(N,N)
lat3 = np.loadtxt(path + '/S_32_50_0' + '.dat', delimiter=' ').astype(int).reshape(N,N)
lat4 = np.loadtxt(path + '/S_32_64_0' + '.dat', delimiter=' ').astype(int).reshape(N,N)

ax[0].imshow(lat1, cmap='tab20', interpolation='nearest')
ax[0].set_title("$k_BT=5$")

ax[1].imshow(lat2, cmap='tab20', interpolation='nearest')
ax[1].set_title("$k_BT=k_BT_c$")

ax[2].imshow(lat3, cmap='tab20', interpolation='nearest')
ax[2].set_title("$k_BT=2$")

ax[3].imshow(lat4, cmap='tab20', interpolation='nearest')
ax[3].set_title("$k_BT=1$")

fig.suptitle("$L=32$ lattice with $J=1$, $B=0$", fontsize = '16')
[axi.set_axis_off() for axi in ax.ravel()]

plt.show()
