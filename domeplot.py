# encoding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

def loadMatrix(src, nr = 600, nc = 600):
	D = np.zeros((nr, nc))
	i = 0
	f = open(src, 'r')
	for line in f.readlines()[4:]:
		D[i] = line.split()
		i += 1
		if(i >= nr):
			break
		
	return D

def plotSR(n, srcFolder = '', nr = 600, nc = 600):
    h_src = srcFolder + '/hdata_' + str(n)
    z_src = srcFolder + '/zdata_' + str(0)
    H = loadMatrix(h_src, nr, nc) # lava thickness
    Z = loadMatrix(z_src, nr, nc) # underground plane
    
    A = np.zeros((nr, nc))
    for i in range(1, 600):
	    A[i, :] = A[i - 1, :] + 1
    Z2 = Z - A # preexisted gaps
    D = H + Z2 # synthesized topography
    print(np.max(H))
    
    row = np.linspace(0, 599, 600)
    col = np.linspace(599, 0, 600)
    xx, yy = np.meshgrid(row, col) # create a grid
    
    # make the norm:  Note the center is offset so that the land has
    # more dynamic range:
    divnorm = colors.TwoSlopeNorm(vmin=-60, vcenter=0, vmax=60)

    fig = plt.figure(figsize = (8,6), dpi=500)
    ax = fig.subplots()
    ax.set_aspect('equal') # equalize the two lengths of figure
    pcm = ax.pcolormesh(xx, yy, D, cmap = plt.cm.RdBu_r, 
                         norm = divnorm, rasterized = True) # pcolormesh the topography
    
    # figure settings
    cb = plt.colorbar(pcm, shrink=0.8, extend='both')
    cb.set_label('Height/m')
    plt.yticks(np.arange(0, D.shape[0], 100), np.arange(0, 300, 50))
    plt.ylabel('Length/km')
    plt.xticks(np.arange(0, D.shape[1], 100), np.arange(0, 300, 50))
    plt.xlabel('Length/km')
    plt.title('t=' + str(n) + 's')
    
def plotSR_T(n, srcFolder = '', nr = 600, nc = 600):
    T_src = srcFolder + '/tdata_' + str(n)
    
    T = loadMatrix(T_src, nr, nc) 
    # T_array=list(list(T))
    # min2=np.sort(T_array)[1]
    print(np.max(T))
    
    row = np.linspace(0, 599, 600)
    col = np.linspace(599, 0, 600)
    xx, yy = np.meshgrid(row, col) # create a grid
    
    # make the norm:  Note the center is offset so that the land has
    # more dynamic range:
    divnorm = colors.TwoSlopeNorm(vmin=1200, vcenter=1300, vmax=1400)

    fig = plt.figure(figsize = (8,6), dpi=300)
    ax = fig.subplots()
    ax.set_aspect('equal') # equalize the two lengths of figure
    pcm = ax.pcolormesh(xx, yy, T, cmap = plt.cm.RdBu_r, 
                         norm = divnorm, rasterized = True) # pcolormesh the topography
    
    # figure settings
    cb = plt.colorbar(pcm, shrink=0.8, extend='both')
    cb.set_label('Temperature/K')
    plt.yticks(np.arange(0, T.shape[1], 100), np.arange(0, 300, 50))
    plt.ylabel('Length/km')
    plt.xticks(np.arange(0, T.shape[0], 100), np.arange(0, 300, 50))
    plt.xlabel('Length/km')
    plt.title('t=' + str(n) + 's')

fd = r'F:/LavaSimulation/Output/SR10'

for i in range(1, 11):
    time = i * 20000
    #later = time + 100000
    plotSR(time, fd)
    # plotSR_T(time, fd)
    #plotFlowArea(time, later, fd)
    # plotSection(time,fd)

    plt.show()
# plot3D(400000,fd)
# plotCross(80000,fd)
# h_src = fd + '/hdata_' + str(100000)
# z_src = fd + '/zdata_' + str(0)

