import matplotlib.pyplot as plt
import numpy as np
import pylab 
from scipy import interpolate

'''
import matplotlib.pylab as pylab
params = {'legend.fontsize': 12.,
          'figure.figsize': (     7.6*0.8     , 4.8),
         'axes.labelsize': 25,
         'axes.titlesize': 25,
         'xtick.labelsize': 20,
         'ytick.labelsize': 20}
pylab.rcParams.update(params)
'''

import os



main_dir = "snapshots_time_evolution_0"
nxyz = np.array([128, 128,  128])   #if you want 1D or 2D, just change this array
#L_ok = np.array([9000., 4500., 4500.]) #you add on top of L an interface which is an absorbing area. L_ok is area without the interface
L = 30000. * np.array([1.5, 1.5, 1.5] ) / 2 #if you wan

LHalf = L/2.

x = np.linspace(-LHalf[0], LHalf[0], nxyz[0], endpoint=False)
y = np.linspace(-LHalf[1], LHalf[1], nxyz[1], endpoint=False)
z = np.linspace(-LHalf[2], LHalf[2], nxyz[2], endpoint=False)
xa, ya = np.meshgrid(x, y, indexing='ij')

timestep = 0
zmaxc = -1.;

minl = np.min(L)
params = {'figure.figsize': (7.8*L[0]/minl, 7.8*L[1]/minl),
 'xtick.labelsize': 22,
         'ytick.labelsize': 22,
         'axes.titlesize': 22, }
pylab.rcParams.update(params)
pylab.xlim(-LHalf[0], LHalf[0])
pylab.ylim(-LHalf[1], LHalf[1])


tid = 1

try:
	os.mkdir("2D_snapshots") 
except:
	print("directory slice already exists !!")
	
	
while(True):
	#print(timestep)
	z1 = np.load(main_dir+"/npy_files/psi_prod_" + str(timestep) + "_xy_c1.npy")
	z2 = np.load(main_dir+"/npy_files/psi_prod_" + str(timestep) + "_xy_c2.npy")
	z = z1
	print(np.shape(x), np.shape(y), np.shape(z))
	
	f = interpolate.interp2d(x, y, np.swapaxes(z, 0, 1), kind='linear')
	xx = np.linspace(-LHalf[0], LHalf[0], 2* nxyz[0], endpoint=False)
	yy = np.linspace(-LHalf[1], LHalf[1], 2 *nxyz[1], endpoint=False)
	
	z = f(xx, yy)
	
	if(zmaxc < 0.):
		zmaxc = np.max(z) / 20
	#z = z * (z > (zmaxc / 5.))
	z_min, z_max = 0, zmaxc #0., 0.002
	fig, ax = plt.subplots()
	c = ax.pcolormesh(xx, yy, z, cmap='RdBu' ,vmin=z_min, vmax=z_max)
	#ax.set_title("Nt=%.0f, v=%.1f,   t = %.2f ms" % (npart_arr[nindx], 0.4 + vindx * 0.4 ,(timestep+1) * 1000 * 40. * 9.23649E-06))
	# set the limits of the plot to the limits of the data
	ax.axis([x.min(), x.max(), y.min(), y.max()])
	#fig.colorbar(c, ax=ax)
	pylab.xlabel("x", fontsize=20)
	pylab.ylabel("y", fontsize=20)
	#pylab.axis('equal')
	
	
	
	
	
	
	pylab.savefig("2D_snapshots/2D_snapshot_%i_xy" % (timestep / tid),  bbox_inches='tight')
	#plt.show()
	pylab.clf()		
	timestep += tid
	plt.close()
	'''
	if(timestep >= 50):
		break
	'''




