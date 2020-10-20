#!/usr/bin/python
import numpy as np
#import cupy as np
#from numpy import linalg
import math, sys, os, time
'''
This is how you comment out a section
'''

# this is how you comment out all on the right side

print("SIM. started at ", time.ctime())
start_time = time.time()


nxyz = np.array([32, 32,  32])  
L = 70. * np.array([1., 1., 1.]) # This means between -L/2 and L/2


t_equil = 50. 		# equilibration time
deltat_equil = 0.2 	# imaginary timestep



printToStdoutEvery = 2 # print to stdout every this steps
printDenEvery      = 2 # print to file every this steps

def potential_external(x, y, z):
	return 0. # npr. x**2 + y**2 + z**2	
	
def energy_density_interaction(den):
	return  0. # energy density
	
def dEps_dPsi(den):
	return 0. #  mean-field-like interaction potential 

def init_psi(x, y, z):
	# početni oblik valne funkcje
	sigma = 1.
	psi_in = np.exp(-0.5* (x**2 + y**2 + z**2) / sigma**2)  + 0.j #Gaussian (!) bitno dodati imaginarni dio
	return psi_in
		
	
#below this point is a code. Above are all the neccesary parameters you need to set
################################################################################
################################################################################
################################################################################
################################################################################



imProp = False
#setting up auxiliary variables
LHalf = L/2.
dx = L / nxyz
dk = (2.*np.pi) / L
d3r = np.prod(dx)
x = np.linspace(-LHalf[0], LHalf[0], nxyz[0], endpoint=False)
kx = np.fft.fftfreq(nxyz[0], dx[0] / ( 2. * np.pi))
d3k = kx[1] - kx[0]
y = np.linspace(-LHalf[1], LHalf[1], nxyz[1], endpoint=False)
ky = np.fft.fftfreq(nxyz[1], dx[1] / ( 2. * np.pi))
d3k *= ky[1] - ky[0]
z = np.linspace(-LHalf[2], LHalf[2], nxyz[2], endpoint=False)
kz = np.fft.fftfreq(nxyz[2], dx[2] / ( 2. * np.pi))
d3k *= kz[1] - kz[0]
#make a mesh
x, y, z = np.meshgrid(x, y, z, indexing='ij')
kx, ky, kz = np.meshgrid(kx, ky, kz, indexing='ij')


sumk2 = (kx**2 + ky**2 + kz**2) / 2
dt_equil = -1.j*deltat_equil
kinprop = np.exp(-1.j * sumk2 * dt_equil) # kinetic propagator
pot_ext = potential_external(x, y, z) 


print(" *** Parameters of the simulation ***")
print("nxyz ", nxyz)
print("L ", L)

print("dx " , dx)
print("dk ", dk)
print("d3r " , d3r)
print("d3k", d3k )
print("\n\n")



#Trotter operator, O(dt^2) global error in energy


def T2_operator(psi, den):
	return psi #
	


	
def energy(psi, den):
	return 0. # energy

output_dir = 'snapshots_time_evolution_0'; num=0
while os.path.exists(output_dir): num+=1; output_dir="snapshots_time_evolution_"+str(num)
if not os.path.exists(output_dir): os.makedirs(output_dir)
if not os.path.exists(output_dir+'/npy_files'): os.makedirs(output_dir+'/npy_files')
file_en_equil = open(output_dir + '/en_equil.dat', 'w', buffering=1)
file_en_prod = open(output_dir + '/en_prod.dat', 'w', buffering=1)


def dft_simulation(t_max,delta_t):	
	global psi
	timestep = 0
	time = 0.
	while time <= t_max:
		den = np.absolute(psi) ** 2
		if(timestep % printToStdoutEvery == 0 or timestep==0):
			en = energy(psi, den)
			out_text = "%.5e %.5e\n" % (time, en)
			file_en_equil.write(out_text)
			print(out_text)
		#log
		if(timestep % printDenEvery == 0 or timestep==0):
			timestep /= printDenEvery
			file=output_dir + '/npy_files/psi_equil_%i' % timestep
			
			#some examples on different estimators
			np.save(file + "_xy", dx[2] * np.sum(den, axis=2) )  			# xy density
			np.save(file + "_x", np.swapaxes(den,0,2)[int(nxyz[2] / 2)][(int)(nxyz[1] / 2)])   	# slice through x
			timestep *= printDenEvery
		psi = T2_operator(psi, den)
		
		
		time += delta_t
		timestep += 1	
		
	return 0






print(" *** ...  EQUILIBRATION ... **** ")


psi  = init_psi(x, y, z)
dft_simulation(t_equil, deltat_equil)


np.save(output_dir + "/psi_final", psi)

print("SIM. ENDED at ", time.ctime())
print( time.time() - start_time)

	
