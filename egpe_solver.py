#!/usr/bin/python
cudaEnable = False # enable CUDA ?
if(cudaEnable):
	import cupy as np
	import numpy
else:
	import numpy
	import numpy as np
import os, time

if(cudaEnable):
	memory_pool = np.cuda.MemoryPool()
	np.cuda.set_allocator(memory_pool.malloc)
	pinned_memory_pool = np.cuda.PinnedMemoryPool()
	np.cuda.set_pinned_memory_allocator(pinned_memory_pool.malloc)
	
	
print("SIM. started at ", time.ctime())
start_time = time.time()

#physics
a12 = -53.3 / 73.3
a22 =  33.9 / 73.3
kappa_3 =  0. * 2.26497E+03 # Kappa_3 coefficient


def N_over_ntilde(a12, a22):
	return (3.*np.sqrt(6.)/(5.*np.pi**2)) * (1. + np.sqrt(a22))**5. / np.abs(a12 + np.sqrt(a22))**(5./2)
	
def healing_length(a12, a22):
		return (8.*np.sqrt(6.)/(5. * np.pi)) * np.sqrt(a22) * np.power(1. + np.sqrt(a22), 3.) / np.abs(a12 + np.sqrt(a22))**(3./2)

def v_over_vtilde(a12, a22):
	return 1. / healing_length(a12, a22)

def num_each_species(nparticles):
	np1 = nparticles * np.sqrt(a22) / (1. + np.sqrt(a22))
	np2 = nparticles - np1
	return np1, np2
	
def rho0(a12, a22):
	#return np.power( -alpha / ( beta * gamma) , 1. / (gamma - 1) )
	return (25*np.pi/1024.) * np.power(a12 + np.sqrt(a22), 2) / (np.power(a22, 3./2) * np.power(1. + np.sqrt(a22), 4))


kappa_3 = 0.53 / 2 / ( rho0(a12, a22)**2 * healing_length(a12, a22)**2)		



	
#B=56.230 QMC functional
'''
alpha = -0.81221541
beta = 5.97406085
gamma = 1.27644308
'''
#Petrov units
x = 1.
alpha = 0.5 * 4*np.pi*(1.+x*x + 2*x*a12/np.sqrt(a22))/np.power(1 + x/np.sqrt(a22), 2);
beta  = 0.5 * (512*np.sqrt(np.pi)/15)*np.power((1+x*np.sqrt(a22)) / (1+x/np.sqrt(a22)), 5./2);
gamma = 3./2;


propagation_order = 2 #choose either 2, 6, 8 or 12th order

npart_fac = 8
npart_arr = np.linspace(50., 300, 10, endpoint=True)

#simulation and physical parameters
nparticles =  N_over_ntilde(a12, a22) * npart_arr[npart_fac] /2 #this will only be used when you do imaginary time propagation
np1, np2 = num_each_species(nparticles)

nparticles = np1 + np2
#np1 = 102. * 1000.
#np2 = 148. * 1000.
print("number of type 1 particles = %d\nnumber of type 2 particles = %d" % (np1, np2))
print("Kappa_3 = %.5f " % kappa_3)


nxyz = numpy.array([128, 128,  128])   #if you want 1D or 2D, just change this array
#L_ok = np.array([9000., 4500., 4500.]) #you add on top of L an interface which is an absorbing area. L_ok is area without the interface
L = 30000. * numpy.array([1.5, 1.5, 1.5] )  #if you want 1D or 2D, just change this array. This means between -L/2 and L/2
				#if you put L_ok=L, this means that no area has absorbing-boundary conditions
L_ok = 0.9 * L
activateAbsorbingWalls = True
fac_time = 1
t_equil = 1.E+06 /fac_time
t_prod =  15.E+06/fac_time #np.inf for infinite      #total propagation time. Count from the real (imaginary) timestep if you propagate in real (imaginary) time
deltat_prod  = 1000. # real timestep
deltat_equil = 1000. # this is absolute value of imaginary part



printToStdoutEvery = 10
printDenEvery = 200
read_from_file = False; infile = 'psi_in.npy'
k_kick_fac = 5
k_kick_ratio = 0.4 + k_kick_fac * 0.2 - (npart_fac % 2) * 0.1
print("k_kick_ratio = %.3e" % k_kick_ratio)
if(k_kick_ratio > 1.11):
	print("velocity too large")
	exit()
#exit()
manipulate_psi = False; k_kick = k_kick_ratio * v_over_vtilde(a12, a22); separation_distance = 36


ylim_max = 750. #when plotting rho(x)


print("alpha = %.3e, beta = %.3e, gamma = %.3e" % (alpha, beta, gamma))

kill_fac = 0.95

init_sigma = 1.3 * (50.*1000) / 73.3


def init_kill_matrix(x, y, z):
	#return 3. - np.exp( - ((x - L_ok[0]/2.)/(L[0]/2. - L_ok[0]/2))**20 )  - np.exp( - ((y - L_ok[1]/2.)/(L[1]/2. - L_ok[1]/2))**2 ) -np.exp(-  ((z - L_ok[2]/2.)/(L[2]/2. - L_ok[2]/2))**2 )
	#return 1 #. if you don't want to kill wavefunction
	matrix = 1.
	if(N_DIM == 1):
		matrix = np.logical_or(np.abs(x) > L_ok[0]/2.)
	if(N_DIM == 2):
		matrix = np.logical_or(np.abs(x) > L_ok[0]/2., np.abs(y) > L_ok[1]/2.)
	if(N_DIM == 3):
		matrix = np.logical_or(np.abs(x) > L_ok[0]/2., np.abs(y) > L_ok[1]/2.)
		matrix = np.logical_or(matrix,   np.abs(z) > L_ok[2]/2.)
	return np.invert(matrix) + matrix*kill_fac # this matrix is 1 inside box, and kill_fac outside


	
def potential_external(x, y, z):
	return 0. #x**2 + y**2 + z**2	


def energy_density_interaction(den1, den2):
	return 2.*np.pi*den1**2 + 2.*np.pi*a22*den2**2 + 4.*np.pi*a12*den1*den2 + (256.*np.sqrt(np.pi)/15.) * np.power(den1 + den2*a22, 5./2)
		
def dEps_dPsi(den1, den2, kappa_3_):
	d12 = (128. *np.sqrt(np.pi)/3.) * np.power(den1 + den2*a22 , 3./2)
	pot_int1 = 4.*np.pi*den1     + 4.*np.pi*a12*den2 + d12 - 1.j * kappa_3_ * den1**2
	pot_int2 = 4.*np.pi*den2*a22 + 4.*np.pi*a12*den1 + d12 * a22
	return pot_int1, pot_int2

def init_psi(x, y, z):
	sigma = init_sigma
	psi1 =  np.exp(-0.5*(x**2/sigma**2 + y**2/sigma**2 + z**2/sigma**2))  + 0.j #Gaussian
	psi2 =  np.exp(-0.5*(x**2/sigma**2 + y**2/sigma**2 + z**2/sigma**2))  + 0.j #Gaussian
	return psi1, psi2
		
def manipulate(psi1, psi2, x):
	psi1_l = np.roll(psi1, -separation_distance, axis=0) 
	psi2_l = np.roll(psi2, -separation_distance, axis=0)
	
	psi1_r = np.roll(psi1, separation_distance, axis=0) 
	psi2_r = np.roll(psi2, separation_distance, axis=0)
	
	psi1 = psi1_l * np.exp(1.j * k_kick * x) + psi1_r * np.exp(-1.j * k_kick * x)
	psi2 = psi2_l * np.exp(1.j * k_kick * x) + psi2_r * np.exp(-1.j * k_kick * x)
	
	return psi1, psi2
	
	
	
	
	

#below this point is a code. Above are all the neccesary parameters you need to set
################################################################################
################################################################################
################################################################################
################################################################################

print("Healing length is %.3e" % healing_length(a12, a22))
print("N over ntilde is %.6e" % N_over_ntilde(a12, a22))
print("v over vtilde is %.6e" % v_over_vtilde(a12, a22))

#exit()

imProp = False
#setting up auxiliary variables
N_DIM = len(nxyz)
LHalf = L/2.
dx = L / nxyz
dk = (2.*np.pi) / L
d3r = np.prod(dx)
#psi = np.copy(x)

x = np.linspace(-LHalf[0], LHalf[0],nxyz[0], endpoint=False)
kx = np.fft.fftfreq(nxyz[0], dx[0] / ( 2. * np.pi))
#print( " ******** x ************ ", x, kx)
print("max(kx) / v=1 = %.5e" % (np.max(kx) / v_over_vtilde(a12, a22)))
#exit()

d3k = kx[1] - kx[0]
if(N_DIM >= 2):
	y = np.linspace(-LHalf[1], LHalf[1], nxyz[1], endpoint=False)
	ky = np.fft.fftfreq(nxyz[1], dx[1] / ( 2. * np.pi))
	d3k *= ky[1] - ky[0]
if(N_DIM >= 3):
	z = np.linspace(-LHalf[2], LHalf[2], nxyz[2], endpoint=False)
	kz = np.fft.fftfreq(nxyz[2], dx[2] / ( 2. * np.pi))
	d3k *= kz[1] - kz[0]
#print( " ******** y ************ ", y, ky)
#make a mesh
if(N_DIM == 2):
	x, y = np.meshgrid(x, y, indexing='ij')
	kx, ky = np.meshgrid(kx, ky, indexing='ij')
elif(N_DIM == 3):
	x, y, z = np.meshgrid(x, y, z, indexing='ij')
	kx, ky, kz = np.meshgrid(kx, ky, kz, indexing='ij')
dt_x = np.ones(nxyz[0]) 
if(N_DIM == 2):	
	dt_y = np.ones(nxyz[1]) 
	dt = np.meshgrid(dt_x, dt_y, indexing='ij')
elif(N_DIM == 3):
	dt_y = np.ones(nxyz[1]) 
	dt_z = np.ones(nxyz[2]) 
	dt_x, dt_y, dt_z = np.meshgrid(dt_x, dt_y, dt_z, indexing='ij')

kill_matrix = init_kill_matrix(x, y, z)



dt_equil = -1.j*dt_x*deltat_equil
dt_prod = dt_x*deltat_prod + 0.j



print(" *** Parameters of the simulation ***")
print("nxyz ", nxyz)
print("L ", L)
print("N_DIM ", N_DIM)
print("dx " , dx)
print("dk ", dk)
print("d3r " , d3r)
print("d3k", d3k )
pot_ext = 0. #default value
#calculate the external potential energy density
if(N_DIM == 1):
	pot_ext = potential_external(x)
if(N_DIM == 2):
	pot_ext = potential_external(x, y)
if(N_DIM == 3):
	pot_ext = potential_external(x, y, z)
kinprop = 1.
sumk2  =0.
if(N_DIM == 3):
	sumk2 = np.array((kx**2 + ky**2 + kz**2) / 2)
if(N_DIM == 2):
	sumk2 = np.array((kx**2 + ky**2) / 2)
if(N_DIM == 1):
	sumk2 = np.array((kx**2) / 2)
#Trotter operator, O(dt^2) global error
def T2_operator(psi1, psi2, dtfac, kappa_3_, den1, den2):
	global dt, pot_ext, imProp, kinprop
	#exp(-1/2 * i dt * V)
	pot_int1, pot_int2 = dEps_dPsi(den1, den2, kappa_3_)
        psi1 *= np.exp(pot_int1 * dt)
	psi2 *= np.exp(pot_int2 * dt)
	#exp(-i * dt * T)
	
	
	psi1 = np.fft.fftn(psi1)
	psi2 = np.fft.fftn(psi2)
	psi1 *= kinprop
	psi2 *= kinprop
	psi1 = np.fft.ifftn(psi1)
	psi2 = np.fft.ifftn(psi2)
	
	
	
	den1 = psi1.real**2 + psi1.imag**2
	den2 = psi2.real**2 + psi2.imag**2
	if(imProp): #normalize		
		n1 = np.sum(den1) * d3r
		n2 = np.sum(den2) * d3r
		psi1 *= np.sqrt(np1/n1)
		psi2 *= np.sqrt(np2/n2)
		den1 *= np1 / n1
		den2 *= np2 / n2
	
	pot_int1, pot_int2 = dEps_dPsi(den1, den2, kappa_3_)	
	#exp(-1/2 * i * dt * V)
	psi1 *= np.exp(pot_int1 * dt)
	psi2 *= np.exp(pot_int2 * dt)
	den1 = psi1.real**2 + psi1.imag**2
	den2 = psi2.real**2 + psi2.imag**2
	if(imProp): #normalize		
		n1 = np.sum(den1) * d3r
		n2 = np.sum(den2) * d3r
		psi1 *= np.sqrt(np1/n1)
		psi2 *= np.sqrt(np2/n2)
		den1 *= np1 / n1
		den2 *= np2 / n2
	return psi1, psi2, den1, den2
	
def repeat(func, n, psi, dtfac):
	psi_tmp = np.copy(psi)
	for i in range(n):
		psi_tmp = func(psi_tmp, dtfac)
	return psi_tmp


	
def energy(psi1, psi2, npart, den1, den2):
	global dt, pot_ext, imProp,  kinprop
	#pot_tot = pot_ext * den + energy_density_interaction(den)
	pot_tot =  energy_density_interaction(den1, den2)
	psi1_k = np.fft.fftn(psi1)
	psi2_k = np.fft.fftn(psi2)
	den1_k = psi1_k.real**2 + psi1_k.imag**2
	den2_k = psi2_k.real**2 + psi2_k.imag**2
	den_k = den1_k + den2_k
	kin_tot =  (sumk2)* den_k 
	return np.sum(pot_tot)*d3r+ np.sum(kin_tot) * npart / np.sum( den_k )

output_dir = 'snapshots_time_evolution_0'; num=0
while os.path.exists(output_dir): num+=1; output_dir="snapshots_time_evolution_"+str(num)
if not os.path.exists(output_dir): os.makedirs(output_dir)
if not os.path.exists(output_dir+'/npy_files'): os.makedirs(output_dir+'/npy_files')
file_en_equil = open(output_dir + '/en_equil.dat', 'w', buffering=1)
file_en_prod = open(output_dir + '/en_prod.dat', 'w', buffering=1)
psi1, psi2 = init_psi(x, y, z)	
psi1 *= np.sqrt(np1) / np.sqrt(d3r * np.sum(psi1.real**2 + psi1.imag**2))
psi2 *= np.sqrt(np2) / np.sqrt(d3r * np.sum(psi2.real**2 + psi2.imag**2))

radius2 = x**2 + y**2 + z**2
radius = np.sqrt(radius2)


sigma = init_sigma

cin_1 = 1.
cin_1 = np.logical_or(np.abs(x) > 2.*sigma, np.abs(y) > 2.*sigma)
cin_1 = np.logical_or(cin_1,   np.abs(z) > 2.*sigma)
cin_1 = np.invert(cin_1) 

cin_2 = 1.
cin_2 = np.logical_or(np.abs(x) > 3.*sigma, np.abs(y) > 3.*sigma)
cin_2 = np.logical_or(cin_2,   np.abs(z) > 3.*sigma)
cin_2 = np.invert(cin_2) 

cin_3 = 1.
cin_3 = np.logical_or(np.abs(x) > 4.*sigma, np.abs(y) > 4.*sigma)
cin_3 = np.logical_or(cin_3,   np.abs(z) > 4.*sigma)
cin_3 = np.invert(cin_3) 



def dft_simulation(t_max,delta_t, kappa_3_):	
	global psi1, psi2, dt, kinprop
	timestep = 0
	time = 0.
	ncalc1 = np1; ncalc2 = np2
	den1 = psi1.real**2 + psi1.imag**2
	den2 = psi2.real**2 + psi2.imag**2
	p = True
	while time <= t_max:
		
		
			
		
		if(timestep % printToStdoutEvery == 0 or timestep==0):			
			if(imProp == False):				
				ncalc1 = d3r * np.sum(den1)
				ncalc2 = d3r * np.sum(den2)
				ncalc = ncalc1 + ncalc2
			else:
				ncalc = nparticles
			en = energy(psi1, psi2, ncalc, den1, den2)	
			'''
			den = den1 + den2
			meanr  = d3r * np.sum(radius  * den) / nparticles
			meanr2 = d3r * np.sum(radius2 * den) / nparticles
			'''
			n1_in = d3r * np.sum(den1 * cin_1)
			n2_in = d3r * np.sum(den2 * cin_1)
			n1_in_2 = d3r * np.sum(den1 * cin_2)
			n2_in_2 = d3r * np.sum(den2 * cin_2)
			n1_in_3 = d3r * np.sum(den1 * cin_3)
			n2_in_3 = d3r * np.sum(den2 * cin_3)
			if(imProp):				
				file_en_equil.write("%.3e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.8e\n" % (time, ncalc1, ncalc2,n1_in, n2_in, n1_in_2, n2_in_2,  n1_in_3, n2_in_3,  en))
			else:
				file_en_prod.write("%.3e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.8e\n" % (time, ncalc1, ncalc2,n1_in, n2_in, n1_in_2, n2_in_2,  n1_in_3, n2_in_3,  en))
			#print("%.3e %.5e %.5e %.8e" % (time, ncalc1, ncalc2, en))
		#log
		if(timestep % printDenEvery == 0 or timestep==0):
			timestep /= printDenEvery
			if(imProp==True):
				file=output_dir + '/npy_files/psi_equil_%i' % timestep
			else:
				file=output_dir + '/npy_files/psi_prod_%i' % timestep
			#if(imProp==False):
			#np.save(file, psi)   # use exponential notation
			np.save(file + "_xy_c1", dx[2] * np.sum(den1, axis=2))   # use exponential notation
			np.save(file + "_xy_c2", dx[2] * np.sum(den2, axis=2))   # use exponential notation
			np.save(file + "_x_slice_c1", np.swapaxes(den1,0,2)[nxyz[2] / 2][nxyz[1] / 2])   # use exponential notation
			np.save(file + "_x_slice_c2", np.swapaxes(den2,0,2)[nxyz[2] / 2][nxyz[1] / 2])   # use exponential notation
			np.save(file + "_y_slice_c1", np.swapaxes(den1,1,2)[nxyz[0] / 2][nxyz[2] / 2])   # use exponential notation
			np.save(file + "_y_slice_c2", np.swapaxes(den2,1,2)[nxyz[0] / 2][nxyz[2] / 2])   # use exponential notation
			timestep *= printDenEvery
		
		psi1, psi2, den1, den2 = T2_operator(psi1, psi2, 1., kappa_3_, den1, den2)
		if(imProp == False and activateAbsorbingWalls == True):
			psi1 *= kill_matrix 
			psi2 *= kill_matrix 
		
		time += delta_t
		timestep += 1	
		'''
		if(time >= 10.E+06 and p==True):
			p=False
			delta_t *= 5
			dt *= 5
			kinprop = np.power(kinprop, 5)
		'''
	return 0

psi1 *= np.sqrt(np1) / np.sqrt(d3r * np.sum(psi1.real**2 + psi1.imag**2))
psi2 *= np.sqrt(np2) / np.sqrt(d3r * np.sum(psi2.real**2 + psi2.imag**2))






imProp = True
dt = dt_equil
kinprop = np.exp(-1j * dt *  sumk2)
dt *= -0.5j
print(" *** EQUILIBRATION ... **** ")
dft_simulation(t_equil, deltat_equil, 0.)

'''
psi1 = np.load("../../../../equil/results/" + str(npart_fac + 1) + "/1/psi_c1.npy")
psi2 = np.load("../../../../equil/results/" + str(npart_fac + 1) + "/1/psi_c2.npy")
print("wavefunctions successfully loaded  !")
'''

psi1 *= np.sqrt(1. / 1.25)


file_sep_dist = open('sep_dist.dat', 'w', buffering=1)
den = psi2.real**2 + psi2.imag**2 + psi1.real**2 + psi1.imag**2
den = np.swapaxes(den,0,2)[nxyz[2] / 2][nxyz[1] / 2]
separation_distance = np.abs(np.argmin(np.abs(  den   -     0.025*np.max(den))) - nxyz[0] / 2) + 1
del den
print(separation_distance)
#separation_distance = (int)(2000./ dx[0])

psi1, psi2 = manipulate(psi1, psi2, x)
imProp = False
dt = dt_prod


den = psi2.real**2 + psi2.imag**2 + psi1.real**2 + psi1.imag**2
den = np.swapaxes(den,0,2)[nxyz[2] / 2][nxyz[1] / 2]
file_sep_dist.write("%d %.5e" % (separation_distance,den[nxyz[0]/2] /   np.max(den)))
print("%d %.5e" % (separation_distance,den[nxyz[0]/2] /   np.max(den)) )
file_sep_dist.close()
del den

kinprop = np.exp(-1j * dt *  sumk2 )
dt *= -0.5j
print(" *** PRODUCTION ... **** ")
tcoll = 14. * separation_distance * dx[0] / k_kick
dft_simulation(tcoll, deltat_prod, kappa_3)

print("SIM. ENDED at ", time.ctime())
print( time.time() - start_time)
'''
np.save("psi1", psi1)
np.save("psi2", psi2)
'''