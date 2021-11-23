
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
	# return np.power( -alpha / ( beta * gamma) , 1. / (gamma - 1) )
	return (25*np.pi/1024.) * np.power(a12 + np.sqrt(a22), 2) / (np.power(a22, 3./2) * np.power(1. + np.sqrt(a22), 4))




def init_psi(x, y, z, init_sigma):
	psi1 = np.exp(-0.5*(x**2/init_sigma**2 + y**2/init_sigma **
	              2 + z**2/init_sigma**2)) + 0.j  # Gaussian
	psi2 = np.exp(-0.5*(x**2/init_sigma**2 + y**2/init_sigma **
	              2 + z**2/init_sigma**2)) + 0.j  # Gaussian
	return psi1, psi2


def manipulate(psi1, psi2, x, x_distance):
    """
    
    """
	psi1_l = np.roll(psi1, -separation_distance, axis=0)
	psi2_l = np.roll(psi2, -separation_distance, axis=0)

	psi1_r = np.roll(psi1, separation_distance, axis=0)
	psi2_r = np.roll(psi2, separation_distance, axis=0)

	psi1 = psi1_l * np.exp(1.j * k_kick * x) + psi1_r * np.exp(-1.j * k_kick * x)
	psi2 = psi2_l * np.exp(1.j * k_kick * x) + psi2_r * np.exp(-1.j * k_kick * x)

	return psi1, psi2




def init_kill_matrix(x, y, z, kill_fac = 0.95):
	# return 3. - np.exp( - ((x - L_ok[0]/2.)/(L[0]/2. - L_ok[0]/2))**20 )  - np.exp( - ((y - L_ok[1]/2.)/(L[1]/2. - L_ok[1]/2))**2 ) -np.exp(-  ((z - L_ok[2]/2.)/(L[2]/2. - L_ok[2]/2))**2 )
	# return 1 #. if you don't want to kill wavefunction
	matrix = 1.
	if(N_DIM == 1):
		matrix = np.logical_or(np.abs(x) > L_ok[0]/2.)
	if(N_DIM == 2):
		matrix = np.logical_or(np.abs(x) > L_ok[0]/2., np.abs(y) > L_ok[1]/2.)
	if(N_DIM == 3):
		matrix = np.logical_or(np.abs(x) > L_ok[0]/2., np.abs(y) > L_ok[1]/2.)
		matrix = np.logical_or(matrix,   np.abs(z) > L_ok[2]/2.)
	# this matrix is 1 inside box, and kill_fac outside
	return np.invert(matrix) + matrix*kill_fac