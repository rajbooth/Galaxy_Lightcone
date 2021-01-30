# Builds a galaxy lighcone from Gadget dark matter snapshots

from astropy import cosmology
from astropy.cosmology import Planck15 as cosmo, z_at_value
from astropy.constants import *
import scipy.integrate as integrate
from scipy.interpolate import UnivariateSpline as spl
from scipy.constants import *
import astropy.units as uni
import numpy as np
import h5py
import math
import MAS_library as MASL
import Pk_library as PKL
import time,sys,os
import datetime
import readgadget
from mpi4py import MPI

print ('All imports OK')

#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#print('My rank is ',rank)

task_ID = int(os.environ['SLURM_ARRAY_TASK_ID'])
print('Task ID =', task_ID)

dims = 2048 # dimensions of density array
verbose = False
DEBUG = 2

# Calculate z for range of comoving distances
z = np.zeros(3000)
d_c = np.zeros(3000)
H = np.zeros(3000)

for d in range(1,3000):
	d_c[d] = d
	z[d] = z_at_value(cosmo.comoving_distance, d * uni.Mpc)
	H[d]  = cosmo.H(z[d]).value
    
# create spine for quick lookup
d2z = spl(d_c,z)
z2H = spl(z, H, s=0)

def read_header(snap):
	snapshot_fname = '/cosma6/data/dp004/dc-smit4/Daemmerung/Planck2013-Npart_2048_Box_3000-Fiducial/run1/snapdir_{0:03d}/Planck2013-L3000-N2048-Fiducial_{0:03d}.0'.format(snap)
	return readgadget.header(snapshot_fname)
    
# read headers for all snapshots to get redshifts
snaps = 24
zz = np.zeros(snaps)
for snap in range(snaps):
	head = read_header(63-snap)
	zz[snap]= head.redshift

# calculate maximum redshift for each snapshot bin boundary
# as mid point between nominal redshift value for snapshot
z_max = zz[:-1:] + (zz[1::] - zz[:-1:])/2
z_max = np.append(z_max, 2*zz[-1] - z_max[-1])

# Calculate maximum comoving distance corresponding to redshift bin maxima
Dc_max = np.zeros(snaps)
for sn in range(snaps):
    Dc_max[sn] = cosmo.comoving_distance(z_max[sn]).value
Dc_min = 0    

#define galaxy datatype
gal = np.dtype([('r', np.float32),('RA', np.float32),('Dec', np.float32),('z', np.float32),('RSD', np.float32),('L', np.float32)])

# define Schechter funciton
alpha = -1.24
p_star = 0.009
def schechter(x, phi_star= p_star, a= alpha): # the luminosity function n(x) with x = L/Lstar 
	return phi_star * x**a * np.exp(-x)
    
# calculate cumulative probability distribution function for Schechter distribution
x = np.logspace(-4,1,1000)
r = np.zeros(1000)
for i in range(1000):
	r[i] = integrate.quad(schechter, x[i], np.inf, args = p_star)[0]

# create spline for looking up luminosity value corresponding to a given p value
P2L = spl(r[::-1], x[::-1], s=0)   

# calculate average particle density and L_min for normalising probability distribution
R = 3000 #Mpc
N = 2048
n = N**3 / R**3
print("Gadget particle density = ",n, " particles per Mpc")
#  we want luminosity value corresponding to this particle density
L_min = P2L(n)
print ("Minimum luminosity for Schechter distribution = ", L_min)

# Load each snapshot
#for snap in range (44,64):
def do_snap(snap):
	#load shapshot 
	print('Loading snapshot', snap)
	
	snapshot_fname = '/cosma6/data/dp004/dc-smit4/Daemmerung/Planck2013-Npart_2048_Box_3000-Fiducial/run1/snapdir_{0:03d}/Planck2013-L3000-N2048-Fiducial_{0:03d}'.format(snap)
	density_fname = '/cosma6/data/dp004/dc-boot5/Lightcone/Density/density_{0:03d}.npy'.format(snap)
	powerspec_path = '/cosma6/data/dp004/dc-boot5/Lightcone/Power_Spectrum/'
	lightcone_path = "/cosma6/data/dp004/dc-boot5/Lightcone/Galaxy_FullSky/"		
	
	#calculate density delta
	# declare the array hosting the density field
	density = np.zeros((dims, dims, dims), dtype=np.float32)

	# read relevant paramaters on the snapshot
	head     = readgadget.header(snapshot_fname)
	BoxSize  = head.boxsize/1e3 #Mpc/h
	Masses   = head.massarr*1e10 #Msun/h
	Nall     = head.nall  
	filenum  = head.filenum
	Omega_m  = head.omega_m
	Omega_l  = head.omega_l
	redshift = head.redshift
	fformat  = head.format
	Hubble   = head.Hubble
	
	Ntotal = np.sum(Nall,dtype=np.int64)

	grid     = 2048  
	MAS      = 'CIC' 
	axis 	 = 0
	do_RSD   = True
	BoxSize = 3000 #Mpc/h
	ptype = 1 #dark matter
	
	# do a loop over all files
	num = 0.0

	for i in range(filenum):
			
		# find the name of the sub-snapshot
		if filenum==1:       snapshot = snapshot_fname
		else:                snapshot = snapshot_fname + '.{0:d}'.format(i)
		if fformat=='hdf5':  snapshot = snapshot + '.hdf5'

		# find the local particles in the sub-snapshot
		head  = readgadget.header(snapshot)
		npart = head.npart

		if verbose:  print ('Sub-snapshot {0:d}, DM particles = {1:d} \n'.format(i, npart[ptype]))
		if (DEBUG>1 and i%100 == 0):  print ('Sub-snapshot {0:d}, DM particles = {1:d}, time = {2:%H:%M:%S} \n'.format(i, npart[ptype], datetime.datetime.now()))

		# read positions in Mpc/h
		pos = readgadget.read_field(snapshot, "POS ", ptype)
		
		# read velocities in km/s 
		if do_RSD:
			vel = readgadget.read_field(snapshot, "VEL ", ptype)
			
		# write galaxy data for each redshift bin into its own file
		fname = lightcone_path + 'galaxy_lightcone.snap{0:02d}'.format(snap)
		
		# open output file
		if (i == 0):
			mode = 'w'  #open file in write mode for first file in snapshot
		else:
			mode = 'r+'  #thereafter open in append mode
		fo = h5py.File(fname,mode)
			
		for oct in range(8): 
			#relocate origin to each corner of the simulation box for each octant
			orig_x = oct%2 * BoxSize
			orig_y = oct//2%2 * BoxSize
			orig_z = oct//4 * BoxSize
			
			# translate particle positions to new origin for each octant
			x = pos[::,0] - orig_x
			y = pos[::,1] - orig_y
			z = pos[::,2] - orig_z
			
			# calculate comoving radial distance, RA and Dec
			r = np.sqrt(x*x + y*y + z*z)
			dec = np.rad2deg(np.arcsin(z/r))
			ra = np.rad2deg(np.arctan2(y,x))
			
			# lookup redshift corresponding to this r
			zz = d2z(r)
			
			if do_RSD:
				# Calculate radial velocity
				vr = np.sqrt(vel[::,0]**2 + vel[::,1]**2 + vel[::,2]**2 ) * np.sign(vel[::,0] + vel[::,1] + vel[::,2])

				# Calculate RSD factor
				# Particle velocities u in internal velocity units (corresponds to km/sec if the default choice for the system of units is adopted). 
				# Peculiar velocities v are obtained by multiplying u with sqrt(a), i.e. v = u * sqrt(a). So v = u / sqrt(1+z)
				f_RSD = np.sqrt(1+zz) * vr / z2H(zz)
			else:
				f_RSD = np.zeros(len(r))

			#Check whether particle within shell max and min
			sn = 63-snap
			if (sn==0):
				F = [(r <= Dc_max[sn])]		
			else:		
				F = [(r > Dc_max[sn-1]) & (r <= Dc_max[sn])]
			f = tuple(F)
			
			# create random luminosity value for each particle
			ngal = len(pos)
			L = P2L(np.random.random(ngal) * n)	

    		# create dataset. Use f to filter only those galaxies within snapshot redshift boundaries
			g = np.array(list(zip(r[f], ra[f], dec[f], zz[f], f_RSD[f], L[f])), dtype= gal)
			
			ds_name = 'octant_{0:01d}'.format(oct)
			if (DEBUG>3): print(ds_name)
			
			# if filenum = 0 then create new datasets for each octant and set dataset atributes
			if (i == 0):
				gals = fo.create_dataset(ds_name, data = g , dtype=gal, maxshape=(None,), chunks = True) # set maxshape = None to make resizeable and chunks = True to enable chunking
				gals.attrs['max_z'] = z_max[sn]
				if (sn==0):
					gals.attrs['min_z'] = 0
				else:
					gals.attrs['min_z'] = z_max[sn-1]
				gals.attrs['snap'] = snap
				gals.attrs['octant'] = oct
				gals.attrs['alpha'] = alpha
				gals.attrs['phi_star'] = p_star
			elif (len(g)>0):
				gals = fo[ds_name]
				gals.resize(gals.shape[0] + len(g), axis=0)   
				gals[-len(g):] = g

				# end of processing for this octant
				
		fo.close()
		if (DEBUG>2): print (fname, " completed, time:", datetime.datetime.now())
		sys.stdout.flush()
		
		# compute density field. 
		MASL.MA(pos, density, BoxSize, MAS) 
		num += pos.shape[0]

	# All files read for snapshot
	if (DEBUG>0): print (fname, " completed, time:", datetime.datetime.now())
	
	# Write density field to file
	rho_avg = np.mean(density, dtype=np.float64)
	density /= rho_avg
	density -= 1.0
	density.tofile(density_fname)
	if verbose: print('Density delta written to file for snap {0:d}, mean density = {1:04f}'.format(snap, rho_avg))

	# Calculate power spectrum from density
	threads = 16

	Pk = PKL.Pk(density, BoxSize, axis, MAS, threads)
	print('Pk calculated')

	#Save power spectra	components in hdf5 file
	fname = 'powerspec_{0:03d}.npy'.format(snap)

    # open output file
	fo = h5py.File(powerspec_path + fname,'w')
	
	# create datasets
	atts = fo.create_dataset("attribs", dtype="f") # empty dataset for holding snapshot attributes
	atts.attrs['z'] = redshift
	atts.attrs['Omega_m'] = Omega_m
	atts.attrs['Omega_l'] = Omega_l
    		
    # 1D P(k)
	dset = fo.create_dataset('k1D', data = Pk.k1D)      
	dset = fo.create_dataset('Pk1D', data = Pk.Pk1D)     
	dset = fo.create_dataset('Nmodes1D', data = Pk.Nmodes1D)  
	
	# 2D P(k)
	dset = fo.create_dataset('kpar', data = Pk.kpar)      
	dset = fo.create_dataset('kper', data = Pk.kper)     
	dset = fo.create_dataset('Pk2D', data = Pk.Pk2D)  
	dset = fo.create_dataset('Nmodes2D', data = Pk.Nmodes2D)  

	# 3D P(k)
	dset = fo.create_dataset('k', data = Pk.k3D)      
	dset = fo.create_dataset('Pk0', data = Pk.Pk[:,0])     
	dset = fo.create_dataset('Pk2', data = Pk.Pk[:,1])  
	dset = fo.create_dataset('Pk4', data = Pk.Pk[:,2]) 
	dset = fo.create_dataset('Nmodes', data = Pk.Nmodes3D)
		  
	fo.close()
		
	print ('Power spectrum data written to file')
	
	# end of do_snap
	
#Run code for this snapshot
do_snap(task_ID + 41)

#do_snap(44)

print('End')