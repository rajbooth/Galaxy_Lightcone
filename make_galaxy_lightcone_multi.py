import numpy as np
import h5py
import math

from astropy import cosmology
from astropy.cosmology import Planck15 as cosmo, z_at_value
from astropy.constants import *
import scipy.integrate as integrate
from scipy.interpolate import UnivariateSpline as spl
from scipy.constants import *
import astropy.units as uni
import datetime
import sys

header3 = np.dtype([('npartFile',np.int32,6),('mass',np.double,6),('time',np.double),('redshift',np.double),
                    ('flag_sfr',np.int32), ('flag_feedback',np.int32),('npartTotal',np.int32,6),
                    ('flag_cooling',np.int32), ('num_files',np.int32),
                    ('BoxSize',np.double),('Omega0',np.double),('OmegaLambda',np.double),('HubbleParam',np.double),  
                     ('w0',np.double), ('wa', np.double), ('flag_stellarage', np.int32), ('flag_metals', np.int32),
                     ('npartTotalHighWord', np.uint32,6),
                    ('flag_entropy_instead_u', np.int32),('flag_doubleprecision',np.int32), ('flag_ic_info', np.int32),
                    ('lpt_scalingfactor',np.float32),
       ('fill',np.int8,32)])
       
part2 = np.dtype([('x',np.float32),('y',np.float32),('z',np.float32)])

def load_snapshot(snap):
	fname = 'Lightcone_{:03d}'.format(snap)  
    #f = open('/lustre/scratch/astro/rb460/Lightcone/'+ fname,'rb')
    #f = open('/cosma6/data/dp004/dc-boot5/Lightcone/'+ fname,'rb')
	f = open('/cosma6/data/dp004/dc-boot5/Lightcone/DM_FullSky/'+ fname,'rb')
	#f = open('/cosma6/data/dp004/dc-boot5/Lightcone/DM_Octant/'+ fname,'rb')
	dummy = np.fromfile(f, dtype=np.int32,count=1)
	headdata=np.fromfile(f, dtype=header3,count=1)
	dummy = np.fromfile(f, dtype=np.int32,count=1)
	pcount = headdata['npartTotal'][0][1]
	print (pcount)
	dummy = np.fromfile(f, dtype=np.int32,count=1)
	parts = np.fromfile(f,dtype=part2, count = pcount)
	f.close()
	return parts

# Calculate z for range of comoving distances
z = np.zeros(3000)
d_c = np.zeros(3000)
for d in range(1,3000):
    d_c[d] = d
    z[d] = z_at_value(cosmo.comoving_distance, d * uni.Mpc)
    
# create spine for quick lookup
d2z = spl(d_c,z)

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

#define galaxy datatype
#gal = np.dtype([('p',part2),('z',np.float32),('L',np.float32)]) #old cartesian version
gal = np.dtype([('r', np.float32),('RA', np.float32),('Dec', np.float32),('z', np.float32),('RSD', np.float32),('L', np.float32)])

for snap in range(53, 64):
	p = load_snapshot(snap)
	ngal = len(p)
	# we want about 1Gbyte per file. Each gal dtype is 5 x float32 = 20 bytes
	nbins = int(np.ceil(ngal * 20 / 1e9))
	print ('Filesize = ', ngal*20 / 1000000, 'Mbytes, bins = ', nbins)
	print(datetime.datetime.now().time())
	sys.stdout.flush()

	# calculate comoving radial distance, ra, and dec    
	r = np.sqrt(p['x']**2  + p['y']**2 + p['z']**2)
	dec = np.rad2deg(np.arccos(p['z']/r))
	ra = np.rad2deg(np.arcsin(p['y']/np.sqrt(p['x']**2  + p['y']**2)))

	# lookup redshift corresponding to this r
	z = d2z(r)
	max_z = np.amax(z)
	min_z = np.amin(z)
	zrange = max_z - min_z
	zinc = zrange /nbins
	print ('Max z = {0:03f}, min z = {1:03f}, z range = {2:03f}, z inc = {3:03f}'.format(max_z, min_z, zrange, zinc))
	print(datetime.datetime.now().time())
	sys.stdout.flush()
		
	# assign each particle to appropriate redshift shell (bin)
	pbin = np.floor((z - min_z)/zinc)
	#print(datetime.datetime.now().time())
	
	# create random luminosity value for each particle
	L = P2L(np.random.random(ngal) * n)	
	#print(datetime.datetime.now().time())
    	
    # write galaxy data for each bin into its own output file
	fpath = "/cosma6/data/dp004/dc-boot5/Lightcone/Galaxy_FullSky/galaxy_lightcone"
	for b in range(nbins):
		fname = fpath + '.snap{0:02d}.shell{1:02d}'.format(snap, b)
    	# open output file
		fo = h5py.File(fname,'w')
    	# create dataset
		g = np.array(list(zip(p[pbin==b],z[pbin==b], L[pbin==b])), dtype= gal)
		gals = fo.create_dataset('gal_data', data = g , dtype=gal)
		gals.attrs['max_z'] = max_z
		gals.attrs['min_z'] = min_z
		gals.attrs['snap'] = snap
		gals.attrs['shells'] = nbins
		gals.attrs['shellnum'] = b
		gals.attrs['alpha'] = alpha
		gals.attrs['phi_star'] = p_star
		fo.close()
		print (fname, " completed, time:", datetime.datetime.now())
		sys.stdout.flush()
		

    