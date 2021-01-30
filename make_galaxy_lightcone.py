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
    f = open('/cosma6/data/dp004/dc-boot5/Lightcone/'+ fname,'rb')
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
def schechter(x, phi_star=1.0, a=-1.24): # the luminosity function n(x) with x = L/Lstar 
    return phi_star * x**a * np.exp(-x)
    
# calculate cumulative probability distribution function for Schechter distribution
x = np.linspace(0.1,10, 1000)
r = np.zeros(1000)
for i in range(1000):
    r[i] = integrate.quad(schechter, x[i], np.inf)[0]

# create spline for looking up luminosity value correspondiong to a given p value
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
gal = np.dtype([('p',part2),('z',np.float32),('L',np.float32)])

depth = 50
 # open output file
#fname = "/lustre/scratch/astro/rb460/Lightcone/Galaxy/galaxy_lightcone"
fname = "/cosma6/data/dp004/dc-boot5/Lightcone/Galaxy/galaxy_lightcone"
fo = h5py.File(fname,'w')
for i in range(43, 64):
	p = load_snapshot(i)
	w = np.where(p['z'] < depth)[0]
	dsname = "snap_{0:2d}".format(i)
	ngal = len(p[w])
	print(dsname, ngal)
	sys.stdout.flush()
	gals = fo.create_dataset(dsname, (ngal,), dtype=gal)
	for j, par  in enumerate(p[w]):
		# calculate comoving radial distance
		r = np.sqrt(par['x']**2  + par['y']**2 + par['z']**2)
        # lookup redshift corresponding to this r
		z = d2z(r)
        # create random luminosity value for this particle
		L = P2L(np.random.random() * n)
        # assign values in gal structure
		gals[j] = [(par,z,L)]
	print ("Snap =", i, " Time:", datetime.datetime.now())
fo.close()
    