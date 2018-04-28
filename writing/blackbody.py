import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import astropy.constants as c
import astropy.units as u
from astropy.visualization import quantity_support; quantity_support()
from astrocail import mcmc, utils

wavs = np.logspace(-0.9, 3.5, 100)*u.micron 
# wavs = np.linspace(0, 30, 5000)*u.micron 
def T_dust(radii):
    L_star = 0.09 * u.Lsun
    Ts = ( L_star / (16 * np.pi * radii**2 * c.sigma_sb))**(1/4)
    return Ts.to(u.K)
    
T_dust(10 * u.au)

def blackbody_flux(temp, wavs):
    nu = c.c / wavs
    coef = 2 * c.h * nu**3 / (c.c**2) 
    dist =  1 / ( np.exp(c.h * nu / (c.k_B * temp)) - 1)
    f_nu = coef * dist
    nu_f_nu  = f_nu * nu
    return nu_f_nu.to(u.W / u.m**2)
    
def disk_SED(disk_bounds, n_annuli):
    # create array of relevant distances
    radii = np.linspace(*disk_bounds, 2 * n_annuli + 1)
    
    annulus_Ts = T_dust(radii[1::2]) # get temperature at middle of each annulus
    annulus_Is = np.array([blackbody_flux(T, wavs) for T in annulus_Ts])
    annulus_areas = np.pi * (radii[2::2]**2 - radii[0:-1:2]**2)
    annulus_fs = np.array([annulus_Is[i,:] * annulus_areas[i] 
        for i in range(n_annuli)])
    
    return annulus_fs * u.W / u.m**2
    
bounds = np.array([10, 42]) * u.au
blackbodies = disk_SED(bounds, 5)

star = blackbody_flux(3500*u.K ,wavs)
comp = star + blackbodies.sum(axis=0)



plt.loglog(wavs, star, ls='--', label='Stellar Conribution')
for blackbody in blackbodies:
    test = plt.loglog(wavs, blackbody, ls=':') 
plt.loglog(wavs, comp, ls='-', label='Observed SED')
# plt.loglog(wavs, blackbody_flux(50*u.K ,wavs))
# star = blackbody(3500, wavs)
# disk = blackbody(500, wavs)
# plt.loglog(wavs, star, label='star')
# plt.loglog(wavs, star + disk, label='comp')
plt.ylim(1e-2, 1e8)
plt.title(r'$f_\lambda$ as a Function of Wavelength')
plt.xlabel('Emitting Wavelength (m)'); plt.ylabel(r'Flux ($\nu F_\nu$)')
plt.legend(); sns.despine(); plt.tight_layout()
plt.savefig('Homework1_blackbody.png'); plt.show()
