import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import astropy.constants as c
import astropy.units as u
from astropy.visualization import quantity_support; quantity_support()
from astrocail import mcmc, utils

wavs = np.logspace(-0.9, 4, 1000)*u.micron 
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
    return f_nu.to(u.Jy)
    
def disk_SED(disk_bounds, n_annuli):
    # create array of relevant distances
    radii = np.linspace(*disk_bounds, 2 * n_annuli + 1)
    
    annulus_Ts = T_dust(radii[1::2]) # get temperature at middle of each annulus
    annulus_Is = np.array([blackbody_flux(T, wavs) for T in annulus_Ts])
    annulus_areas = np.pi * (radii[2::2]**2 - radii[0:-1:2]**2).to('m**2')
    annulus_fs = np.array([annulus_Is[i,:] * annulus_areas[i] 
        for i in range(n_annuli)])
    
    return annulus_fs * u.Jy
    
bounds = np.array([41.9, 42]) * u.au
blackbody_fluxes = disk_SED(bounds, 1)
star_flux = blackbody_flux(3500*u.K ,wavs) * 4 * np.pi * (0.85 * u.Rsun.to('m'))**2 
combined_fluxes = star_flux + blackbody_fluxes.sum(axis=0)
distance_coeff = 1 / (4 * np.pi * 9.91*u.pc.to('m')**2)

for flux in blackbody_fluxes:
    test = plt.loglog(wavs, flux * distance_coeff, 
        ls=':', c='xkcd:cerulean', label='Disk Contribution') 
plt.loglog(wavs, star_flux * distance_coeff, 
    ls='--', c='xkcd:crimson', label='Stellar Conribution')
plt.loglog(wavs, combined_fluxes * distance_coeff, 
    ls='-', color='black', label='Observed SED')
plt.ylim(1e-2, 1e2); plt.xlim(1e-1, 9e3)
plt.xlabel(r'Emitting Wavelength ($\mu$m)'); plt.ylabel(r'Flux Density (Jy)')
plt.legend(); sns.despine(); plt.tight_layout()
plt.axvline(1.4e3, c='r', ls=':', lw=1.5, ymax=0.27)
plt.savefig('figures/disk_SED.pdf', dpi=1200); plt.show()
