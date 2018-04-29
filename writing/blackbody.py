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
    
bounds = np.array([15.92,16]) * u.au
blackbody_fluxes = disk_SED(bounds, 1)
star_flux = blackbody_flux(3500*u.K ,wavs) * 4 * np.pi * (0.85 * u.Rsun.to('m'))**2 *3
combined_fluxes = star_flux + blackbody_fluxes.sum(axis=0)
distance_coeff = 1 / (4 * np.pi * (9.91*u.pc.to('m'))**2)
wavs

for flux in blackbody_fluxes:
    test = plt.loglog(wavs, flux * distance_coeff, 
        ls=':', c='xkcd:cerulean', label='Disk Contribution') 
plt.loglog(wavs, star_flux * distance_coeff, 
    ls='--', c='xkcd:crimson', label='Stellar Conribution')
plt.loglog(wavs, combined_fluxes * distance_coeff, 
    ls='-', color='black', alpha=0.6, label='Observed SED')
index = (np.abs(wavs-1.4*u.mm)).argmin()
plt.plot(wavs[index], combined_fluxes[index]*distance_coeff, 
    'r*', markersize=15, label=r'ALMA')
plt.ylim(1e-3, 1e1); plt.xlim(1e-1, 9e3)
plt.xlabel(r'Emitting Wavelength ($\mu$m)'); plt.ylabel(r'Flux Density (Jy)')
plt.legend(); sns.despine(); plt.tight_layout()
plt.savefig('figures/disk_SED.pdf', dpi=1200); 
plt.show()
