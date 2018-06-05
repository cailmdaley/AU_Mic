import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import astropy.constants as c
import astropy.units as u
import pandas as pd
from uncertainties import ufloat, umath
import uncertainties.unumpy as unp
from numba import jit

M_dust = unp.uarray( 
    [[10**(-7.54), 0,  10**(-7.54)           ], 
     [10**(-7.54),  0, 10**(-7.54)           ]],
    [[(10**(-7.54)  - 10**(-7.546)) / 2, 0, 0],
     [(10**(-7.534) - 10**(-7.54))  / 2, 0, 0]]
    ) * u.Msun.to('g')
r_in =  unp.uarray( 
    [[24.0, 8.8, 5.2], 
     [24.0, 8.8, 5.2]],
    [[0.9, 1.,   0  ],
     [0.6, 11.,  0  ]]
    ) * u.au.to('cm')
r_out = unp.uarray( 
    [[42.3, 40.3, 40.6], 
     [42.3, 40.3, 40.6]], 
    [[0.5,  0.4,  0   ],
     [0.4,  0.4,  0   ]]
    ) * u.au.to('cm')
p = unp.uarray( 
    [[0.8, 2.3, 2.8], 
     [0.8, 2.3, 2.8]], 
    [[0.4, 0.31, 0],
     [0.5, 0.21, 0]]
    )
macgregor_flux_values = unp.uarray([7.14, 7.14], [0.25, 0.12])
r_in_errs = np.array([ [0.9, 0.6], [1, 11], [0,0] ]) * u.au.to('cm')
r_out_errs = np.array([ [0.5, 0.4], [0.4, 0.4], [0,0] ]) * u.au.to('cm')
r = np.array([np.linspace(
    r_in[0,i].n  - r_in[0,i].std_dev, 
    r_out[1,i].n + r_out[1,i].std_dev, 
    1000) for i in range(3)])
r_au = r * u.cm.to('au')

@jit
def T_dust(radii):
    L_star = 0.09 * u.Lsun
    Ts = ( L_star / (16 * np.pi * radii**2 * c.sigma_sb))**(1/4)
    return Ts.to(u.K)
    
@jit
def blackbody_flux(Ts, wav, area=1, distance=1):
    nu = c.c / wav
    coef = 2 * c.h * nu**3 / (c.c**2) 
    distribution =  1 / ( np.exp(c.h * nu / (c.k_B * Ts)) - 1)
    spectral_radiance = coef * distribution
    f_nu = spectral_radiance * ( area / (4 * np.pi * distance**2 ) )
    return f_nu.to(u.Jy)

def macgregor_surface_density(n_annuli):
    radii = np.array([
        np.linspace(r_in[i,1].n, r_out[i,1].n, 2 * n_annuli + 1) 
        for i in range(r_in.shape[0])]) * u.cm
    areas = np.array(
        [np.pi * (radii[i, 2::2]**2 - radii[i, 0:-1:2]**2) 
        for i in range(radii.shape[0])])*u.cm**2
    
    blackbody_fluxes = np.array([
        blackbody_flux(T_dust(radii[i, 1::2]), 
        wav=1.28*u.mm, 
        area=areas[i], 
        distance=9.91*u.pc).to_value('mJy')
        for i in range(radii.shape[0])])
        
    macgregor_fluxes = np.array([radii[i,1::2].to_value('au')**p[i,1] 
        for i in range(radii.shape[0])] )
    norm = (  macgregor_flux_values / macgregor_fluxes.sum(axis=1) ) # mJy
    macgregor_fluxes *= norm.reshape(2,1) 
    
    m_dust_grain = (np.pi * (1.28*u.mm)**2 
        / ( 2.3*u.cm**2/u.g ) ).to_value('g')
    surface_densities = m_dust_grain * macgregor_fluxes / blackbody_fluxes 
    return surface_densities
    
@jit
def our_surface_density():
    sigma_c= M_dust * (p + 2) / (2 * np.pi * (r_out**(p+2) - r_in**(p+2) ))
    sigma = sigma_c.reshape(2,3,1) * np.array([r]*2)**p.reshape(2,3,1)
    return sigma
    
sigma = our_surface_density()
sigma[:, 1] = macgregor_surface_density(1000)

fig, ax = plt.subplots()
labels = ['This work', 'MacGregor et al. (2013)', r'$r_{min} \leq 15$ au best fit']
ls = ['-', '--', ':']
colors = [sns.color_palette()[0], 'crimson', 'k']
for i in range(1, -1, -1):
    ax.semilogy(r_au[i], unp.nominal_values(sigma[0,i]), 
        label=labels[i], ls=ls[i], color=colors[i])
        
    inner = np.where(r_au[i] < r_au[i,0]  
        + unp.std_devs(r_in[:,i]).sum() * u.cm.to('au'))
    outer = np.where(r_au[i] > r_au[i,-1] 
        - unp.std_devs(r_out[:,i]).sum() * u.cm.to('au'))
    ax.fill_between(r_au[i][inner], 
        (unp.nominal_values(sigma[0,i,inner]) - unp.std_devs(sigma[0,i,inner])).flatten(), 
        (unp.nominal_values(sigma[1,i,inner]) + unp.std_devs(sigma[1,i,inner])).flatten(), 
        alpha=0.15, color=colors[i])
    ax.fill_between(r_au[i][outer], 
        (unp.nominal_values(sigma[0,i,outer]) - unp.std_devs(sigma[0,i,outer])).flatten(), 
        (unp.nominal_values(sigma[1,i,outer]) + unp.std_devs(sigma[1,i,outer])).flatten(), 
        alpha=0.15, color=colors[i])
    ax.fill_between(r_au[i][inner[0][-1]:outer[0][0]],
        unp.nominal_values(sigma[0, i, inner[0][-1]:outer[0][0]]) 
            - unp.std_devs(sigma[0, i, inner[0][-1]:outer[0][0]]), 
        unp.nominal_values(sigma[1, i, inner[0][-1]:outer[0][0]]) 
            + unp.std_devs(sigma[1, i, inner[0][-1]:outer[0][0]]), 
        alpha=0.6, color=colors[i])
            
    start = (np.abs(r_au[i]-r_in[0, i].n*u.cm.to('au'))).argmin()
    ax.plot(r_au[i, start], sigma[0, i, start].n, 'o', c=colors[i], markersize=4)
    
    end = (np.abs(r_au[i]-r_out[0, i].n*u.cm.to('au'))).argmin()
    ax.plot(r_au[i, end], sigma[0, i, end].n, 'o', c=colors[i], markersize=4)
ax.semilogy(r_au[2], unp.nominal_values(sigma[0,2]), label=labels[2], ls=ls[2], c=colors[2])
plt.xlabel('$r$ (au)'); plt.ylabel(r'Surface density (g/cm$^2$)')
plt.legend(loc=4); sns.despine(); plt.savefig('../../../writing/figures/surface_density', dpi=1000)
plt.show()

# ax1.fill_between(ra_range, SE_amps-rms, SE_amps+rms, alpha=0.4)
