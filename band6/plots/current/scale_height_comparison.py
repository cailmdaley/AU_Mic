import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

sns.set_style('whitegrid')
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,3), sharey=True)
c = sns.color_palette('deep', 6)
fontsize=14

# Optical & Near-IR

## Metchev
rs = np.linspace(20, 40, 4)
HWHMs = np.linspace(2.5,4.2,4)/2
ax1.plot(rs, HWHMs,  ls=':', label='Metchev et al. (2005)', color=c[2])
ax1.annotate(s='1', xy=(rs[-1]-0.5, HWHMs[-1]-0), color=c[2], fontsize=fontsize)

## Schuppler
rs = np.linspace(37.5,42.5, 4)
HWHMs = rs * 0.1/2
ax1.plot(rs, HWHMs,  ls='-', color=c[0],
    label='Schuppler et al. (2015), short wavelengths')
ax1.annotate(s='2', xy=(rs[-1]+2, HWHMs[-1]-0.05), color=c[0], fontsize=fontsize)
 
## Krist, apparent
rs = np.linspace(16.5, 50, 100)
HWHMs = np.linspace(2.5,3.5,100) / 2
ax1.plot(rs, HWHMs,  ls=':', label='Krist et al. (2005), apparent', color=c[1])
ax1.annotate(s='3', xy=(rs[-1]-1, HWHMs[-1]+0.05), color=c[1], fontsize=fontsize)

## Schneider
rs = np.linspace(20, 50, 4)
HWHMs = rs*0.03
ax1.plot(rs, HWHMs,  ls=':', color=c[3],
    label='Schneider et al. (2014) + Schuppler et al. (2015)')
ax1.annotate(s='4', xy=(rs[-1]-1, HWHMs[-1]), color=c[3], fontsize=fontsize)
    
## Krist, modeled
rs = np.linspace(11.7, 48.8, 100)
HWHMs = 1.73*(rs/20)**0.074 / 2
ax1.plot(rs, HWHMs,  ls='-', color=c[5], 
    label='Krist et al. (2005), modeled')
ax1.annotate(s='5', xy=(rs[-1], HWHMs[-1]+0.08), color=c[5], fontsize=fontsize)

# Millimeter

## Macgregor
rs = np.linspace(8.8,40.3, 4)
HWHMs = rs*0.05
ax2.errorbar(rs, HWHMs, yerr=0.1, ls=':', uplims=True, 
    label='Macgregor et al. (2013) + Schuppler et al. (2015)', color=c[2])
ax2.annotate(s='6', xy=(rs[-1]+1, HWHMs[-1]), color=c[2], fontsize=fontsize)

# This work
rs = np.linspace(24,42.3, 3)
HWHMs = rs*0.031*2.355/2
ax2.plot(rs, HWHMs, ls='-', label='This work', color=c[0])
ax2.fill_between(rs, rs*0.026*2.35/2, rs*0.036*2.35/2, 
    alpha=0.3, color=c[0])
ax2.annotate(s='7', xy=(rs[-1]+1, HWHMs[-1]), color=c[0], fontsize=fontsize)
    
rs = np.linspace(37.5,42.5, 4)
HWHMs = rs*0.03/2
ax2.plot(rs, HWHMs,  ls='-', color=c[1],
    label='Schuppler et al. (2015), reference model')
ax2.annotate(s='8', xy=(rs[-1]+1, HWHMs[-1]), color=c[1], fontsize=fontsize)
    
rs = np.linspace(37.5,42.5, 4)
HWHMs = rs*0.015/2
ax2.plot(rs, HWHMs,  ls='-', color=c[3],
    label='Schuppler et al. (2015), long wavelengths')
ax2.annotate(s='9', xy=(rs[-1]+1, HWHMs[-1]), color=c[3], fontsize=fontsize)
    
ax1.set_xlim(0,50); 
ax1.invert_xaxis()
ax2.set_xlim(0,50)
ax1.set_ylim(0, 2.3)
ax2.tick_params(axis='y', labelright=True)

ax1.set_ylabel('HWHM (au)')
fig.text(0.5125, 0.0175, 'Radius (au)', ha='center', fontsize=11)
ax1.set_title('Optical & Near-IR')
ax2.set_title('Millimeter')

plt.subplots_adjust(wspace=0, bottom=0.15)

plt.savefig('scale_height_comparison.png', dpi=1000)
plt.savefig('../../../writing/figures/scale_height_comparison.png', dpi=1000)
plt.show()
