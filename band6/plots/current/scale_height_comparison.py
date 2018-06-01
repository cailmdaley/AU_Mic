import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

sns.set_style('whitegrid')
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,3), sharey=True)
optical = sns.color_palette()[0] 
mm = 'crimson'

rs = np.linspace(24,42.3, 3)
Hs = rs*0.031
ax2.plot(rs, rs*0.031*2.35/2, ls='-', label='This work')
ax2.fill_between(rs, rs*0.026*2.35/2, rs*0.036*2.35/2, 
    alpha=0.3)

rs = np.linspace(8.8,40.3, 4)
ax2.errorbar(rs, rs*0.05, yerr=0.1, ls='--', uplims=True, 
    label='Macgregor et al. (2013) + Schüppler et al. (2015)')
    
rs = np.linspace(37.5,42.5, 4)
ax2.plot(rs, rs*0.03/2,  ls=':',
    label='Schüppler et al. (2015), reference model')
    
ax2.plot(rs, rs*0.015/2,  ls='-.',
    label='Schüppler et al. (2015), long wavelengths')
    
ax1.plot(rs, rs*0.03,  ls='-',
    label='Schneider et al. (2014) + Schüppler et al. (2015)')
    
ax1.plot(rs, rs*0.1/2,  ls='--',
    label='Schüppler et al. (2015), short wavelengths')
    
rs = np.linspace(11.7, 48.8, 100)
ax1.plot(rs, 1.73*(rs/20)**0.074 / 2,  ls=':', 
    label='Krist et al. (2005), modeled')
    
rs = np.linspace(16.5, 50, 100)
ax1.plot(rs, np.linspace(2.5,3.5,100) / 2,  ls='-.', label='Krist et al. (2005), apparent')

rs = np.linspace(20, 40, 4)
ax1.plot(rs, np.linspace(2.5,4.2,4)/2,  ls=':', label='Metchev et al. (2005)')

ax1.set_xlim(0,50); ax2.set_xlim(0,50)
ax1.legend(bbox_to_anchor=(0, -0.05, 0.9, -0.05)); ax2.legend(bbox_to_anchor=(0, -0.05, 1, -0.05))
ax2.tick_params(axis='y', labelright=True)
ax1.invert_xaxis()
plt.subplots_adjust(wspace=0)
plt.show()
