import numpy as np
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from matplotlib.colors import ListedColormap
from numba import jit

rs = np.linspace(32, 48, 1000)
dist_from_ring = np.abs(40 - rs)

M_star = 0.5 * u.Msun.to('Mearth')
min_Mp = 3.2e20*u.kg.to('Mearth')
M_ps = np.linspace(min_Mp, 15, rs.shape[0])
es = np.linspace(0, 1, rs.shape[0])
np.round(M_ps,2)
dist_from_ring

mass_map = np.zeros((rs.shape[0], es.shape[0]))
mass_map[:] = np.nan

def loop():
    for i, e in enumerate(es):
        for M_p in M_ps[::-1]:
            r_Hs = rs * (1 - e) * (1/3 * M_p/M_star)**(1/3)
            if r_Hs.all() !=0:
                indices = np.abs(5*r_Hs - dist_from_ring).argsort()[:4]
            if M_p == M_ps[-1]:
                mass_map[i, indices.min():indices.max()] = min_Mp
            mass_map[i, indices.min()] = M_p
            mass_map[i, indices.max()] = M_p
loop()

sns.set_style('ticks')
cmap = sns.color_palette('Blues_d', 10)
cmap.reverse()
xticks = np.linspace(32,48,9, dtype=int)
yticks = np.linspace(0,1,3)
ax = sns.heatmap(mass_map[::-1,:], cmap=cmap, cbar_kws={'pad' : 0.01}, xticklabels=xticks, yticklabels=yticks)

cbar = ax.collections[0].colorbar
cbar.set_ticks([min_Mp, 3, 6, 9, 12, 15])
cbar.set_ticklabels(['{:.0g}'.format(min_Mp), 3, 6, 9, 12, 15])
cbar.set_label(r'Perturber mass (M$_\oplus$)')
cbar.ax.tick_params(axis='y', right='false', pad=0.1)

ax.set_xticks(np.linspace(0,1,9) * ax.get_xlim()[1])
ax.set_yticks(np.linspace(0,1,3)[::-1] * ax.get_ylim()[0])
ax.set_xlabel(r' Perturber separation (au)'); ax.set_ylabel(r'Perturber Eccentricity')

for direction, spine in ax.spines.items():
    if direction in ['left', 'bottom']:
        spine.set_visible(True)
plt.savefig('perturber_plot', dpi=800); plt.show()
