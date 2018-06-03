import numpy as np
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from matplotlib.colors import ListedColormap
from numba import jit


rs = np.linspace(20, 60, 1000)
dist_from_ring = np.abs(40 - rs)

M_star = 0.5 * u.Msun.to('Mearth')
min_Mp = 3.2e20*u.kg.to('Mearth')
M_ps = np.linspace(min_Mp, 15, rs.shape[0])

M_ps = (3 * M_star) * (dist_from_ring / (5 * rs) )**3
M_ps[M_ps <= min_Mp] = min_Mp

sns.set_style('ticks')
ax1 = plt.gca(); ax2 = ax1.twinx()
ax1.semilogy(rs, M_ps)
ax1.set_xlabel(r'Perturber semi-major axis (au)'); ax1.set_ylabel(r'Perturber mass ($M_\oplus$)')
ax1.set_title('Minimum perturber mass')

ax2.set_ylim(ax1.get_ylim()); ax2.set_yscale('log')
ax2.set_yticks([min_Mp, 0.002192, 1.5, 17.1478, 317.828])
ax2.set_yticklabels([r'Inefficient Damping', 'Pluto', 'Efficient Damping', 'Neptune', 'Jupiter'])
ax2.minorticks_off(); plt.tight_layout()

plt.savefig('perturber_plot_1D', dpi=800); 
plt.show()