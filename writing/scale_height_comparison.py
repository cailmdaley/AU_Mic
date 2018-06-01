import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import astropy.constants as c
import astropy.units as u
import pandas as pd
from uncertainties import ufloat, umath
import uncertainties.unumpy as unp
from numba import jit

fig, ax = plt.subplots(figsize=(10,3))
rs = np.linspace(24.0, 42.3, 2)
Hs = rs*ufloat(0.031, 0.005)
ax.plot(rs, unp.nominal_values(Hs), label='This work')
plt.fill_between(rs, unp.nominal_values(Hs)+unp.std_devs(Hs), unp.nominal_values(Hs)-unp.std_devs(Hs), alpha=0.4)

rs = np.linspace(8.8, 40.3, 5)
Hs = rs*0.05
ax.plot(rs, Hs)
ax.errorbar(rs, Hs, yerr=0.25, uplims=True, label='MacGregor et al. (2013) + Schüppler et al. (2015))')

rs = np.linspace(-40.3, -8.8, 5)
Hs = np.abs(rs*0.03)
ax.plot(rs, Hs)
ax.errorbar(rs, Hs, yerr=0.25, uplims=True, label='Schneider et al. (2013) + Schüppler et al. (2015))')

ax.plot(-40, 4/2, 'o', label='Metchev et al. (2005)')


rs = np.linspace(-48.6, -11.7, 1000)rs[800]
Hs = 1.73 * (np.abs(rs)/20)**0.074

ax.fill_between(rs, 2.5, 3.5, alpha=0.4, label='Krist et al. (2005, apparent)')
ax.plot(rs, Hs, label='Krist et al. (2005, modeled)')

ax.set_ylim(0,4); 
plt.legend()
plt.show()
