from scipy.ndimage.interpolation import rotate
from scipy.optimize import curve_fit
from astropy.io import fits
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# image = 'aumic_18aug_usermask_natural'
# image = 'aumic_24jun_usermask_natural'
image = 'aumic_usermask_natural'

# Read the header from the observed FITS continuum image:
head = fits.getheader(image + ".fits")
print(head)
# Read in images and rotate so that disk is horizontal
angleSE = 129.5 - 90
im = rotate(fits.getdata(image + ".fits").squeeze(),
            angleSE, reshape=False) * 1e6

# Generate x and y axes: offset position in arcsec
nx = head['NAXIS1']
xpix = head['CRPIX1']
xdelt = head['CDELT1']

ny = head['NAXIS2']
ypix = head['CRPIX2']
ydelt = head['CDELT2']

# Convert from degrees to arcsecs
ra = ((np.arange(nx) - xpix + 1) * xdelt) * 3600
dec = ((np.arange(ny) - ypix + 1) * ydelt) * 3600
xaxis = ra[257:]
print(xaxis[0])
print(dec[257])

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
sns.set_style('whitegrid')

SE_midplane_intensity = np.zeros(len(xaxis))
NW_midplane_intensity = np.zeros(len(xaxis))
for i in range(len(xaxis)):
    NW_midplane_intensity[i] = np.mean(im[255:258, 257+i])
    SE_midplane_intensity[i] = np.mean(im[255:258, 256-i])

ax1.set_title("Midplane intensity")
ax1.set_ylabel(r"Surface brightness ($\mu$Jy/beam)...")
ax1.set_xlim(0, 6)
ax1.set_ylim(-50, 350)
ax1.plot(xaxis, NW_midplane_intensity, label='NW')
ax1.plot(xaxis, SE_midplane_intensity, label='SE')
legend1 = ax1.legend()


# Define model Gaussian function:
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

p0 = [100., 0., 0.3]
bounds = ([-50, -0.5, 0], [400, 0.5, 1])
NW_amps, NW_mus, NW_sigmas = np.zeros(len(xaxis)), np.zeros(len(xaxis)), np.zeros(len(xaxis))
SE_amps, SE_mus, SE_sigmas = np.zeros(len(xaxis)), np.zeros(len(xaxis)), np.zeros(len(xaxis))

for i in range(len(xaxis)):
    (NW_amps[i], NW_mus[i], NW_sigmas[i]), NW_var_matrix = curve_fit(gauss, dec[256-35:256+35], im[256-35:256+35,257+i], p0=p0)
    (SE_amps[i], SE_mus[i], SE_sigmas[i]), SE_var_matrix = curve_fit(gauss, dec[256-35:256+35], im[256-35:256+35,256-i], p0=p0)

# ax2.set_title("Spine elevation from midplane")
ax2.set_xlabel("Projected seperation from star (\")")
ax2.set_ylabel("Spine elevation from midplane")
ax2.set_xlim(0, 6)
ax2.plot(xaxis, NW_mus, label='NW')
ax2.plot(xaxis, SW_mus, label='SE')
legend2 = ax2.legend()

ax3.set_xlabel("Projected seperation from star (\")")
ax2.set_ylabel("Disk FWHM")
ax2.set_xlim(0, 6)
ax2.plot(xaxis, NW_sigmas, label='NW')
ax2.plot(xaxis, SW_sigmas, label='SE')
legend3 = ax3.legend()

fig.savefig("boccaletti_plots.png")
fig.show()
