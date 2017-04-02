from scipy.ndimage.interpolation import rotate
from scipy.optimize import curve_fit
from astropy.io import fits
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# image = '../cleans/aumic_composite_usermask_natural'
image = '../cleans/aumic_ctrpix_test_usermask_natural'

# Read the header from the observed FITS continuum image:
head = fits.getheader(image + ".fits")
# Read in images and rotate so that disk is horizontal
angleSE = 129.5 - 90
im = rotate(fits.getdata(image + ".fits").squeeze(),
            angleSE, reshape=False) * 1e6

# Generate x and y axes: offset position in arcsec
nx = head['NAXIS1']
xpix = head['CRPIX1']
xdelt = head['CDELT1']
# print(head['CRVAL1'])
# print(head['CRPIX1'])

ny = head['NAXIS2']
ypix = head['CRPIX2']
ydelt = head['CDELT2']

# Convert from degrees to arcsecs
ra = ((np.arange(nx) - xpix + 1) * xdelt) * 3600
dec = ((np.arange(ny) - ypix + 1) * ydelt) * 3600
ra[256]
xaxis = ra[1:256][::-1]
# print(ra[0],ra[-1])
# print(dec[256])

SE_midplane_intensity = np.zeros(len(xaxis))
NW_midplane_intensity = np.zeros(len(xaxis))
for i in range(len(xaxis)):
    NW_midplane_intensity[i] = np.mean(im[255:258, 256 + i])
    SE_midplane_intensity[i] = np.mean(im[255:258, 256 - i])


# Define model Gaussian function:
def gauss(x, *p):
    A, mu, sigma = p
    return A * np.exp(-(x - mu)**2 / (2. * sigma**2))

p0 = [100., 0., 0.3]
bounds = ([-50, -0.5, 0], [400, 0.5, 1])
NW_amps, NW_mus, NW_sigmas = np.zeros(
    len(xaxis)), np.zeros(len(xaxis)), np.zeros(len(xaxis))
SE_amps, SE_mus, SE_sigmas = np.zeros(
    len(xaxis)), np.zeros(len(xaxis)), np.zeros(len(xaxis))

for i in range(len(xaxis)):
    (SE_amps[i], SE_mus[i], SE_sigmas[i]), SE_var_matrix = curve_fit(gauss, dec[
        256 - 35:256 + 35], im[256 - 35:256 + 35, 256 - i], p0=p0, bounds=bounds)
    (NW_amps[i], NW_mus[i], NW_sigmas[i]), NW_var_matrix = curve_fit(gauss, dec[256-35:256+35], im[256-35:256+35,256+i], p0=p0, bounds=bounds)

# Subtract beam FWHM
bmin = head['bmin'] * 3600.
bmaj = head['bmaj'] * 3600.
bpa = head['bpa']
b_FWHM = 2* bmin * bmaj / \
    np.sqrt((bmin * np.cos(bpa + angleSE - 90))**2 +
            (bmaj * np.sin(bpa + angleSE - 90))**2)
print(b_FWHM)

SE_disk_FWHM = np.sqrt((2.35482004503*SE_sigmas)**2 - b_FWHM**2)
NW_disk_FWHM = np.sqrt((2.35482004503*NW_sigmas)**2 - b_FWHM**2)

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(8, 12), sharex=True)
sns.set_style('whitegrid')

ax1.set_title("Midplane intensity")
ax1.set_ylabel(r"Surface brightness ($\mu$Jy/beam)")
ax1.set_xlim(0, 5)
ax1.set_ylim(-50, 350)
ax1.plot(xaxis, SE_midplane_intensity, label='SE')
ax1.plot(xaxis, NW_midplane_intensity, '--', label='NW')
legend1 = ax1.legend()

ax2.set_title("Spine Intensity")
ax2.set_ylabel(r"Surface brightness ($\mu$Jy/beam)")
ax2.plot(xaxis, SE_amps, label='SE')
ax2.plot(xaxis, NW_amps, '--', label='NW')
ax2.set_ylim(-50, 350)
legend2 = ax2.legend()

ax3.set_title("Spine deviation from midplane")
ax3.set_ylabel("Elevation (\")")
ax3.plot(xaxis, SE_mus, label='SE')
ax3.plot(xaxis, NW_mus, '--', label='NW')
ax3.set_ylim(-0.5, 0.5)
legend3 = ax3.legend()

ax4.set_xlabel("Projected seperation from star (\")")
ax4.set_title("Disk Height")
ax4.set_ylabel(r'Fit FWHM (")')
ax4.plot(xaxis, SE_sigmas*2.3548200450, label='SE')
ax4.plot(xaxis, NW_sigmas*2.3548200450, '--', label='NW')
ax4.set_ylim(0, 1)
legend4 = ax4.legend()
plt.suptitle("Composite")
fig.savefig("boccaletti_plots_ctrpix_test_composite.png")
plt.show()
