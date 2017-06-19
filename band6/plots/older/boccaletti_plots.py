from scipy.ndimage.interpolation import rotate
from scipy.optimize import curve_fit
from astropy.io import fits
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# Important parameters
angleSE = 128.6 - 90

image = '../cleans/aumic_composite_natural'
# image = '../cleans/aumic_ctrpix_test_usermask_natural'

# Read the header from the observed FITS continuum image:
head = fits.getheader(image + ".fits")
# Read in images and rotate so that disk is horizontal
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
xaxis = ra[2:int(xpix)][::-1]

# Define y extent of gaussian
y_extent = 1 #arcsec
y_below = int(np.round((ypix - (y_extent / 3600) / ydelt)))
y_above = int(np.round((ypix + (y_extent / 3600) / ydelt)))
y_below, y_above

#Find intensity of midplane
# SE_midplane_intensity = np.zeros(len(xaxis))
# NW_midplane_intensity = np.zeros(len(xaxis))
# for i in range(len(xaxis)):
#     NW_midplane_intensity[i] = np.mean(
#         im[y_below:y_above, int(xpix)+i])
#     SE_midplane_intensity[i] = np.mean(
#         im[y_below:y_above, int(xpix)-i])


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
    print(i)
    (NW_amps[i], NW_mus[i], NW_sigmas[i]), NW_var_matrix = \
        curve_fit(gauss, dec[y_below:y_above],
        im[y_below:y_above, int(xpix)+i], p0=p0, bounds=bounds)

    (SE_amps[i], SE_mus[i], SE_sigmas[i]), SE_var_matrix = \
        curve_fit(gauss, dec[y_below:y_above],
        im[y_below:y_above, int(xpix)-i], p0=p0, bounds=bounds)

# Subtract beam FWHM
bmin = head['bmin'] * 3600. / 2
bmaj = head['bmaj'] * 3600. / 2
bpa = head['bpa']; bpa

theta = (-bpa + angleSE) * np.pi/180
b_FWHM_y = 2* bmin * bmaj / np.sqrt((bmin * np.cos(theta))**2 +
            (bmaj * np.sin(theta))**2)
b_FWHM_y
SE_disk_FWHM = np.sqrt((2.35482004503 * SE_sigmas)**2 - b_FWHM_y**2)
NW_disk_FWHM = np.sqrt((2.35482004503 * NW_sigmas)**2 - b_FWHM_y**2)

sns.set_style("whitegrid")
sns.set_context("talk")
fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 10), sharex=True)

# ax1.set_title("Midplane intensity")
# ax1.set_ylabel(r"Surface brightness ($\mu$Jy/beam)")
# ax1.set_ylim(-50, 350)
# ax1.plot(xaxis, SE_midplane_intensity, label='SE')
# ax1.plot(xaxis, NW_midplane_intensity, '--', label='NW')
# legend1 = ax1.legend()

ax1.set_title("Spine Intensity")
ax1.set_ylabel(r"Surface brightness ($\mu$Jy/beam)")
ax1.plot(xaxis, SE_amps, label='SE')
ax1.plot(xaxis, NW_amps, '--', label='NW')
# ax1.axvline(1.02, ls=':', color='m')
# ax1.axvline(1.70, ls=':', color='m')
# ax1.axvline(2.96, ls=':', color='m')
# ax1.axvline(4.10, ls=':', color='m')
ax1.set_xlim(0, 4.5)
ax1.set_ylim(-50, 350)
ax1.set_xticks([0,1,2,3,4])
legend2 = ax1.legend(loc=4)

ax2.set_title("Spine deviation from midplane")
ax2.set_ylabel("Elevation (\")")
ax2.plot(xaxis, SE_mus, label='SE')
ax2.plot(xaxis, NW_mus, '--', label='NW')
ax2.set_ylim(-0.25, 0.25)
legend3 = ax2.legend(loc=4)

ax3.set_xlabel("Projected seperation from star (\")")
ax3.set_title("Disk FWHM")
ax3.set_ylabel(r'Beam-subtracted FWHM (")')
ax3.plot(xaxis, SE_disk_FWHM, label='SE')
ax3.plot(xaxis, NW_disk_FWHM, '--', label='NW')
# ax3.plot(xaxis, SE_sigmas * 2.3548200450, label='SE')
# ax3.plot(xaxis, NW_sigmas * 2.3548200450, '--', label='NW')
ax3.set_ylim(0, 0.5)
legend4 = ax3.legend(loc=4)
plt.suptitle("Composite")
fig.savefig("boccaletti_plots_composite.png")
plt.show()
