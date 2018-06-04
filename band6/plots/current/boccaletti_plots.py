from scipy.ndimage.interpolation import rotate
from astropy.io import fits
from astropy.modeling import models, fitting
from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set_style("ticks")
sns.set_context("talk")

# Important parameters
angleSE = 128.2 - 90
sigma_to_FWHM = 2.35482004503

image = '../../cleans/current/band6_star_all.natural_clean'

# Read the header from the observed FITS continuum image:
head = fits.getheader(image + ".fits")
# Read in images and rotate so that disk is horizontal
im = rotate(fits.getdata(image + ".fits").squeeze(), angleSE, 
    reshape=False) * 1e6

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
distance = 9.725

# Define y extent of gaussian
SE_xpix_range = np.where((ra >  0) & (ra <= 5))[0][::-1]
NW_xpix_range = np.where((ra <= 0) & (ra >= -5))[0]
ra_range = ra[SE_xpix_range] * distance #au
ypix_range = np.where(abs(dec) < 10/distance)[0]



NW_amps, NW_mus, NW_sigmas = np.zeros(
    len(NW_xpix_range)), np.zeros(len(NW_xpix_range)), np.zeros(len(NW_xpix_range))
SE_amps, SE_mus, SE_sigmas = np.zeros(
    len(SE_xpix_range)), np.zeros(len(SE_xpix_range)), np.zeros(len(SE_xpix_range))

fitter = fitting.LevMarLSQFitter()
for i in range(len(NW_xpix_range)):
    
    NW_slice = im[ypix_range, NW_xpix_range[i]]
    
    #Fit gaussian to data
    NW_gauss_init = models.Gaussian1D(amplitude=np.max(NW_slice), mean=0, stddev=0.3)
    NW_gauss_fit = fitter(NW_gauss_init, dec[ypix_range], NW_slice)
    NW_amps[i], NW_mus[i], NW_sigmas[i] = NW_gauss_fit.parameters
    
for i in range(len(SE_xpix_range)):
    
    SE_slice = im[ypix_range, SE_xpix_range[i]]
    
    #Fit gaussian to data
    SE_gauss_init = models.Gaussian1D(amplitude=np.max(SE_slice), mean=0, stddev=0.3)
    SE_gauss_fit = fitter(SE_gauss_init, dec[ypix_range], SE_slice)
    SE_amps[i], SE_mus[i], SE_sigmas[i] = SE_gauss_fit.parameters
    
print(ra_range[np.where(SE_amps < 14.9 * 3)])    
print(ra_range[np.where(NW_amps < 14.9 * 3)[0][:-1]])    
    

# Subtract beam FWHM
bmin = head['bmin'] * 3600. / 2
bmaj = head['bmaj'] * 3600. / 2
bpa = head['bpa']; bpa

theta = (-bpa + angleSE) * np.pi/180
b_FWHM_y = 2* bmin * bmaj / np.sqrt((bmin * np.cos(theta))**2 +
            (bmaj * np.sin(theta))**2)
b_FWHM_x = 2* bmin * bmaj / np.sqrt((bmin * np.cos(theta + np.pi/4))**2 +
            (bmaj * np.sin(theta + np.pi/4))**2)
b_sigma_x = b_FWHM_x/sigma_to_FWHM
            
SE_disk_FWHM = np.sqrt((sigma_to_FWHM * SE_sigmas)**2 - b_FWHM_y**2)
NW_disk_FWHM = np.sqrt((sigma_to_FWHM * NW_sigmas)**2 - b_FWHM_y**2)

fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 10), sharex=True)
color = 'crimson'
# ax1.set_ylabel(r"Surface brightness ($\mu$Jy/beam)")
# ax1.set_ylim(-50, 350)
# ax1.plot(xaxis, SE_midplane_intensity, label='SE')
# ax1.plot(xaxis, NW_midplane_intensity, '--', label='NW')
# legend1 = ax1.legend()
ax1.set_ylabel(r"Surface brightness ($\mu$Jy/beam)")

SE = ax1.plot(ra_range, SE_amps, label='SE')
NW = ax1.plot(ra_range, NW_amps[:-1], '--', label='NW', color=color)
rms = 14.9
ax1.fill_between(ra_range, SE_amps-rms, SE_amps+rms, alpha=0.4)
ax1.fill_between(ra_range, NW_amps[:-1]-rms, NW_amps[:-1]+rms, alpha=0.4, color=color)

inset = zoomed_inset_axes(ax1, zoom = 1, loc=3)
xs = np.arange(-.7*10,.7*10, 0.01 * 10)
gauss = 100*np.exp(-xs**2/(2*(b_sigma_x*10)**2))
gauss_inset = inset.plot(xs, gauss, 'k:', lw=0.7)
inset.set_axis_off()

# ax1.axvline(1.02, ls=':', color='m')
# ax1.axvline(1.70, ls=':', color='m')
# ax1.axvline(2.96, ls=':', color='m')
# ax1.axvline(4.10, ls=':', color='m')
# ax1.axhline(3*14.9)
ax1.set_xlim(0, 4.5*10)
ax1.set_xticks(np.array([0,1,2,3,4])*10)
legend1 = ax1.legend(loc='upper right')

ax2.set_ylabel("Elevation (au)")
ax2.plot(ra_range, SE_mus*distance, label='SE')
ax2.plot(ra_range, NW_mus[:-1]*distance, '--', label='NW', color=color)
# ax2.axvline(1.02, ls=':', color='m')
# ax2.axvline(1.70, ls=':', color='m')
# ax2.axvline(2.96, ls=':', color='m')
# ax2.axvline(4.10, ls=':', color='m')
ax2.set_ylim(-0.15*10, 0.15*10)
# legend3 = ax2.legend(loc='upper left')

ax3.set_xlabel("Projected separation from star (au)")
ax3.set_ylabel(r'Beam-subtracted FWHM (au)')
ax3.plot(ra_range, SE_disk_FWHM * distance, label='SE')
ax3.plot(ra_range, NW_disk_FWHM[:-1] * distance, '--', label='NW', color=color)
# ax3.plot(ra_range, SE_sigmas * sigma_to_FWHM * distance, ':', label='SE unsubtracted')
# ax3.plot(ra_range, NW_sigmas[:-1] * sigma_to_FWHM * distance, '-.', label='NW unsubtracted', color=color)

# legend3 = ax3.legend(loc='upper right')
ax3.set_ylim(0, 0.7*10)
ax3.set_yticklabels([0, 2, 4, 6])
# legend4 = ax3.legend(loc='upper left')
# plt.suptitle("Composite")

plt.draw()
ylabel_x = ax2.yaxis.label.get_position()[0]
for ax in [ax1,ax2,ax3]:
    ylabel_y = ax.yaxis.label.get_position()[1]
    ax.yaxis.label.set_position((ylabel_x, ylabel_y))
    ax.yaxis._autolabelpos=False

fig.savefig("../../../writing/figures/boccaletti_plots.pdf", dpi=1000)
plt.show()

#Check if slice profile is Gaussian
# xpix = np.where((ra >= 2) & (ra <= 2.03))[0]
# # xpix = np.where((ra <= -3) & (ra >= -3.03))[0]
# ypix_range = np.where(abs(dec) < 1)[0]
# NW_gauss_init = models.Gaussian1D(amplitude=np.max(im[ypix_range, xpix]), mean=0, stddev=0.3)
# gauss_fit = fitter(gauss_init, dec[ypix_range], im[ypix_range, xpix])
# plt.plot(dec[ypix_range], im[ypix_range, xpix])
# plt.plot(dec[ypix_range], gauss_fit(dec[ypix_range]))
# plt.show()

# find offset of two (inner) local intensity maxima
# slice_indices = (5 < ra_range) & (ra_range < 12)
# slice_range =  ra_range[slice_indices]
# slice_range[SE_amps[slice_indices].argmax()]
# slice_range[NW_amps[:-1][slice_indices].argmax()]

b_FWHM_x
b_FWHM_y
