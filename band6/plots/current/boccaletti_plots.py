from scipy.ndimage.interpolation import rotate
from astropy.io import fits
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set_style("ticks")
sns.set_context("talk")

# Important parameters
angleSE = 128.48 - 90
sigma_to_FWHM = 2.35482004503

image = '../../cleans/current/band6_star_all.natural_clean'

# Read the header from the observed FITS continuum image:
head = fits.getheader(image + ".fits")
# Read in images and rotate so that disk is horizontal
im = rotate(fits.getdata(image + ".fits").squeeze(), angleSE, 
    reshape=False) * 1e6
plt.imshow(im, origin='lower')

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

# Get beam FWHM
bmin = head['bmin'] * 3600. / 2
bmaj = head['bmaj'] * 3600. / 2
bpa = head['bpa']; bpa
theta = -(-bpa + angleSE) * np.pi/180

(-bpa + angleSE)
(bpa - angleSE)

b_FWHM_y = 2* bmin * bmaj / np.sqrt((bmin * np.cos(theta))**2 +
            (bmaj * np.sin(theta))**2)
b_FWHM_x = 2* bmin * bmaj / np.sqrt((bmin * np.cos(theta + np.pi/4))**2 +
            (bmaj * np.sin(theta + np.pi/4))**2)
b_sigma_x = b_FWHM_x/sigma_to_FWHM
b_sigma_y = b_FWHM_y/sigma_to_FWHM

b_FWHM_x
b_FWHM_y

# Define y extent of gaussian
SE_xpix_range = np.where((ra >  0) & (ra <= 5))[0][::-1]
NW_xpix_range = np.where((ra <= 0) & (ra >= -5))[0]
au_range = ra[SE_xpix_range] * distance #au
ypix_range = np.where(abs(dec) < 5/distance)[0]

def gaussian(xs, amplitude, mean, sigma):
    return amplitude * np.exp(- 0.5 * (xs - mean) ** 2 / sigma ** 2)

NW_fits_list, SE_fits_list = [ [], [] ], [ [], [] ]
for i in range(len(NW_xpix_range)):
    #Fit gaussian to data
    NW_slice = im[ypix_range, NW_xpix_range[i]]
    NW_fit = curve_fit(f=gaussian, 
        xdata=dec[ypix_range]*distance, ydata=NW_slice, 
        sigma=[20]*len(NW_slice), absolute_sigma=True, 
        p0=(np.max(NW_slice), 0, 1.2))
    NW_fits_list[0].append(NW_fit[0]) # append best-fit A, mu, sigma
    NW_fits_list[1].append([NW_fit[1][0,0], NW_fit[1][1,1], NW_fit[1][2,2]]) # append variances
        
    if i != len(SE_xpix_range):
        SE_slice = im[ypix_range, SE_xpix_range[i]]
        SE_fit = curve_fit(f=gaussian, xdata=dec[ypix_range]*distance, ydata=SE_slice, 
            sigma=[15]*len(SE_slice), absolute_sigma=True, 
            p0=(np.max(SE_slice), 0, 1.2))
        SE_fits_list[0].append(SE_fit[0]) # append best-fit A, mu, sigma
        SE_fits_list[1].append([SE_fit[1][0,0], SE_fit[1][1,1], SE_fit[1][2,2]]) # append variances
    
NW_fits = np.array([NW_fits_list[0], np.sqrt(NW_fits_list[1])])
SE_fits = np.array([SE_fits_list[0], np.sqrt(SE_fits_list[1])])
            
import uncertainties.unumpy as unp
# Remove beam component
SE_square_difference = (sigma_to_FWHM * unp.uarray(*SE_fits[:,:,2]) )**2 - (b_FWHM_y*distance)**2
SE_square_difference[SE_square_difference < 0] = np.nan
SE_disk_FWHM = unp.sqrt(SE_square_difference)

NW_square_difference = (sigma_to_FWHM * unp.uarray(*NW_fits[:,:,2]) )**2 - (b_FWHM_y*distance)**2
NW_square_difference[NW_square_difference < 0] = np.nan
NW_disk_FWHM = unp.sqrt(NW_square_difference)

fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 10), sharex=False)
color = 'crimson'
# ax1.set_ylabel(r"Surface brightness ($\mu$Jy/beam)")
# ax1.set_ylim(-50, 350)
# ax1.plot(xaxis, SE_midplane_intensity, label='SE')
# ax1.plot(xaxis, NW_midplane_intensity, '--', label='NW')
# legend1 = ax1.legend()
ax1.set_ylabel(r"Surface brightness ($\mu$Jy/beam)")


SE = ax1.plot(au_range, SE_fits[0,:,0], label='SE')
NW = ax1.plot(au_range, NW_fits[0,:-1,0], '--', label='NW', color=color)
# ax1.fill_between(au_range, SE_fits[0,:,0] - SE_fits[1,:,0], SE_fits[0,:,0] - SE_fits[1,:,0], alpha=0.4)
ax1.fill_between(au_range, NW_fits[0,:-1,0] - NW_fits[1,:-1,0], NW_fits[0,:-1,0] + NW_fits[1,:-1,0], alpha=0.4, color=color)

# ax1.axvline(1.02, ls=':', color='m')
# ax1.axvline(1.70, ls=':', color='m')
# ax1.axvline(2.96, ls=':', color='m')
# ax1.axvline(4.10, ls=':', color='m')
# ax1.axhline(3*14.9)
# ax1.set_xticks(np.array([0,1,2,3,4])*10)
legend1 = ax1.legend(loc='upper right')

ax2.set_ylabel("Elevation (au)")
SE = ax2.plot(au_range, SE_fits[0,:,1], label='SE')
NW = ax2.plot(au_range, NW_fits[0,:-1,1], '--', label='NW', color=color)
ax2.fill_between(au_range, SE_fits[0,:,1] - SE_fits[1,:,1], SE_fits[0,:,1] + SE_fits[1,:,1], alpha=0.4)
ax2.fill_between(au_range, NW_fits[0,:-1,1] - NW_fits[1,:-1,1], NW_fits[0,:-1,1] + NW_fits[1,:-1,1], alpha=0.4, color=color)
# ax2.axvline(1.02, ls=':', color='m')
# ax2.axvline(1.70, ls=':', color='m')
# ax2.axvline(2.96, ls=':', color='m')
# ax2.axvline(4.10, ls=':', color='m')

ax3.set_xlabel("Projected separation from star (au)")
ax3.set_ylabel(r'Beam-subtracted FWHM (au)')

SE_disk_FWHM_val = unp.nominal_values(SE_disk_FWHM)
SE_disk_FWHM_sigma = unp.std_devs(SE_disk_FWHM)
NW_disk_FWHM_val = unp.nominal_values(NW_disk_FWHM)
NW_disk_FWHM_sigma = unp.std_devs(NW_disk_FWHM)
SE_disk_FWHM_val

ax3.plot(au_range, SE_disk_FWHM_val, label='SE')
ax3.plot(au_range, NW_disk_FWHM_val[:-1], '--', label='NW', color=color)
ax3.fill_between(au_range, SE_disk_FWHM_val - SE_disk_FWHM_sigma, SE_disk_FWHM_val + SE_disk_FWHM_sigma, alpha=0.4)
ax3.fill_between(au_range, NW_disk_FWHM_val[:-1] - NW_disk_FWHM_sigma[:-1], NW_disk_FWHM_val[:-1] + NW_disk_FWHM_sigma[:-1], color=color, alpha=0.4)
# ax3.plot(au_range, SE_sigmas * sigma_to_FWHM * distance, ':', label='SE unsubtracted')
# ax3.plot(au_range, NW_sigmas[:-1] * sigma_to_FWHM * distance, '-.', label='NW unsubtracted', color=color)

ax1.set_xlim(0, 43); ax2.set_xlim(0, 43); ax3.set_xlim(0, 43)
ax1.set_ylim(0, 500); ax2.set_ylim(-2.5, 2.5); ax3.set_ylim(0, 5)
ax1.set_xticklabels([]); ax2.set_xticklabels([])

# radial inset
inset_height = 1/4
inset = zoomed_inset_axes(ax1, zoom = 1, loc=3)
xs = np.arange(-7, 7, 0.1)
gauss = ax1.get_ylim()[1]*inset_height*np.exp(-xs**2/(2*(b_sigma_y*9.725)**2))
gauss_inset = inset.plot(xs, gauss, 'k:', lw=0.7)
inset.set_axis_off()

# vertical inset
inset2 = zoomed_inset_axes(ax3, zoom = 1, loc=3)
dec_inds = np.abs(dec) * distance < 7
xs = dec[dec_inds] * distance
A = ax3.get_ylim()[1] * inset_height
gauss = A*np.exp(-xs**2/(2*(b_sigma_x*distance)**2))

inds = np.where(np.abs(ra)*distance < 20)
for i in inds[0]:
    # i = np.random.randint(0, len(inds[0]))
    profile = im[dec_inds, i]
    profile *= A/profile.max()
    inset2.plot(xs, profile, 'k:', lw=0.7, alpha=0.7)
inset2.plot(xs, gauss, 'r-', lw=1)
inset2.set_axis_off()


plt.tight_layout()
plt.draw()
ylabel_x = ax2.yaxis.label.get_position()[0]
for ax in [ax1,ax2,ax3]:
    ylabel_y = ax.yaxis.label.get_position()[1]
    ax.yaxis.label.set_position((ylabel_x-13, ylabel_y))
    ax.yaxis._autolabelpos=False

fig.savefig("../../../writing/figures/boccaletti_plots.pdf", dpi=1000)
plt.show()


# find offset of two (inner) local intensity maxima
# slice_indices = (5 < au_range) & (au_range < 12)
# slice_range =  au_range[slice_indices]
# slice_range[SE_amps[slice_indices].argmax()]
# slice_range[NW_amps[:-1][slice_indices].argmax()]

# find 3-sigma extent
print(au_range[np.where(SE_fits[0,:,0] < 15 * 3)])    
print(au_range[np.where(NW_fits[0,:,0] < 15 * 3)[0][:-1]])    

# find average FWHM interior to 30
inds = np.where(au_range < 35)
FWHMs = np.array([SE_disk_FWHM[inds], NW_disk_FWHM[:-1][inds]])
FWHM_mean = np.mean(FWHMs)
FWHM_std = unp.sqrt( np.sum( (FWHMs - np.mean(FWHMs))**2 ) / (FWHMs.size-1) )
print(FWHM_mean, FWHM_std)
