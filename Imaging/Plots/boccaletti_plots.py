from scipy.ndimage.interpolation import rotate
from scipy.optimize import curve_fit
from astropy.io import fits
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


image = 'aumic_usermask_natural'
intensity=True
spine=False

# Read the header from the observed FITS continuum image:
head = fits.getheader(image + ".fits")
# Read in images and rotate so that disk is horizontal
angleSE = 129.5 - 90
im = rotate(fits.getdata(image + ".fits").squeeze(),
            angleSE, reshape=False) * 1e6
im[:, :40] = 0

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

raSE = ra[:257]
print(ra[:15])

if intensity == True:
    midplane_intensity = np.zeros(len(raSE))
    for i in range(len(raSE)):
        midplane_intensity[i] = np.mean(im[255:258, i])

    sns.set_style('whitegrid')
    fig, ax = plt.subplots()
    ax.set_title("Midplane intensity")
    ax.set_xlabel("Projected seperation from star (\")")
    ax.set_ylabel(r"Surface brightness ($\mu$Jy)...")
    ax.set_xlim(0, 7)
    ax.set_ylim(-50, 350)
    ax.plot(ra[:257], midplane_intensity)
    # ax.invert_xaxis()
    fig.savefig('midplane_intensity.png')
    plt.show()

if spine == True:
    # Define model Gaussian function:
    def gauss(x, *p):
        A, mu, sigma = p
        return A*np.exp(-(x-mu)**2/(2.*sigma**2))

    p0 = [100., 0., 0.3]
    bounds = ([-50, -0.5, 0], [400, 0.5, 1])
    amps, mus, sigmas = np.zeros(len(raSE)), np.zeros(len(raSE)), np.zeros(len(raSE))

    for i in range(256):
        (amps[i], mus[i], sigmas[i]), var_matrix = curve_fit(gauss, dec[256-35:256+35], im[256-35:256+35,len(raSE)-i], p0=p0)

    sns.set_style('whitegrid')
    fig, ax = plt.subplots()
    ax.set_title("Gaussian fit statistics")
    ax.set_xlabel("Projected seperation from star (\")")
    ax.set_ylabel(r"")
    ax.set_xlim(0, 6)
    ax.plot(raSE[::-1], sigmas)
    ax.invert_xaxis()
    fig.savefig("gauss_fit.png")
    plt.show()
