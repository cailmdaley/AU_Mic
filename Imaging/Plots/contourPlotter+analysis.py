from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Ellipse
from matplotlib.ticker import MultipleLocator, LinearLocator, AutoMinorLocator
from scipy.ndimage.interpolation import rotate
from scipy.ndimage import zoom
from astropy.io import fits
import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np


def return_axis(ax, image, cpal, cbmin, cbmax, cbtmj, cbtmn, rms, cont_levs, text=None, residuals=None, axislabels=True):

    angleSE = 129.5 - 90
    angleNW = 311.2 - 90
    # Read the header from the observed FITS continuum image:
    head = fits.getheader(image + ".fits")
    # Read in images and rotate so that disk is horizontal
    im = rotate(fits.getdata(image + ".fits").squeeze(),
                angleSE, reshape=False)

    if residuals:
        resid = fits.getdata(residuals + ".fits").squeeze()
    np.shape(im)

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

    # Set axes limits
    xmin = -5.5
    xmax = 5.5
    ymin = -5.5
    ymax = 5.5
    ax.set_xlim(xmax, xmin)
    ax.set_ylim(ymin, ymax)
    ax.grid(False)

    # Set x and y major and minor tics
    majorLocator = MultipleLocator(1)
    ax.xaxis.set_major_locator(majorLocator)
    ax.yaxis.set_major_locator(majorLocator)

    minorLocator = MultipleLocator(0.2)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    # Set x and y labels
    if axislabels == True:
        ax.set_xlabel(r'$\Delta \alpha$ (")', fontsize=15)
        ax.set_ylabel(r'$\Delta \delta$ (")', fontsize=15)
        ax.xaxis.set_ticklabels(
            ['', '', '-4', '', '-2', '', '0', '', '2', '', '4', ''], fontsize=13)
        ax.yaxis.set_ticklabels(
            ['', '', '-4', '', '-2', '', '0', '', '2', '', '4', ''], fontsize=13)
    else:
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])

    # Set physical range of colour map
    cxmin = ra[0]
    cxmax = ra[-1]
    cymin = dec[0]
    cymax = dec[-1]

    # Set limits and tics of colorbar - velocity scale
    print(np.min(im))
    print(np.max(im))

    # Plot image as a colour map in units of micro Jy
    im = im * 1e6
    cmap = ax.imshow(im,
                     extent=[cxmin, cxmax, cymin, cymax],
                     vmin=np.min(im),
                     vmax=np.max(im),
                     origin='lower',
                     cmap=cpal)

    # Scale countour levels to micro Jy and plot contours
    cont_levs = cont_levs * rms * 1e6

    if residuals:
        ax.contour(resid,
                   levels=cont_levs,
                   colors='k',
                   linewidths=0.75,
                   linestyles='solid')
        ax.contour(resid,
                   levels=-1 * np.flip(cont_levs, axis=0),
                   colors='k',
                   linewidths=0.75,
                   linestyles='dashed')
    else:
        ax.contour(ra, dec, im,
                   colors='k',
                   levels=cont_levs,
                   linewidths=0.75,
                   linestyles='solid')
        ax.contour(ra, dec, im,
                   levels=-1 * np.flip(cont_levs, axis=0),
                   colors='k',
                   linewidths=0.75,
                   linestyles='dashed')

    # Create the colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="8%", pad=0.0)
    cbar = fig.colorbar(cmap, ax=ax, cax=cax, orientation='horizontal')
    cbar.ax.xaxis.set_tick_params(direction='out',
                                  length=3,
                                  which='major',
                                  bottom='off',
                                  top='on',
                                  labelsize=8,
                                  pad=-2,
                                  labeltop='on',
                                  labelbottom='off')
    cbar.ax.xaxis.set_tick_params(direction='out',
                                  length=2,
                                  which='minor',
                                  bottom='off',
                                  top='on')
    minorLocator = AutoMinorLocator(cbtmj / cbtmn)
    cbar.ax.xaxis.set_minor_locator(minorLocator)
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),
                            rotation=45, fontsize=11)
    cbar.set_ticks(np.arange(cbmin, cbmax, cbtmj))

    # Colorbar label
    cbar.ax.text(0.432, 0.365, r'$\mu Jy / bm$', fontsize=13,
                 path_effects=[PathEffects.withStroke(linewidth=2, foreground="w")])

    # Overplot the beam ellipse
    beam_ellipse_color = 'k'
    bmin = head['bmin'] * 3600.
    bmaj = head['bmaj'] * 3600.
    bpa = head['bpa']
    el = Ellipse(xy=[-4.6, -4.6],
                 width=bmin,
                 height=bmaj,
                 angle=-bpa + angleSE,
                 edgecolor='k',
                 hatch='/////',
                 facecolor='none',
                 zorder=10)
    ax.add_artist(el)

    # Plot the scale bar
    x = 3.97
    y = -5
    ax.plot([x, x - 1], [y, y], '-', linewidth=2, color='k')
    ax.text(x + 0.25, y + 0.15, "10 au", fontsize=17,
            path_effects=[PathEffects.withStroke(linewidth=2, foreground="w")])

    # Plot a cross at the source position
    ax.plot([0.0], [0.0], '+', markersize=307, markeredgewidth=2, color='k')

    # Add PA lines and legend
    ax.plot([5, -5], [0, 0], '-', linewidth=1,
            color='k', label='Boccaletti SE PA')
    h = 5 * np.sin(np.radians(angleNW - angleSE))
    ax.plot([5, -5], [-h, h], '--', linewidth=1,
            color='blue', label='Boccaletti NW PA')
    legend = ax.legend()

    # Add cardinal directions:
    y = -5
    ax.arrow(0, y, np.cos(np.radians(angleSE)), np.sin(np.radians(angleSE)),
             width=0.04, head_width=0.2, head_length=0.3, fc='k', ec='k')
    ax.text(0 + np.cos(np.radians(angleSE)), y + np.sin(np.radians(angleSE)) + 0.15,
            'E', path_effects=[PathEffects.withStroke(linewidth=2, foreground="w")])
    ax.arrow(0, y, np.cos(np.radians(90 + angleSE)), np.sin(np.radians(90 +
                                                                       angleSE)), width=0.04, head_width=0.2, head_length=0.3, fc='k', ec='k')
    ax.text(0 + np.cos(np.radians(90 + angleSE)) + 0.4, y + np.sin(np.radians(90 + angleSE)
                                                                   ) + 0.2, 'N', path_effects=[PathEffects.withStroke(linewidth=2, foreground="w")])

    # Add figure text
    # if text:
    #     for t in text:
    #         ax.text(*t, fontsize=18,
    # path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])


# Set seaborn plot styles
sns.set_style("ticks",
              {"xtick.direction": "in",
               "ytick.direction": "in"})
sns.set_context("talk")


# Set colorpalette
from my_colormaps import *
cpals = [cubehelix_1, cubehelix_2, cubehelix_3, jesse_reds]
colorsaves = ['cubehelix_1', 'cubehelix_2', 'cubehelix_3', 'jesse_reds']
which_cpal = 3

# Create figure
fig, (natural_ax, taper_ax) = plt.subplots(
    1, 2, sharex=False, sharey=False, figsize=(11.6, 6.2))

# Plot subplots on seperate axes
return_axis(ax=natural_ax,
            image = '../cleans/aumic_ctrpix_test_usermask_natural',
            axislabels=True,
            cpal=cpals[which_cpal],
            cbmin=-50,
            cbmax=351,
            cbtmj=50,
            cbtmn=10,
            rms=1.77794445335e-05,
            cont_levs=np.arange(2, 40, 2),
            text=[(4.8, 4.4, 'AU Mic ALMA 1.4mm'),
                  (4.43, 3.95, 'natural weighting')])


return_axis(ax=taper_ax,
            image='../cleans/aumic_junmar_usermask_200klam',
            axislabels=False,
            cpal=cpals[which_cpal],
            cbmin=-50,
            cbmax=601,
            cbtmj=100,
            cbtmn=20,
            rms=2.23016886594e-05,
            cont_levs=np.arange(2, 40, 2),
            text=[(4.8, 4.4, 'AU Mic ALMA 1.4mm'),
                  (4, 3.95, r'200k$\lambda$ taper')])

plt.subplots_adjust(wspace=0)

# Save and show figure
plt.savefig('AU_mic_composite_ctrpix_test_PA_lines.png')
plt.show()
