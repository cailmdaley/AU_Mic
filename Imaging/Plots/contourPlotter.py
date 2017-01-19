from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Ellipse
from matplotlib.ticker import MultipleLocator, LinearLocator, AutoMinorLocator
import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
from astropy.io import fits
import seaborn as sns
import numpy as np


def plotcmd(image, cpal, cbmin, cbmax, cbtmj, cbtmn, rms, cont_levs, text=None, residuals=None, save=None, show=False):

    # Set seaborn plot styles
    sns.set_style("ticks",
                  {"xtick.direction": "in",
                   "ytick.direction": "in"})
    sns.set_context("talk")

    # Read the header from the observed FITS continuum image:
    head = fits.getheader(image + ".fits")
    im = fits.getdata(image + ".fits").squeeze()
    if residuals:
        resid = fits.getdata(residuals + ".fits").squeeze()

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

    # Create figure
    fig, ax = plt.subplots(figsize=(6, 6))

    # Set axes limits
    xmin = -4.0
    xmax = 4.0
    ymin = -4.0
    ymax = 4.0
    ax.set_xlim(xmax, xmin)
    ax.set_ylim(ymin, ymax)
    ax.grid(False)

    # Set x and y labels
    ax.set_xlabel(r'$\Delta \alpha$ (")')
    #   labelpad=15, fontsize=18)
    ax.set_ylabel(r'$\Delta \delta$ (")')
    #   labelpad=15, fontsize=18)

    # Set x and y major and minor tics
    majorLocator = MultipleLocator(1)
    ax.xaxis.set_major_locator(majorLocator)
    ax.yaxis.set_major_locator(majorLocator)

    minorLocator = MultipleLocator(0.2)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    # Set physical range of colour map
    cxmin = ra[0]
    cxmax = ra[-1]
    cymin = dec[-1]
    cymax = dec[0]

    # Set limits and tics of colorbar - velocity scale
    print(np.min(im))
    print(np.max(im))

    # Plot image as a colour map in units of micro Jy
    im = im * 1e6
    cmap = plt.imshow(im,
                      extent=[cxmin, cxmax, cymin, cymax],
                      vmin=cbmin,
                      vmax=np.max(im),
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
    cbar = fig.colorbar(cmap, cax=cax, orientation='horizontal')
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
    cbar.set_ticks(np.arange(cbmin, cbmax, cbtmj))

    # Colorbar label
    fig.text(0.477, 0.862, r'$\mu Jy / bm$', fontsize=10,
             path_effects=[PathEffects.withStroke(linewidth=2, foreground="w")])

    # Overplot the beam ellipse
    beam_ellipse_color = 'k'
    bmin = head['bmin'] * 3600.
    bmaj = head['bmaj'] * 3600.
    bpa = head['bpa']
    el = Ellipse(xy=[-3.5, -3.5],
                 width=bmin,
                 height=bmaj,
                 angle=-bpa,
                 edgecolor='k',
                 hatch='/////',
                 facecolor='none',
                 zorder=10)
    ax.add_artist(el)

    # Plot the scale bar
    x = 3
    y = -3.7
    ax.plot([x, x - 1], [y, y], '-', linewidth=1, color='k')
    ax.text(x - 0.14, y + 0.06, "10 au", fontsize=11,
            path_effects=[PathEffects.withStroke(linewidth=2, foreground="w")])

    # Plot a cross at the source position
    ax.plot([0.0], [0.0], '*', markersize=9, markeredgewidth=1, color='k')

    # Add figure text
    if text:
        for t in text:
            ax.text(*t, fontsize=11,
                    path_effects=[PathEffects.withStroke(linewidth=2, foreground="w")])

    # Save and show figure
    if save:
        fig.savefig(save)
    if show:
        plt.show()

# Set colorpalette
from my_colormaps import *
cpals = [cubehelix_1, cubehelix_2, cubehelix_3, jesse_reds]
colorsaves = ['cubehelix_1', 'cubehelix_2', 'cubehelix_3', 'jesse_reds']

i = 3
plotcmd('aumic_usermask_natural',
        cpal=cpals[i],
        cbmin=-50,
        cbmax=301,
        cbtmj=50,
        cbtmn=25,
        rms=1.4753316463611554e-05,
        cont_levs=np.arange(2, 40, 2),
        text=[(1.3, -3.4, 'AU Mic ALMA 1.4mm'),
              (1.04, -3.7, 'natural weighting')],
        save='aumic_usermask_natural' + '_countour_image_' +
        colorsaves[i] + '.png',
        show=True)

plotcmd('aumic_usermask_natural_200klam',
        cpal=cpals[i],
        cbmin=-50,
        cbmax=551,
        cbtmj=100,
        cbtmn=50,
        rms=1.9399181837798096e-05,
        cont_levs=np.arange(2, 40, 2),
        text=[(1.3, -3.4, 'AU Mic ALMA 1.4mm'),
              (0.72, -3.7, r'200k$\lambda$ taper')],
        save='aumic_usermask_natural_200klam' +
        '_countour_image_' + colorsaves[i] + '.png',
        show=False)
