from astropy.io import fits
from astropy.modeling import models, fitting
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Ellipse
from matplotlib.ticker import MultipleLocator, LinearLocator, AutoMinorLocator
from my_colormaps import *
import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set seaborn plot styles and color pallete
sns.set_style("ticks",
              {"xtick.direction": "in",
               "ytick.direction": "in"})
sns.set_context("talk")
cpal = jesse_reds


class Observation:
    def __init__(self, filename, rms, fig=None, pos=(0,1), **kwds):

        self.file = filename
        self.rms = rms
        self.fig = fig
        self.pos = pos
        self.__dict__.update(kwds)


    def get_fits(self):
        self.head = fits.getheader(self.file)
        self.im = fits.getdata(self.file).squeeze()

        # change units to micro Jy
        self.im *= 1e6
        self.rms *= 1e6

        # Read in header spatial info to create ra
        nx = self.head['NAXIS1']
        xpix = self.head['CRPIX1']
        xval = self.head['CRVAL1']
        self.xdelt = self.head['CDELT1']

        ny = self.head['NAXIS2']
        ypix = self.head['CRPIX2']
        yval = self.head['CRVAL2']
        self.ydelt = self.head['CDELT2']

        # Convert from degrees to arcsecs
        self.ra_offset = np.array(((np.arange(nx) - xpix + 1) * self.xdelt) * 3600)
        self.dec_offset = np.array(((np.arange(ny) - ypix + 1) * self.ydelt) * 3600)

        # self.ra_abs = ((np.arange(nx) - xpix + 1) * xdelt + xval) * 3600
        # self.dec_abs = ((np.arange(ny) - ypix + 1) * ydelt + yval) * 3600
        return


    def make_axis(self):
        if not self.fig:
            self.fig, self.ax = plt.subplots()
        else:
            self.ax = self.fig.axes[self.pos[0]]

        xmin = -5.0
        xmax = 5.0
        ymin = -5.0
        ymax = 5.0
        self.ax.set_xlim(xmax, xmin)
        self.ax.set_ylim(ymin, ymax)
        self.ax.grid(False)

        # Set x and y major and minor tics
        majorLocator = MultipleLocator(1)
        self.ax.xaxis.set_major_locator(majorLocator)
        self.ax.yaxis.set_major_locator(majorLocator)

        minorLocator = MultipleLocator(0.2)
        self.ax.xaxis.set_minor_locator(minorLocator)
        self.ax.yaxis.set_minor_locator(minorLocator)

        # Set x and y labels
        self.ax.set_xlabel(r'$\Delta \alpha$ (")', fontsize=15)
        self.ax.set_ylabel(r'$\Delta \delta$ (")', fontsize=15)
        self.ax.xaxis.set_ticklabels(
            ['', '', '-4', '', '-2', '', '0', '', '2', '', '4', ''], fontsize=13)
        self.ax.yaxis.set_ticklabels(
            ['', '', '-4', '', '-2', '', '0', '', '2', '', '4', ''], fontsize=13)
        self.ax.tick_params(which='both', right='on')

        # Set labels depending on position in figure
        if self.pos[0] == 0: #left
            self.ax.tick_params(axis='y', labelright='off')
        elif self.pos[0] == self.pos[1]-1: #right
            self.ax.set_xlabel('')
            self.ax.set_ylabel('')
            self.ax.tick_params(axis='y', labelleft='off', labelright='on')
        else: #middle
            self.ax.tick_params(axis='y', labelleft='off')
            self.ax.set_xlabel('')
            self.ax.set_ylabel('')

        # Set physical range of colour map
        self.extent = [self.ra_offset[0], self.ra_offset[-1],
                      self.dec_offset[-1], self.dec_offset[0]]

        return


    def fill_axis(self):
        # Plot image as a colour map
        cmap = self.ax.imshow(self.im,
                             extent=self.extent,
                             vmin=np.min(self.im),
                             vmax=np.max(self.im),
                             cmap=cpal)

        # Set contour levels
        cont_levs = np.arange(2, 40, 2) * self.rms

        # add residual contours if resdiual exists; otherwise, add image contours
        try:
            self.ax.contour(self.resid,
                           levels=cont_levs,
                           colors='k',
                           linewidths=0.75,
                           linestyles='solid')
            self.ax.contour(self.resid,
                           levels=-1 * np.flip(cont_levs, axis=0),
                           colors='k',
                           linewidths=0.75,
                           linestyles='dashed')
        except AttributeError:
            self.ax.contour(self.ra_offset, self.dec_offset, self.im,
                           colors='k',
                           levels=cont_levs,
                           linewidths=0.75,
                           linestyles='solid')
            self.ax.contour(self.ra_offset, self.dec_offset, self.im,
                           levels=-1 * np.flip(cont_levs, axis=0),
                           colors='k',
                           linewidths=0.75,
                           linestyles='dashed')

        # Create the colorbar
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("top", size="8%", pad=0.0)
        cbar = self.fig.colorbar(cmap, ax=self.ax, cax=cax, orientation='horizontal')
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
        # If major/minor spacings exist, use them; otherwise, use default value
        try:
            self.cbtmj
            self.cbtmn
        except AttributeError:
            self.cbtmj = 100
            self.cbtmn = 20

        minorLocator = AutoMinorLocator(self.cbtmj / self.cbtmn)
        cbar.ax.xaxis.set_minor_locator(minorLocator)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),
                                rotation=45, fontsize=11)
        cbar.set_ticks(np.arange(-10000, 10000, self.cbtmj))

        # Colorbar label
        cbar.ax.text(0.432, 0.365, r'$\mu Jy / bm$', fontsize=13,
                     path_effects=[PathEffects.withStroke(linewidth=2, foreground="w")])

        # Overplot the beam ellipse
        beam_ellipse_color = 'k'
        bmin = self.head['bmin'] * 3600.
        bmaj = self.head['bmaj'] * 3600.
        bpa = self.head['bpa']

        el = Ellipse(xy=[-4.2, -4.2],
                     width=bmin,
                     height=bmaj,
                     angle=-bpa,
                     edgecolor='k',
                     hatch='/////',
                     facecolor='none',
                     zorder=10)

        self.ax.add_artist(el)

        # Plot the scale bar
        x = 3.97
        y = -4.7
        self.ax.plot(
            [x, x - 1],
            [y, y],
            '-', linewidth=2, color='k')
        self.ax.text(
            x + 0.2, y + 0.15, "10 au",
            fontsize=17,
            path_effects=[PathEffects.withStroke(linewidth=2, foreground="w")])

        # Plot a cross at the source position
        self.ax.plot([0.0], [0.0], '*', markersize=9, markeredgewidth=1, color='k')

        # Add figure text
        try:
            for t in self.text:
                self.ax.text(*t, fontsize=18,
                    path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])

        except AttributeError:
            pass


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

num = 3
plot_width=11.6/2
fig = plt.subplots(1, num, sharex=False, sharey=False,
    figsize=(plot_width * num, 6.2))[0]
plt.subplots_adjust(wspace=-0.1)

mar = Observation('../cleans/aumic_mar_centered_200klam.fits', 2.92637287203e-05,
    fig=fig, pos=(0, num),
    text=[[4.7, 4.4, 'AU Mic ALMA 1.4mm'],
          [4.7, 3.4, 'March 2014']])
aug = Observation('../cleans/aumic_aug_centered_natural.fits', 3.83357546525e-05,
    fig=fig, pos=(1, num),
    text=[[4.7, 4.4, 'AU Mic ALMA 1.4mm'],
          [4.7, 3.4, 'August 2014']])

# jun = Observation('../cleans/aumic_jun_timing_centered_natural.fits', 2.80185104202e-05,
#     fig=fig, pos=(2, num),
#     text=[[4.7, 4.4, 'AU Mic ALMA 1.4mm'],
#           [4.7, 3.4, 'June 2015']])


all_natural = Observation('../../cleans/current/aumic_all_natural.fits',
    1.4494822607957758e-05, fig=fig, pos=(0, num),
    text=[[4.7, 4.4, 'AU Mic ALMA 1.4mm'],
          [4.7, 3.4, 'Natural']])
all_briggs = Observation('../../cleans/current/aumic_all_briggs.fits', 
    1.578679439262487e-05, fig=fig, pos=(1, num), 
    text=[[4.7, 4.4, 'AU Mic ALMA 1.4mm'],
          [4.7, 3.4, 'Robust=0.5']])
          
all_uniform = Observation('../../cleans/current/aumic_all_uniform.fits', 
    3.4449552913429216e-05, fig=fig, pos=(2, num), 
    text=[[4.7, 4.4, 'AU Mic ALMA 1.4mm'],
          [4.7, 3.4, 'Uniform']])


for obs in [all_natural, all_briggs, all_uniform]:
    obs.get_fits()
    obs.make_axis()
    obs.fill_axis()

plt.savefig('aumic_triptych_all_different_weightings.png')
plt.show()

# plt.show()