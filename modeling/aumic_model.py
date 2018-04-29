import numpy as np
import argparse
import subprocess as sp; import os
from astropy.io import fits
from collections import OrderedDict
import matplotlib.pyplot as plt 
import matplotlib.colors as colors

from astrocail import fitting, plotting, mcmc
from disk_model import debris_disk_python3, raytrace
import aumic_fitting

#initialize an instance of class Disk with specified parameters
# theta = np.fromstring(string[3:-1], sep='   ')
param_dict = OrderedDict([
    ('temp_index',        -0.5),
    ('m_disk',            -7.544664), # log
    ('sb_law',            0.87462),
    ('r_in',              23.86895),
    ('d_r',               18.519918),
    ('r_crit',            150.0),
    ('inc',               88.337326),
    ('m_star',            0.31),
    ('co_frac',           0.0001),
    ('v_turb',            0.081),
    ('Zq',                70.0),
    ('column_densities', [0.79, 1000]),
    ('abundance_bounds', [50, 500]),
    ('hand',              -1),
    ('rgrid_size',        1000),
    ('zgrid_size',        500),
    ('l_star',            0.09),
    ('scale_factor',      0.1),
    ('pa',                128.483764),
    ('mar_starflux',      0.000366),
    ('aug_starflux',      0.000141),
    ('jun_starflux',      0.000211),
    ('annulus_r_in',      10.840851),
    ('annulus_r_out',     .05), # annulus_r_in gets added to this later
    ('annulus_mass',      0.000167)])
    

def make_fits(model, disk_params, annulus_params):
    structure_params = disk_params[:-1]
    PA = disk_params[-1]

    model_disk = debris_disk.Disk(
        structure_params,
        obs=[300, 351, 300, 5],
        annulus=annulus_params)
    raytrace.total_model(model_disk,
        distance=9.91, # pc
        imres=0.03, # arcsec/pix
        xnpix=512, #image size in pixels
        freq0=230.56459558740002,
        PA=PA,
        offs=[0.0,0.0], # offset from image center
        nchans=1, # continum
        isgas=False, # continuum!
        includeDust=True, #continuuum!!
        extra=0, # ?
        modfile = model.root + model.name)
    return model_disk

def fix_fits(model, obs, starflux):
    # open fits and add starflux,
    model_fits = fits.open(model.path + '.fits')
    crpix = int(model_fits[0].header['CRPIX1'])
    model_im = model_fits[0].data[0]
    try:
        model_im[crpix, crpix] = model.crpix_diskflux + starflux
    except AttributeError:
        model.crpix_diskflux = model_im[crpix, crpix]
        model_im[crpix, crpix] += starflux

    model_fits[0].header['CRVAL1'] = obs.ra
    model_fits[0].header['CRVAL2'] = obs.dec

    model_fits.writeto(model.path + '.fits', overwrite=True)
# def fix_fits(model, starflux):
#     # open fits and add starflux,
#     model_fits = fits.open(model.path + '.fits')
#     crpix = int(model_fits[0].header['CRPIX1'])
#     model_im = model_fits[0].data[0]
#     try:
#         model_im[crpix, crpix] = model.crpix_diskflux + starflux
#     except AttributeError:
#         model.crpix_diskflux = model_im[crpix, crpix]
#         model_im[crpix, crpix] += starflux
# 
#     model_fits[0].header['CRVAL1'] = obs.ra
#     model_fits[0].header['CRVAL2'] = obs.dec
#     model_fits.writeto(model.path + '.fits', overwrite=True)

def make_fit_plot(concise=False):
    param_dict['m_disk'] = 10**param_dict['m_disk']
    param_dict['d_r'] += param_dict['r_in']
    param_dict['annulus_r_out'] += param_dict['annulus_r_in']

    disk_params    = list(param_dict.values())[:-6]
    starfluxes     = list(param_dict.values())[-6:-3]
    annulus_params = list(param_dict.values())[-3:]

    # intialize model and make fits image
    print('Making model...')
    model = fitting.Model(observations=aumic_fitting.band6_observations,
        root='model_files/', name='bestfit')
    make_fits(model, disk_params, annulus_params)
    
    # disk_image = fits.open(model.path + '.fits')[0].data
    # plt.figure(figsize=(10,10))
    # plt.imshow(disk_image[0,:,:], origin="lower", cmap="inferno",
    #     norm=colors.LogNorm(vmin=disk_image.max()*1e-4, vmax=disk_image.max()))
    # plt.colorbar()
    # plt.show()

    print('Sampling and cleaning...')
    paths = []
    for pointing, rms, starflux in zip(model.observations, aumic_fitting.band6_rms_values[:-1], starfluxes):
        ids = []
        for obs in pointing:
            fix_fits(model, obs, starflux)

            ids.append('_' + obs.name[12:20])
            model.obs_sample(obs, ids[-1])
            model.make_residuals(obs, ids[-1])

        cat_string1 = ','.join([model.path+ident+'.vis' for ident in ids])
        cat_string2 = ','.join([model.path+ident+'.residuals.vis' for ident in ids])
        paths.append('{}_{}'.format(model.path, obs.name[12:15]))

        sp.call(['uvcat', 'vis={}'.format(cat_string2), 'out={}.residuals.vis'.format(paths[-1])], stdout=open(os.devnull, 'wb'))
        sp.call(['uvcat', 'vis={}'.format(cat_string1), 'out={}.vis'.format(paths[-1])], stdout=open(os.devnull, 'wb'))

        model.clean(paths[-1] + '.residuals', rms, show=False)
        model.clean(paths[-1], rms, show=False)


    cat_string1 = ','.join([path + '.vis' for path in paths])
    cat_string2 = ','.join([path + '.residuals.vis' for path in paths])

    sp.call(['uvcat', 'vis={}'.format(cat_string1), 'out={}_all.vis'.format(model.path)], stdout=open(os.devnull, 'wb'))
    sp.call(['uvcat', 'vis={}'.format(cat_string2), 'out={}_all.residuals.vis'.format(model.path)], stdout=open(os.devnull, 'wb'))

    model.clean(model.path+'_all', aumic_fitting.band6_rms_values[-1], show=False)
    model.clean(model.path+'_all.residuals', aumic_fitting.band6_rms_values[-1], show=False)

    paths.append('{}_all'.format(model.path))

    print('Making figure...')
    if concise:
        fig = plotting.Figure(
            layout=(1,3),
            paths=[
                aumic_fitting.band6_fits_images[-1],
                paths[-1] + '.fits',
                paths[-1] + '.residuals.fits'],
            rmses=3*[aumic_fitting.band6_rms_values[-1]],
            texts=[
                [[4.6, 4.0, 'Data']],
                [[4.6, 4.0, 'Model']],
                [[4.6, 4.0, 'Residuals']]
                ],
            title=None, #r'Run 6 Global Best Fit Model & Residuals',
            savefile='test' + '_bestfit_concise.pdf')
    else:
        fig = plotting.Figure(layout=(4,3),
           paths=[[obs, path + '.fits', path + '.residuals.fits']
               for obs, path in zip(aumic_fitting.band6_fits_images, paths)],
           rmses=[3*[rms] for rms in aumic_fitting.band6_rms_values],
           texts=[
               [[[4.6, 4.0, date]],
               [[4.6, 4.0, 'rms={}'.format(np.round(rms*1e6))]], None]
               for date, rms in zip(['March', 'August', 'June', 'All'],
               aumic_fitting.band6_rms_values)],
           title= r'Global Best Fit Model & Residuals',
          savefile='test' + '_bestfit_global.pdf')
          
# make_fit_plot()
          
# param_dict['m_disk'] = 10**param_dict['m_disk']
# param_dict['d_r'] += param_dict['r_in']
# param_dict['annulus_r_out'] += param_dict['annulus_r_in']
# 
# disk_params    = list(param_dict.values())[:-6]
# starfluxes     = list(param_dict.values())[-6:-3]
# annulus_params = list(param_dict.values())[-3:]
# 
# # intialize model and make fits image
# model = fitting.Model(observations=None, root='./', name=name)
# simple_disk = make_fits(model, disk_params, annulus_params)
# # fix_fits(model, starfluxes[2])
# # print(json.dumps(param_dict, indent=4))
#         # model.obs_sample(obs)
#         # model.get_chi(obs)
# 
# #Here's some hacky code for trying to make a contour plot
# #of the dust temperature and surface density
# #---------------------------------------------------------
# 
# # let's check out our continuum image now and see how the code did
# hdulist = fits.open(model.path + '.fits')
# disk_image = hdulist[0].data
# hdulist.close()
# 
# # disk_image[0,256,498]
# 
# plt.clf()
# plt.figure(figsize=(10,10))
# plt.imshow(disk_image[0,:,:], origin="lower", cmap="inferno",
#     norm=colors.LogNorm(vmin=disk_image.max()*1e-4, vmax=disk_image.max()))
# cbar = plt.colorbar()
# cbar.set_label(r"$I_{\nu}$ [ Jy pix$^{-1}$]")
# plt.xlabel("Pixel Column Number")
# plt.ylabel("Pixel Row Number")
# plt.title("AU Mic Model, $\Sigma$ F = %.3f Jy" % np.sum(disk_image[0,:,:] * 512.))
# plt.show()


#get dust temperature for every point in the disk grid
# temp_grid = simple_disk.temperature()
# #do the same for the dust density
# den_grid = simple_disk.dust_density()
# 
# plt.clf()
# plt.figure(figsize=(10,10))
# plot_temp = plt.contour(simple_disk.rf/simple_disk.AU, simple_disk.zf/simple_disk.AU, temp_grid.T, colors="k", linestyles=":")
# plot_den = plt.contour(simple_disk.rf/simple_disk.AU, simple_disk.zf/simple_disk.AU, den_grid.T, 200)
# # plt.clabel(plot_temp, fontsize=10, fmt="%.2f [K]")
# plt.colorbar(plot_den, label=r"$\rho_D$ [g cm$^{-3}$]")
# plt.title(r"AU Mic Model, $T$ and $\rho_D$")
# plt.xlabel("R [au]")
# plt.ylabel("Z [au]")
# plt.ylim(0, 15.)
# plt.xlim(5, 50)
# plt.show()
# 


# #initialize all AU Mic variables, beginning with the constants
# L_STAR_LSOL    = 0.09 # - parent star luminosity (in L_solar)
# STAR_MASS_MSOL = 0.31 # - parent star mass (in solar masses)
# T_INDEX        = -0.5 # - qq parameter
# 
# # NOTE: Here are the free parameters we should be testing
# DISK_MASS_MSOL = 3.68e-8 # - disk mass (in solar masses)
# RIN_AU         = 20     # - disk inner radius (in au)
# ROUT_AU        = 40    # - disk outer radius (in au)
# INCLINATION    = 88.5     # - disk inclination angle (in deg)
# P_INDEX        = 2.3     # - radial power law index (pp, for density)
# SH_PARAM       = .025      # - scale height parameter, h (numerical constant, H(R) = h * R)
# 
# #initialize model parameters (size of the cylindrical grid)
# R_GRIDSIZE = 500 #number of gridpoints in the radial direction
# Z_GRIDSIZE = 500 #number of gridpoints quitin the vertical direction
# 
# simple_disk = debris_disk_python3.Disk(params=[T_INDEX,DISK_MASS_MSOL,P_INDEX,RIN_AU,ROUT_AU,150.,
#                                INCLINATION,STAR_MASS_MSOL,1e-4,0.081,70.,[.79,1000],[50,500],
#                                -1, R_GRIDSIZE, Z_GRIDSIZE, L_STAR_LSOL, SH_PARAM],
#                                obs=[300, 351, 300, 5], annulus=[6,7.5,10e-5])
# 
# #run the radiative transfer and make a continuum image of the dust
# rt.total_model(simple_disk, distance=9.91, imres=0.03, xnpix=512, freq0=230.5, PA=128.41, 
# 	           offs=[0.0,0.0], nchans=1, modfile=save_path+save_root, isgas=False,
# 	           includeDust=True, extra=0.)
