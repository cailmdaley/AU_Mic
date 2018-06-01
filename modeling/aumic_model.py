import numpy as np
import argparse
import subprocess as sp; import os
from astropy.io import fits
from collections import OrderedDict
import matplotlib.pyplot as plt 
import matplotlib.colors as colors
import copy

from astrocail import fitting, plotting, mcmc
from disk_model import debris_disk, raytrace
import aumic_fitting

#initialize an instance of class Disk with specified parameters
# theta = np.fromstring(string[3:-1], sep='   ')
# default parameter dict:
param_dict = OrderedDict([
    ('temp_index',        -0.5),
    ('m_disk',            -7.54), # log
    ('sb_law',             2),
    ('r_in',               20),
    ('d_r',                22),
    ('r_crit',            150.0),
    ('inc',               88.6),
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
    ('scale_factor',      0.031),
    ('pa',                128.49),
    ('mar_starflux',      3.90e-4),
    ('aug_starflux',      1.50e-4),
    ('jun_starflux',      2.20e-4)])
    

def make_fits(model, disk_params):
    structure_params = disk_params[:-1]
    PA = disk_params[-1]

    model_disk = debris_disk.Disk(
        structure_params,
        obs=[300, 351, 300, 5],
        annulus=None)
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
          
def lnprob(name):
    temp_param_dict = copy.copy(param_dict)
    temp_param_dict['m_disk'] = 10**temp_param_dict['m_disk']
    temp_param_dict['d_r'] += temp_param_dict['r_in']

    disk_params = temp_param_dict.values()[:-3]
    starfluxes = temp_param_dict.values()[-3:]

    # intialize model and make fits image
    model = fitting.Model(observations=aumic_fitting.band6_observations,
        root='model_files/',
        name=name)
    make_fits(model, disk_params)
    for pointing, starflux in zip(model.observations, starfluxes):
        for obs in pointing:
            fix_fits(model, obs, starflux)
            model.obs_sample(obs)
            model.get_chi(obs)

    # model.delete()

    return model
    
def lnprob_new(name):
    temp_param_dict = copy.copy(param_dict)
    temp_param_dict['m_disk'] = 10**temp_param_dict['m_disk']
    temp_param_dict['d_r'] += temp_param_dict['r_in']

    disk_params = temp_param_dict.values()[:-3]
    starfluxes = temp_param_dict.values()[-3:]

    # intialize model and make fits image
    model = fitting.Model(observations=aumic_fitting.band6_observations,
        root='model_files/',
        name=name)
    make_fits(model, disk_params)
    for pointing, starflux in zip(model.observations, starfluxes):
        for obs in pointing:
            fix_fits(model, obs, starflux)
            model.get_chi_galario(obs)

    # model.delete()

    return model

old_good = sum(lnprob('old_good').chis)
new_good = sum(lnprob_new('new_good').chis)
print('old_good:', old_good, 'new_good:', new_good)

param_dict['scale_factor'] = 1

old_bad = sum(lnprob('old_bad').chis)
new_bad = sum(lnprob_new('new_bad').chis)
print('old_bad:', old_bad, 'new_bad:', new_bad)
print('')
print('bad - good:')
print('old:', old_bad - old_good)
print('new:', new_bad - new_good)
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
