import numpy as np
import argparse
import subprocess as sp; import os
from astropy.io import fits
from collections import OrderedDict

from astrocail import fitting, plotting, mcmc
from disk_model import debris_disk, raytrace
import aumic_fitting


# default parameter dict:
param_dict = OrderedDict([
    ('temp_index',        -0.5),
    ('m_disk',            -8.),
    ('sb_law',            2.3),
    ('r_in',              8.8),
    ('d_r',                31.5),
    ('r_crit',            150.0),
    ('inc',               89.5),
    ('m_star',            0.31),
    ('co_frac',           0.0001),
    ('v_turb',            0.081),
    ('Zq',                70.0),
    ('column_densities', [0.79, 1000]),
    ('abundance_bounds', [50, 500]),
    ('hand',              -1),
    ('rgrid_size',        500),
    ('zgrid_size',        500),
    ('l_star',            0.09),
    ('scale_factor',      0.1),
    ('pa',                128.41),
    ('starflux',      2.50e-4)])

def main():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter, description= '''Python commands associated with emcee run7, which has 18 walkers and varies the following parameters:
    1) disk mass
    2) surface brightness power law exponent
    3) scale factor, multiplied by radius to get scale height
    4) inner radius
    5) outer radius (really inner radius + dr)
    6) inclination
    7) position angle
    8) star flux (a single value is used for all three observations)
    This run implements a single starflux across all three dates.''')
    
    parser.add_argument('-r', '--run', action='store_true', 
        help='begin or resume eemcee run.')
        
    parser.add_argument('-a', '--analyze', action='store_true', 
        help='analyze sampler chain, producing an evolution plot, corner plot, and image domain figure.')
    parser.add_argument('-b', '--burn_in', default=0, type=int, 
        help='number of steps \'burn in\' steps to exclude')
    parser.add_argument('-bf', '--best_fit', action='store_true', 
        help='generate best fit model images and residuals')
    parser.add_argument('-c', '--corner', action='store_true', 
        help='generate corner plot')
    parser.add_argument('-e', '--evolution', action='store_true', 
        help='generate walker evolution plot.')
    parser.add_argument('-kde', '--kernel_density', action='store_true', 
        help='generate kernel density estimate (kde) of posterior distribution')
    args=parser.parse_args()
    
    
    if args.run:
        mcmc.run_emcee(run_name='run7', nsteps=10000, nwalkers=18, 
            lnprob = lnprob, to_vary = [
            ('m_disk',        -7.55,         0.05,       (-np.inf, np.inf)),
            ('sb_law',        2.3,           2,          (-5.,     10.)), 
            ('scale_factor',  0.05,          0.03,       (0,       np.inf)),
            ('r_in',          8.8,           5,          (0,       np.inf)),
            ('d_r',           31.5,          10,         (0,       np.inf)),
            ('inc',           90,            1,          (0,       np.inf)),
            ('pa',            128.48,        0.1,        (0,       360)),
            ('starflux',      2.50e-4,       1e-4,      (0,       np.inf))])
    else:
        run = mcmc.MCMCrun('run7', nwalkers=18, burn_in = args.burn_in)
        
    if args.analyze:
        make_best_fits(run)
        
        aumic_fitting.label_fix(run)
        run.evolution()
        run.kde()
        run.corner()
        
        return 
        
        
    if args.best_fit:
        make_best_fits(run)
        
    aumic_fitting.label_fix(run)
    if args.evolution:
        run.evolution()
    if args.kernel_density:
        run.kde()
    if args.corner:
        run.corner()
    
    return
    
def make_fits(model, disk_params):
    structure_params = disk_params[:-1]
    PA = disk_params[-1]

    model_disk = debris_disk.Disk(structure_params, obs=[300, 131, 300, 20])
    raytrace.total_model(model_disk,
        distance=9.91, # pc
        imres=0.03, # arcsec/pix
        xnpix=512, #image size in pixels
        freq0=model.observations[0][0].uvf[0].header['CRVAL4']*1e-9, # obs frequeency
        PA=PA,
        offs=[0.0,0.0], # offset from image center
        nchans=1, # continum
        isgas=False, # continuum!
        includeDust=True, #continuuum!!
        extra=0.0, # ?
        modfile = model.root + model.name)

    
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
        
# define likelehood functions
def lnprob(theta, run_name, to_vary):
    """
    For each parameter to be varied, return -infinity if the proposed value lies
    outside the bounds for the parameter; if all parameters are ok, return 0.
    """
    
    # reassign parameters to vary to those input by emcee
    for i, free_param in enumerate(to_vary):
        lower_bound, upper_bound = free_param[-1]
        
        if lower_bound < theta[i] < upper_bound:
            param_dict[free_param[0]] = theta[i]
        else: return -np.inf

    param_dict['m_disk'] = 10**param_dict['m_disk']
    param_dict['d_r'] += param_dict['r_in']
    
    disk_params = param_dict.values()[:-1]
    starflux = param_dict.values()[-1]
    
    # intialize model and make fits image 
    model = fitting.Model(observations=aumic_fitting.band6_observations,
        root=run_name + '/model_files/', 
        name='model' + str(np.random.randint(1e10)))
    make_fits(model, disk_params)
    for obs in aumic_fitting.band6_observations.flatten():
        fix_fits(model, obs, starflux)
        model.obs_sample(obs)
        model.get_chi(obs)
    
    model.delete()
    return -0.5 * sum(model.chis)


def make_best_fits(run, concise=False):
    print('Starting to make model image and residuals...')
    # subset_df = run.main[run.main['r_in'] < 15]
    subset_df = run.main
    model_params = subset_df[subset_df['lnprob'] == subset_df['lnprob'].max()].drop_duplicates() # best fit
    print('Model parameters:')
    print(model_params.to_string())
    print('')
        
    for param in model_params.columns[:-1]:
        param_dict[param] = model_params[param].values
        
        
    param_dict['m_disk'] = 10**param_dict['m_disk']
    param_dict['d_r'] += param_dict['r_in']

    disk_params = param_dict.values()[:-1]
    starflux = param_dict.values()[-1]
    
    # intialize model and make fits image 
    print('Making model...')
    model = fitting.Model(observations=aumic_fitting.band6_observations,
        root=run.name + '/model_files/', 
        name=run.name + '_bestfit')
    make_fits(model, disk_params)
    
    print('Sampling and cleaning...')
    paths = []
    for pointing, rms in zip(aumic_fitting.band6_observations, aumic_fitting.band6_rms_values):
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
            savefile=run.name+'/' + run.name + '_bestfit_concise.pdf')
    else:
        fig = plotting.Figure(layout=(4,3),
           paths=[[obs, path + '.fits', path + '.residuals.fits']
               for obs, path in zip(aumic_fitting.band6_fits_images, paths)],
           rmses=[3*[rms] for rms in aumic_fitting.band6_rms_values],
           texts=[
               [[[4.6, 4.0, date]],
               [[4.6, 4.0, 'rms={}'.format(np.round(rms*1e6))]],
               None]
               for date, rms in zip(['March', 'August', 'June', 'All'],
               aumic_fitting.band6_rms_values)
               ],
           title= run.name + r'Global Best Fit Model & Residuals',
          savefile=run.name+'/' + run.name + '_bestfit_global.pdf')
        
if __name__ == '__main__':
    main()
