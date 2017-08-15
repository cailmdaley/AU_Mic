import numpy as np
import argparse
import subprocess as sp; import os
from astropy.io import fits
from collections import OrderedDict

from astrocail import fitting, plotting, mcmc
from disk_model import debris_disk, raytrace
from aumic_observations import band6_observations, band6_rms_values, band6_fits_images

def main():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter, description= '''Python commands associated with emcee run8, which has 50 walkers and varies the following parameters:
    1)  disk mass
    2)  surface brightness power law exponent
    3)  inner radius
    4)  outer radius (really inner radius + dr)
    5)  inclination
    6)  position angle
    7)  march starflux
    8)  august starflux
    9)  june starflux
This run fixes the scale factor at 0.03, in order to investigate how sensitive our modeling is to the puffiness of the disk. If such a low scale factor yields significantly worse chi squareds and residuals, we can conclude that we are able to detect differences and scale factor and can thus measure it. However, if the chi squareds and residuals do not worsen very much, this indicates we do not have the sensitivity to measure the scale factor, and can only report an upper limit.''')
    
    parser.add_argument('-r', '--run', action='store_true', 
        help='begin or resume eemcee run.')
        
    parser.add_argument('-a', '--analyze', action='store_true', 
        help='analyze sampler chain, producing an evolution plot, corner plot, and image domain figure.')
    parser.add_argument('-b', '--burn_in', default=0,
        help='number of steps \'burn in\' steps to exclude')
    parser.add_argument('-c', '--corner', action='store_true', 
        help='generate corner plot')
    parser.add_argument('-e', '--evolution', action='store_true', 
        help='generate walker evolution plot.')
    parser.add_argument('-kde', '--density', action='store_true', 
        help='generate kernel density estimate (kde) of posterior distribution')
    args=parser.parse_args()
    
    
    if args.run:
        mcmc.run_emcee(run_name='run9', nsteps=10000, nwalkers=50,
            lnprob = lnprob, to_vary = [
            ('m_disk',            -7.55,        0.05,      (-np.inf, np.inf)),
            ('sb_law',            1.5,          1.5,       (-5.,     10.)),
            ('scale_factor',      0.2,          0.2,       (0,       np.inf)),
            ('r_in',              10,           10,        (0,       np.inf)),
            ('d_r',               31.5,         10,        (0,       np.inf)),
            ('inc',               89,           1,         (0,       90)),
            ('pa',                128.48,       0.1,       (1,       360)),
            ('mar_starflux',      3.94e-4,      5e-5,      (0,       np.inf)),
            ('aug_starflux',      1.50e-4,      1e-5,      (0,       np.inf)),
            ('jun_starflux',      2.30e-4,      1e-5,      (0,       np.inf))])
    else:
        run = mcmc.MCMCrun('run9', nwalkers=50, burn_in=args.burn_in)
        
    if args.analyze:
        make_best_fits(run)
        run.evolution()
        run.kde()
        run.corner()
        
    if args.evolution:
        run.evolution()
    if args.density:
        run.kde()
    if args.corner:
        run.corner()

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
    ('scale_factor',      0.001),
    ('pa',                128.41),
    ('mar_starflux',      2.50e-4),
    ('aug_starflux',      2.50e-4),
    ('jun_starflux',      2.50e-4)])

def make_fits(model, disk_params):
    structure_params = disk_params[:-1]
    PA = disk_params[-1]

    model_disk = debris_disk.Disk(structure_params, obs=[300, 131, 300, 20], sh_relation='const')
    raytrace.total_model(model_disk,
        distance=9.91, # pc
        imres=0.03, # arcsec/pix
        xnpix=512, #image size in pixels
        freq0=model.observations[0][0].uvf[0].header['CRVAL4']*1e-9, # obs frequency
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

    disk_params = param_dict.values()[:-3]
    starfluxes = param_dict.values()[-3:]

    # intialize model and make fits image
    model = fitting.Model(observations=band6_observations,
        root=run_name + '/model_files/',
        name='model' + str(np.random.randint(1e10)))
    make_fits(model, disk_params)
    for pointing, starflux in zip(model.observations, starfluxes):
        for obs in pointing:
            fix_fits(model, obs, starflux)
            model.obs_sample(obs)
            model.get_chi(obs)

    model.delete()

    return -0.5 * sum(model.chis)


def make_best_fits(run):
<<<<<<< HEAD
    # subset_df = run.main[run.main['r_in'] < 15]
    subset_df = run.main
=======
    subset_df = run.main[run.main['r_in'] < 15]
    # subset_df = run.main
>>>>>>> 54909d30a650cdf10296fdd8dc4748c028066764
    model_params = subset_df[subset_df['lnprob'] == subset_df['lnprob'].max()].drop_duplicates() # best fit
    print('Model parameters:')
    print(model_params.to_string())
    print('')

    for param in model_params.columns[:-1]:
        param_dict[param] = model_params[param].values


    param_dict['m_disk'] = 10**param_dict['m_disk']
    param_dict['d_r'] += param_dict['r_in']

    disk_params = param_dict.values()[:-3]
    starfluxes = param_dict.values()[-3:]

    # intialize model and make fits image
    print('Making model...')
    model = fitting.Model(observations=band6_observations,
        root=run.name + '/model_files/',
        name=run.name + '_bestfit')
    make_fits(model, disk_params)

    print('Sampling and cleaning...')
    visibilities = []
    residuals = []
    for pointing, rms, starflux in zip(model.observations, band6_rms_values, starfluxes):
        ids = []
        for obs in pointing:
            fix_fits(model, obs, starflux)

            ids.append('_' + obs.name[12:20])
            model.obs_sample(obs, ids[-1])
            model.make_residuals(obs, ids[-1])

        # make model visibilities
        vis_files = ','.join([model.path+ident+'.vis' for ident in ids])
        visibilities.append('{}_{}'.format(model.path, obs.name[12:15]))
        sp.call(['uvcat', 'vis={}'.format(vis_files), 'out={}.vis'.format(visibilities[-1])], stdout=open(os.devnull, 'wb'))
        model.clean(visibilities[-1], rms, show=False)

        # make residuals
        res_files = ','.join([model.path+ident+'.residuals.vis' for ident in ids])
        residuals.append('{}_{}.residuals'.format(model.path, obs.name[12:15]))
        sp.call(['uvcat', 'vis={}'.format(res_files), 'out={}.vis'.format(residuals[-1])], stdout=open(os.devnull, 'wb'))
        model.clean(residuals[-1], rms, show=False)

    print('Making figure...')
    fig = plotting.Figure(
        paths=[[band6_fits_images[i], visibilities[i]+'.fits', residuals[i]+'.fits']
            for i in range(len(visibilities))],
        rmses=[3*[rms] for rms in band6_rms_values],
        texts=[[[[4.6, 4.0, date]], [[4.6, 4.0, 'rms={}'.format(np.round(rms*1e6))]], None]
            for date, rms in zip(['March', 'August', 'June'], band6_rms_values)],
<<<<<<< HEAD
        # savefile=run.name+'/run8_bestfit_small_r_in.pdf', title=r'Run 8 Best Fit Model & Residuals for $r_{in} < 15$')
        savefile=run.name+'/run8_bestfit_global.pdf', title=r'Run 8 Global Best Fit Model & Residuals')
=======
        savefile=run.name+'/run9_bestfit_small_r_in.pdf', title=r'Run 9 Best Fit Model & Residuals for $r_{in} < 15$')
        # savefile=run.name+'/run9_bestfit_global.pdf', title=r'Run 9 Global Best Fit Model & Residuals')
>>>>>>> 54909d30a650cdf10296fdd8dc4748c028066764


    
if __name__ == '__main__':
    main()
