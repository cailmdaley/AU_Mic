import numpy as np
import argparse
import subprocess as sp; import os
from astropy.io import fits
from collections import OrderedDict
import matplotlib.pyplot as plt; plt.switch_backend('agg')
import copy


from astrocail import fitting, plotting, mcmc
from disk_model import debris_disk, raytrace
import aumic_fitting

run_name='run26'
def main():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter, description= '''Python commands associated with emcee run26, which has 50 walkers and varies the following parameters:
    1)  disk mass
    2)  surface brightness power law exponent
    3)  inner radius
    4)  outer radius (really inner radius + dr)
    5)  inclination
    6)  position angle
    7)  march starflux
    8)  august starflux
    9) june starflux
This run is a redo of the previous 'skinny disk' model run17, and the aspect ratio is fixed at 0.003. This run uses the new Gaia distance to AU Mic and a new gridding resolution.''')

    parser.add_argument('-r', '--run', action='store_true',
        help='begin or resume eemcee run.')

    parser.add_argument('-a', '--analyze', action='store_true',
        help='analyze sampler chain, producing an evolution plot, corner plot, and image domain figure.')
    parser.add_argument('-b', '--burn_in', default=0, type=int,
        help='number of steps \'burn in\' steps to exclude')
    parser.add_argument('-bf', '--best_fit', action='store_true',
        help='generate best fit model images and residuals')
    parser.add_argument('-con', '--concise', action='store_true',
        help='concise best fit')
    parser.add_argument('-c', '--corner', action='store_true',
        help='generate corner plot')
    parser.add_argument('-cvars', '--corner_vars', default=None, nargs='+',
        help='concise best fit')
    parser.add_argument('-e', '--evolution', action='store_true',
        help='generate walker evolution plot.')
    parser.add_argument('-kde', '--kernel_density', action='store_true',
        help='generate kernel density estimate (kde) of posterior distribution')
    args=parser.parse_args()

    if args.run:
        mcmc.run_emcee(run_name=run_name, nsteps=10000, nwalkers=50,
            lnprob = lnprob, to_vary = [
            ('m_disk',            -7.54,        0.06,      (-np.inf, np.inf)),
            ('sb_law',            1,            0.5,       (-np.inf, np.inf)),
            ('r_in',              18,           5,         (0,       np.inf)),
            ('d_r',               24.3,         5,         (0,       np.inf)),
            ('inc',               88.6,         0.4,       (0,       90.)),
            ('pa',                128.49,       0.07,      (0,       360)),
            ('mar_starflux',      3.9e-4,       0.2e-4,    (0,       np.inf)),
            ('aug_starflux',      1.5e-4,       0.2e-4,    (0,       np.inf)),
            ('jun_starflux',      2.2e-4,       0.2e-4,    (0,       np.inf))])
    else:
        run = mcmc.MCMCrun(run_name, nwalkers=50, burn_in=args.burn_in)
        # old_nsamples = run.groomed.shape[0]
        # run.groomed = run.groomed[run.groomed['r_in'] + run.groomed['d_r'] > 20]
        # print('{} samples removed.'.format(old_nsamples - run.groomed.shape[0]))
        
        if args.corner_vars:
            cols = list(run.groomed.columns)
            col_indices = [cols.index(col) for col in args.corner_vars]

        if args.analyze or args.best_fit: make_best_fits(run, concise=args.concise)
        aumic_fitting.label_fix(run)
        if args.corner_vars: args.corner_vars = run.groomed.columns[col_indices]
        
        if args.analyze or args.evolution: run.evolution()
        if args.analyze or args.kernel_density: run.kde()
        if args.analyze or args.corner: run.corner(variables=args.corner_vars)

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
    ('scale_factor',      0.003),
    ('pa',                128.49),
    ('mar_starflux',      3.90e-4),
    ('aug_starflux',      1.50e-4),
    ('jun_starflux',      2.20e-4)])

def make_fits(model, disk_params):
    structure_params = disk_params[:-1]
    PA = disk_params[-1]

    model_disk = debris_disk.Disk(structure_params, obs=[300, 351, 300, 5])
    raytrace.total_model(model_disk,
        distance=9.725, # pc
        imres=0.015, # arcsec/pix
        xnpix=1024, #image size in pixels
        freq0=model.observations[0][0].uvf[0].header['CRVAL4']*1e-9, # obs frequency
        PA=PA,
        offs=[0.0,0.0], # offset from image center
        nchans=1, # continum
        isgas=False, # continuum!
        includeDust=True, #continuuum!!
        extra=0, # ?
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

    param_dict_copy = copy.copy(param_dict)
    param_dict_copy['m_disk'] = 10**param_dict_copy['m_disk']
    param_dict_copy['d_r'] += param_dict_copy['r_in']

    disk_params = param_dict_copy.values()[:-3]
    starfluxes = param_dict_copy.values()[-3:]

    # intialize model and make fits image
    model = fitting.Model(observations=aumic_fitting.band6_observations,
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


def make_best_fits(run, concise=False):
    subset_df = run.main#[run.main['r_in'] < 15]
    model_params = subset_df[subset_df['lnprob'] == subset_df['lnprob'].max()].iloc[0] # best fit
    print('Model parameters:')
    print(model_params.to_string())
    print('')

    for param in model_params.index[:-1]:
        param_dict[param] = model_params[param]


    param_dict['m_disk'] = 10**param_dict['m_disk']
    param_dict['d_r'] += param_dict['r_in']

    disk_params = param_dict.values()[:-3]
    starfluxes = param_dict.values()[-3:]

    # intialize model and make fits image
    print('Making model...')
    model = fitting.Model(observations=aumic_fitting.band6_observations,
        root=run.name + '/model_files/',
        name=run.name + '_bestfit')
    make_fits(model, disk_params)

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
        plotter = plotting.Plotter(
            paths=[[
                aumic_fitting.band6_fits_images[-1],
                paths[-1] + '.fits',
                paths[-1] + '.residuals.fits']],
            rmses=[3*[aumic_fitting.band6_rms_values[-1]]],
            global_colorscale=False)
        plotter.plot_all()

        # Subplot labels
        for i, text in enumerate(['Data', 'Model', 'Residuals']):
            plotter.add_text(subplot_index=i, x=4.6, y=4.0, text=text)

        # Scale bar
        x = -3.015; y = -4.7
        plotter.axes[0].plot([x, x - 10/9.725], [y, y], ls='-', linewidth=2, color='k')
        plotter.add_text(subplot_index=0, x=x+0.32, y=y+0.15, text="10 au")

        # For residual subplot: plot disk PA/extent as dashed line 
        # and stellar location as star
        proj_angle = np.radians(model_params['pa'] + 90)
        r_out_arcsec = model_params[['r_in', 'd_r']].sum() / 9.725
        x_proj = r_out_arcsec * -np.cos( proj_angle ) # neg b/c RA is left-positive
        y_proj = r_out_arcsec *  np.sin( proj_angle ) 
        plotter.axes[-1].plot([-x_proj, x_proj], [-y_proj, y_proj], c='k', ls=':', lw=2)
        plotter.axes[-1].plot(0, 0, marker='*', markersize=7, markeredgewidth=1, color='k')
        
        plotter.save(run.name+'/' + run.name + '_bestfit_concise.pdf', dpi=1000)
        # plotter.show()
    
    # else:
    #     fig = plotting.Figure(layout=(4,3),
    #        paths=[[obs, path + '.fits', path + '.residuals.fits']
    #            for obs, path in zip(aumic_fitting.band6_fits_images, paths)],
    #        rmses=[3*[rms] for rms in aumic_fitting.band6_rms_values],
    #        texts=[
    #            [[[4.6, 4.0, date]],
    #            [[4.6, 4.0, 'rms={}'.format(np.round(rms*1e6))]],
    #            None]
    #            for date, rms in zip(['March', 'August', 'June', 'All'],
    #            aumic_fitting.band6_rms_values)
    #            ],
    #        title= run.name + r'Global Best Fit Model & Residuals',
    #       savefile=run.name+'/' + run.name + '_bestfit_global.pdf')
          

def sample_disk_flux(run, best_fit = False):
    if best_fit:
        subset_df = run.main#[run.main['r_in'] < 15]
        model_params = subset_df[subset_df['lnprob'] == subset_df['lnprob'].max()].drop_duplicates() # best fit
    else:
        model_params = run.main.sample()
    
    for param in model_params.columns[:-1]:
        param_dict[param] = model_params[param].values


    param_dict['m_disk'] = 10**param_dict['m_disk']
    param_dict['d_r'] += param_dict['r_in']

    disk_params = param_dict.values()[:-3]
    starfluxes = param_dict.values()[-3:]

    # intialize model and make fits image
    model = fitting.Model(observations=aumic_fitting.band6_observations,
        root=run.name + '/model_files/',
        name=run.name + '_bestfit')
    make_fits(model, disk_params)

    model_im = fits.open(model.path + '.fits')[0].data[0]
    disk_flux = model_im.sum()
    print(disk_flux)
    return disk_flux
    
def disk_flux_stats(n_samples):
    run = mcmc.MCMCrun(run_name, nwalkers=50, burn_in=-2000)
    
    disk_fluxes = [sample_disk_flux(run) for i in range(n_samples)]
    print(disk_fluxes)
    
    bf = sample_disk_flux(run, best_fit = True)
    print('Best-fit disk flux: {}'.format(bf))
    print('Standard deviation of {} samples: {}'.format(n_samples, np.std(disk_fluxes)))
    

if __name__ == '__main__':
    main()
    
