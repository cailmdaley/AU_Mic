import numpy as np
import subprocess as sp
from astropy.io import fits
from collections import OrderedDict

from astrocail import fitting
from disk_model import debris_disk, raytrace
from aumic_observations import band6_observations


# default parameter dict:
param_dict = OrderedDict([
    ('temp_index',        -0.5),
    ('m_disk',            3.67e-08),
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
    model_fits = fits.open(model.root + model.name + '.fits')
    crpix = int(model_fits[0].header['CRPIX1'])
    model_im = model_fits[0].data[0]
    model_im[crpix, crpix] = starflux
    
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
    model = fitting.Model(observations=band6_observations,
        root=run_name + '/model_files/', 
        name='model' + str(np.random.randint(1e10)))
    make_fits(model, disk_params)
    for obs in band6_observations.flatten():
        fix_fits(model, obs, starflux)
        model.obs_sample(obs)
        model.get_chi(obs)
    
    model.delete()
    return -0.5 * sum(model.chis)


def make_best_fit(run):
    best_fit = run.chain[run.chain['lnprob'] == run.chain['lnprob'].max()]
        
    for param in best_fit.index[:-1]:
        param_dict[param] = best_fit[param]
        
    param_dict['m_disk'] = 10**param_dict['m_disk']
    param_dict['d_r'] += param_dict['r_in']

    disk_params = param_dict.values()[:-1]
    starflux = param_dict.values()[-1]
    
    # intialize model and make fits image 
    model = fitting.Model(observations=band6_observations,
        root=run.name + '/model_files/', 
        name='best_fit')
    model.delete()
    make_fits(model, disk_params)
    
    for pointing in band6_observations:
        ids = []
        for obs in pointing:
            fix_fits(model, obs, starflux)
            
            ids.append('_' + obs.name[12:20])
            model.obs_sample(obs, ids[-1])
            model.make_residuals(obs, ids[-1])
        # make model visibilities
        vis_files = ','.join([run.name+'/model_files/'+model.name+i+'.vis' for i in ids])
        sp.call(['uvcat', 'vis={}'.format(vis_files), 'out={}/model_files/best_fit_{}.vis'.format(run.name, obs.name[12:15])])
        
        # make residuals
        res_files = ','.join([run.name+'/model_files/'+model.name+i+'.residuals.vis' for i in ids])
        sp.call(['uvcat', 'vis={}'.format(res_files), 'out={}/model_files/best_fit_{}.residuals.vis'.format(run.name, obs.name[12:15])])
        
    # fig = plt.subplots(1, 2, sharex=False, sharey=False,
    #     figsize=(11.6, 6.5))[0]
    # plt.subplots_adjust(wspace=-0.0)

    # for pointing in model.observations:
    #     for i, obs in enumerate(pointing):
    #         model.obs_sample(obs)
    #         sp.call(['mv', model.path + '.vis', model.path + str(i) + '.vis'])
    #         
    #         
    # 
    # # cat_files = ','.join([run_name+ 'model_files/' + model.name'.vis' for obs in band6_observations[2]])
    # sp.call(['uvcat', 'vis={}'.format(cat_files), 'out=model_data/best_fit.vis'])
    # 
    # # Clean down to half the observation rms
    # sp.call(['invert',
    #     'vis=model_data/{}.vis'.format('best_fit'),
    #     'map=model_data/{}.mp'.format('best_fit'),
    #     'beam=model_data/{}.bm'.format('best_fit'),
    #     'cell=0.03arcsec', 'imsize=512', 'options=systemp,mfs', 'robust=2'])
    # sp.call(['clean',
    #     'map=model_data/{}.mp'.format('best_fit'),
    #     'beam=model_data/{}.bm'.format('best_fit'),
    #     'out=model_data/{}.cl'.format('best_fit'),
    #     'niters=100000', 'cutoff={}'.format(obs.rms/2)])
    # sp.call(['restor',
    #     'map=model_data/{}.mp'.format('best_fit'),
    #     'beam=model_data/{}.bm'.format('best_fit'),
    #     'model=model_data/{}.cl'.format('best_fit'),
    #     'out=model_data/{}.cm'.format('best_fit')])
    # sp.call(['fits', 'op=xyout',
    #     'in=model_data/{}.fits'.format('best_fit'),
    #     'out=model_data/{}.im'.format('best_fit')], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
