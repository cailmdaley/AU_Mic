from astrocail import casa
from glob import glob
import pandas as pd
import numpy as np
from astropy.io import fits

clean_dir = 'cleans/'
data_dir = '../visibility_processing/24jun2015_flare_reduced/'
ms_dir = 'measurement_sets/'
mask_path = '../cleans/current/aumic_mask_jun_natural.region'
rms_path = '../cleans/current/aumic_rms.region'

def concat_flare():
    # Concat phase-shifted June flare files into a single flare file
    mses=glob(data_dir + '*.fixvis')
    casa.concat(mses, 'band6_24jun_flare.fixvis.ms')
    
def concat_noflare():
    # Concat phase-shifted June flare files into a single flare file
    mses=glob('../final_visibilities/band6_star_jun*.ms')
    casa.concat(mses, ms_dir + 'band6_24jun_noflare.fixvis.ms')
    
def concat_mar_and_aug():
    # Concat phase-shifted June flare files into a single flare file
    for month in ['aug', 'mar']:
        mses=glob('../final_visibilities/band6_star_{}*.ms'.format(month))
        casa.concat(mses, 'band6_{}.fixvis.ms'.format(month))

def clean(filename): # Clean flare file 
    casa.obs_clean(
        filename=filename, 
        rms=rms_path,
        mask=mask_path,
        input_dir = ms_dir,
        output_dir='cleans/',
        datacolumn='data') # Clean rms is 4.33658331298e-05
        
def imfit_flare_location():
    casa.pipe([
        "from pprint import pprint",
        '''noflare = imfit( 
        imagename = '{}' + 'band6_24jun_noflare.fixvis.natural_clean.image/',
        region    = 'circle [ [256pix, 256pix], 20 pix]',
        rms       = 2.09217450873e-05,
        summary   = '{}' + 'band6_24jun_noflare.fixvis.natural_clean.imfit_summary')
        '''.format(clean_dir, clean_dir),
        
        "print('No flare:')",
        "pprint(noflare['deconvolved']['component0'])",
        
        '''flare = imfit( 
        imagename = '{}' + 'band6_24jun_flare.fixvis.natural_clean.image/',
        region    = 'circle [ [256pix, 256pix], 20 pix]',
        rms       = 4.47719756884e-05,
        summary   = '{}' + 'band6_24jun_flare.fixvis.natural_clean.imfit_summary')
        '''.format(clean_dir, clean_dir),
        
        "print('Flare:')",
        "pprint(flare['deconvolved']['component0'])"])
        
def clean_flare_timebins():
    # For timebin n (in range 3-9 inclusive):
    #   1. Concat the four spws ending with `n.fixvis`
    #   2. Clean resulting file.
    for i in range(3, 10):
        mses = glob(data_dir + '*{}.fixvis'.format(i))
        casa.concat(mses, ms_dir + 'band6_24jun_flare_{}.concat.ms'.format(i))
        casa.pipe(["fixvis(vis='{}band6_24jun_flare_{}.concat.ms', outputvis='{}band6_24jun_flare_{}.concat.fixvis.ms', phasecenter='J2000 -0.85013027rad -0.54702896rad')".format(ms_dir, i, ms_dir, i)])
        casa.obs_clean(
            filename='band6_24jun_flare_{}.concat.fixvis'.format(i), 
            rms=rms_path,
            mask=mask_path,
            input_dir=ms_dir,
            output_dir='cleans/',
            datacolumn='data')

# def imfit_flare_timebins():
#     rms_values=[0.000151195035669, 0.000120755589768, 0.000104318845558, 
#         0.000101897044979, 0.000107725165907, 9.96594866115e-05, 
#         0.000102482076988]
#     casa.pipe(['''imfit( 
#         imagename = '{}' + 'band6_24jun_flare_{}.fixvis.concat.natural_clean.image/',
#         region    = 'circle [ [20:45:09.867700,-31.20.32.89000], 0.5 arcsec ]',
#         rms       = {},
#         summary   = '{}' + 'band6_24jun_flare_{}.fixvis.concat.natural_clean.imfit_summary')
#         '''.format(clean_dir, i, rms_values[i-3], clean_dir, i)
#         for i in range(3,10)])
# 
# def calculate_flare_fluxes():
#     summary_file = clean_dir + 'band6_24jun_noflare.fixvis.natural_clean.imfit_summary'
#     imfit_stats = pd.read_table(summary_file, header=1, skiprows=[2], 
#         index_col=0, delim_whitespace=True, ).iloc[:,1:]
#     imfit_stats.loc['units'] = pd.read_table(summary_file, header=0, 
#         skipfooter=None, index_col=1, delim_whitespace=True).columns[1:]
#     imfit_stats.loc['no flare'] = np.loadtxt(summary_file)[1:]
# 
#     for i in range(3, 10):
#         summary_file=clean_dir + 'band6_24jun_flare_{}.fixvis.concat.natural_clean.imfit_summary'.format(i)
#         imfit_stats.loc[i] = np.loadtxt(summary_file)[1:]
#     return imfit_stats    
#     imfit_stats['I-disk'] = imfit_stats['I'] - 0.001
    
def imstat_flare_timebins():
    commands = [
        "foo = imstat('{}' + 'band6_24jun_noflare.fixvis.natural_clean.image', region='circle [ [256pix, 256pix], 20 pix]')".format(clean_dir),
        "stats_array = [ ['no_flare', [foo['maxpos'][0], foo['maxpos'][1]], foo['max'][0] ] ]"]
    for i in range(3,10): 
        commands.append("foo = imstat('{}' + 'band6_24jun_flare_{}.concat.fixvis.natural_clean.image/', region='circle [ [256pix, 256pix], 20 pix]')".format(clean_dir, i))
#         region='circle [ [20:45:09.867700,-31.20.32.89000], 0.5 arcsec ]',
        commands.append("stats_array.append( [{}, [foo['maxpos'][0], foo['maxpos'][1]], foo['max'][0]] )".format(i))
    commands.append("np.savetxt('flare_imstats.txt', np.array(stats_array, dtype='object'), fmt='%s', delimiter='\t')")
    casa.pipe(commands)
        
def calculate_flare_fluxes():
    df = pd.read_table('flare_imstats.txt', index_col=0, names=['pixel', 'flux'])
    df['rms'] = [2.09217450873e-05, 0.000151195035669, 0.000120755589768, 
        0.000104318845558, 0.000101897044979, 0.000107725165907, 
        9.96594866115e-05, 0.000102482076988]
    model_image = fits.getdata('../../modeling/run27/model_files/run27_bestfit_un_.fits').squeeze()
    for i in df.index:
        df.loc[i, 'flux'] -= model_image[tuple(eval(df.loc['3', 'pixel']))]
    df.loc[:, ['flux', 'rms']] *= 1e3 # mJy
    return df.round(2)
    
# # Concatenate and clean flare vs no-flare files; compare imfit locations
# concat_flare()
# clean('band6_24jun_flare.fixvis')
# 
# concat_noflare()
# clean('band6_24jun_noflare.fixvis')
# imfit_flare_location()

# # Image flare in each time bin, get brightest pixel, and subtract best-fit convolved disk flux in corresponding pixel
# clean_flare_timebins()
# imstat_flare_timebins()
# calculate_flare_fluxes()
