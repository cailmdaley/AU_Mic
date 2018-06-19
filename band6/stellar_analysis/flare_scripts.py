from astrocail import casa
from glob import glob
import pandas as pd
import numpy as np

clean_dir = 'cleans/'
data_dir = '../visibility_processing/24jun2015_flare_reduced/'
mask_path = '../cleans/current/aumic_mask_jun_natural.region'
rms_path = '../cleans/current/aumic_rms.region'

def concat_flare():
    # Concat phase-shifted June flare files into a single flare file
    mses=glob(data_dir + '*.fixvis')
    casa.concat(mses, 'band6_24jun_flare.fixvis.ms')
    
def concat_noflare():
    # Concat phase-shifted June flare files into a single flare file
    mses=glob('../final_visibilities/band6_star_jun*.ms')
    casa.concat(mses, 'band6_24jun_noflare.fixvis.ms')
    
def concat_mar_and_aug():
    # Concat phase-shifted June flare files into a single flare file
    for month in ['aug', 'mar']:
        mses=glob('../final_visibilities/band6_star_{}*.ms'.format(month))
        casa.concat(mses, 'band6_{}.fixvis.ms'.format(month))
concat_mar_and_aug()

def clean(filename): # Clean flare file 
    casa.obs_clean(
        filename=filename, 
        rms=rms_path,
        mask=mask_path,
        output_dir='cleans/',
        datacolumn='data') # Clean rms is 4.33658331298e-05
        
def imfit_flare_location():
    casa.pipe([
        '''noflare = imfit( 
        imagename = '{}' + 'band6_24jun_noflare.fixvis.natural_clean.image/',
        region    = 'circle [ [20:45:09.867700,-31.20.32.89000], 0.5 arcsec ]',
        rms       = 2.09217450873e-05,
        summary   = '{}' + 'band6_24jun_noflare.fixvis.natural_clean.imfit_summary')
        '''.format(clean_dir, clean_dir),
        
        "print('No flare:')",
        "print(noflare['deconvolved']['component0']['shape']['direction'])",
        
        '''flare = imfit( 
        imagename = '{}' + 'band6_24jun_flare.fixvis.natural_clean.image/',
        region    = 'circle [ [20:45:09.867700,-31.20.32.89000], 0.5 arcsec ]',
        rms       = 4.47719756884e-05,
        summary   = '{}' + 'band6_24jun_flare.fixvis.natural_clean.imfit_summary')
        '''.format(clean_dir, clean_dir),
        
        "print('Flare:')",
        "print(flare['deconvolved']['component0']['shape']['direction'])"])
        
        
def clean_flare_timebins():
    # For timebin n (in range 3-9 inclusive):
    #   1. Concat the four spws ending with `n.fixvis`
    #   2. Clean resulting file.
    for i in range(3, 10):
        mses = glob(data_dir + '*{}.fixvis'.format(i))
        casa.concat(mses, data_dir + 'band6_24jun_flare_{}.fixvis.concat.ms'.format(i))
        casa.obs_clean(
            filename='band6_24jun_flare_{}.fixvis.concat'.format(i), 
            rms=rms_path,
            mask=mask_path,
            input_dir=data_dir,
            output_dir='cleans/',
            datacolumn='data')

def imfit_flare_timebins():
    rms_values=[0.000151195035669, 0.000120755589768, 0.000104318845558, 
        0.000101897044979, 0.000107725165907, 9.96594866115e-05, 
        0.000102482076988]
    casa.pipe(['''imfit( 
        imagename = '{}' + 'band6_24jun_flare_{}.fixvis.concat.natural_clean.image/',
        region    = 'circle [ [20:45:09.867700,-31.20.32.89000], 0.5 arcsec ]',
        rms       = {},
        summary   = '{}' + 'band6_24jun_flare_{}.fixvis.concat.natural_clean.imfit_summary')
        '''.format(clean_dir, i, rms_values[i-3], clean_dir, i)
        for i in range(3,10))]
            
def get_imfit_stats():
    summary_file = clean_dir + 'band6_24jun_noflare.fixvis.natural_clean.imfit_summary'
    imfit_stats = pd.read_table(summary_file, header=1, skiprows=[2], 
        index_col=0, delim_whitespace=True, ).iloc[:,1:]
    imfit_stats.loc['units'] = pd.read_table(summary_file, header=0, 
        skipfooter=None, index_col=1, delim_whitespace=True).columns[1:]
    imfit_stats.loc['no flare'] = np.loadtxt(summary_file)[1:]
        
    for i in range(3, 10):
        summary_file=clean_dir + 'band6_24jun_flare_{}.fixvis.concat.natural_clean.imfit_summary'.format(i)
        imfit_stats.loc[i] = np.loadtxt(summary_file)[1:]
    return imfit_stats
    
# concat_flare()
# clean('band6_24jun_flare.fixvis')
# concat_noflare()
# clean('band6_24jun_noflare.fixvis')
# imfit_flare_location()



# clean_flare_timebins()
# imfit_flare_timebins()
get_imfit_stats()
