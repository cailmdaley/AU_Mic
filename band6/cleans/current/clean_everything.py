from astrocail import casa
from glob import glob
import numpy as np

#To concat measurement sets
# mses=glob('../../final_visibilities/band6_star*.ms')
# aug = mses[:4]; jun = mses[4:8]; mar = mses[8:12]
# casa.concat(mses, 'band6_star_all')
# for date, outfile in zip([aug, jun, mar, mses], to_clean):
#     print(date)
#     print(outfile[0])
#     casa.concat(date, outfile[0])



# rms values--either from imstat on clean images or from residuals
band6_star_mar_allspws_natural_rms = 2.58940304301e-05
band6_star_aug_allspws_natural_rms = 2.70553631515e-05
band6_star_jun_allspws_natural_rms = 2.09217450873e-05
band6_star_all_natural_rms         = 1.47444691066e-05


band6_star_mar_allspws_200klambda_rms = 2.81487308808e-05
band6_star_all_200klambda_rms         = 1.91655869923e-05 

            
# Tuples with cleaning info
natural_cleans = [
    ('band6_star_aug_allspws', 'aumic_rms.region', 'aumic_mask_aug_natural.region',
        'natural', None),
    ('band6_star_jun_allspws', 'aumic_rms.region', 'aumic_mask_jun_natural.region',
        'natural', None),
    ('band6_star_mar_allspws', 'aumic_rms.region', 'aumic_mask_mar_taper.region', 
        'natural', None),
    ('band6_star_all', band6_star_all_natural_rms, 'aumic_mask_all_natural.region', 
        'natural', None)]
taper_cleans = [
    ('band6_star_mar_allspws', 'aumic_rms.region', 'aumic_mask_mar_taper.region', 
        'natural', ['200klambda']),
    ('band6_star_all', 'aumic_rms.region', 'aumic_mask_all_taper.region', 
        'natural', ['200klambda'])]
all_cleans = [
    ('band6_star_all', 'aumic_rms.region', 'aumic_mask_all_natural.region', 
        'natural', None),
    ('band6_star_all', 'aumic_rms.region', 'aumic_mask_all_taper.region', 
        'natural', ['200klambda'])]
            
            
            
rms_strings = []
casa.obs_clean(*all_cleans[-1])
    # rms_strings.append(obs[0] + '_rms = ' + raw_input('please enter clean rms: '))
# print(np.array(rms_strings))
    
