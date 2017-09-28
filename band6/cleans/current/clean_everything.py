from astrocail import casa
from glob import glob
import numpy as np

natural_cleans = [
    # ('band6_star_aug_allspws', 'aumic_mask_aug_natural.region'),
    # ('band6_star_jun_allspws', 'aumic_mask_jun_natural.region'),
    # ('band6_star_mar_allspws', 'aumic_mask_mar_natural.region'),
    ('band6_star_all', 'aumic_rms.region', 'aumic_mask_all_natural.region', 
        'natural',  None)]
taper_cleans = [
    # ('band6_star_mar_allspws', 'aumic_mask_mar_taper.region'),
    ('band6_star_all', 'aumic_rms.region', 'aumic_mask_all_taper.region', 
        'natural', ['200klambda'])]
            
            
# mses=glob('../../final_visibilities/band6_star*.ms')
# aug = mses[:4]; jun = mses[4:8]; mar = mses[8:12]
# 
# for date, outfile in zip([aug, jun, mar, mses], to_clean):
#     print(date)
#     print(outfile[0])
#     casa.concat(date, outfile[0])

            
            
rms_strings = []
for obs in taper_cleans:
    casa.obs_clean(*obs)
    # rms_strings.append(obs[0] + '_rms = ' + raw_input('please enter clean rms: '))
    # casa.to_fits(obs[0])
print(np.array(rms_strings))
    
band6_star_aug_allspws_rms = 2.70553631515e-05
band6_star_jun_allspws_rms = 2.09217450873e-05
band6_star_mar_allspws_rms = 4.79727064596e-05

band6_star_all200klambda_rms = 1.91655869923e-05
band6_star_allnatural_rms = 1.47444691066e-05
