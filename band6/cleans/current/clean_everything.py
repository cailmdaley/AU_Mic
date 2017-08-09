from astrocail import casa
from glob import glob
import numpy as np

to_clean = [('band6_star_aug_allspws', 'aumic_mask_aug_natural.region'),
            ('band6_star_jun_allspws', 'aumic_mask_jun_natural.region'),
            ('band6_star_mar_allspws', 'aumic_mask_mar_natural.region'),
            ('band6_star_all',         'aumic_mask_all_natural.region')]
            
            
            
# mses=glob('../../final_visibilities/band6_star*.ms')
# aug = mses[:4]; jun = mses[4:8]; mar = mses[8:12]
# 
# for date, outfile in zip([aug, jun, mar, mses], to_clean):
#     print(date)
#     print(outfile[0])
#     casa.concat(date, outfile[0])

            
            
rms_strings = []
for obs in to_clean:
    casa.obs_clean(obs[0], 'aumic_rms.region', obs[1])
    rms_strings.append(obs[0] + '_rms = ' + raw_input('please enter clean rms: '))
print(np.array(rms_strings))
    
band6_star_aug_allspws_rms = 2.70553631515e-05
band6_star_jun_allspws_rms = 2.09217450873e-05
band6_star_mar_allspws_rms = 4.79727064596e-05
band6_star_all_rms = 1.47444691066e-05