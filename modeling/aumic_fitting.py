from astrocail.fitting import Observation
import numpy as np

mar0 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_mar_spw0_FINAL', rms=6.5e-05)
mar1 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_mar_spw1_FINAL', rms=6.124e-05)
mar2 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_mar_spw2_FINAL', rms=6.068e-05)
mar3 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_mar_spw3_FINAL', rms=6.468e-05)
aug0 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_aug_spw0_FINAL', rms=5.879e-05)
aug1 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_aug_spw1_FINAL', rms=5.336e-05)
aug2 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_aug_spw2_FINAL', rms=6.092e-05)
aug3 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_aug_spw3_FINAL', rms=5.558e-05)
jun0 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_jun_spw0_FINAL', rms=5.369e-05)
jun1 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_jun_spw1_FINAL', rms=4.658e-05)
jun2 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_jun_spw2_FINAL', rms=5.083e-05)
jun3 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_jun_spw3_FINAL', rms=5.559e-05)

band6_observations=np.array(
    [[mar0, mar1, mar2, mar3],
    [aug0, aug1, aug2, aug3],
    [jun0, jun1, jun2, jun3]])
    
band6_rms_values = np.array([
2.58940304301e-05,
2.70553631515e-05,
2.09217450873e-05,
1.47444691066e-05])

band6_fits_images = [
    '../band6/cleans/current/band6_star_mar_allspws.natural_clean.fits', 
    '../band6/cleans/current/band6_star_aug_allspws.natural_clean.fits', 
    '../band6/cleans/current/band6_star_jun_allspws.natural_clean.fits',
    '../band6/cleans/current/band6_star_all.natural_clean.fits']
    

def label_fix(run):
    for df in [run.main, run.groomed]:
    
        df.loc[:,'d_r'] += df.loc[:, 'r_in']
        try:
            df.loc[:,'starflux'] *= 1e6
        except:
            df.loc[:,'mar_starflux'] *= 1e6
            df.loc[:,'jun_starflux'] *= 1e6
            df.loc[:,'aug_starflux'] *= 1e6
            
        df.loc[:,'inc'].where(df.loc[:,'inc'] < 90, 180-df.loc[:,'inc'], inplace=True)
        
        df.rename(inplace=True, columns={
            'm_disk' : r'Disk Mass ($M_{\odot}$)', 
            'sb_law' : r'SB Law',
            'scale_factor' : r'Scale Factor', 
            'r_in' : r'$r_{in}$ (au)',
            'd_r' : r'$r_{out}$ (au)',
            'inc' : r'$i$ ($\degree$)',
            'pa'  : r'PA  ($\degree$)',
            'mar_starflux' : r'March $F_{*}$ ($\mu$Jy)',
            'aug_starflux' : r'August $F_{*}$ ($\mu$Jy)',
            'jun_starflux' : r'June $F_{*}$ ($\mu$Jy)'})
    
