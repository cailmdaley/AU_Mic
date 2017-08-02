from astrocail.fitting import Observation

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

band6_observations=[mar0, mar1, mar2, mar3,
                    aug0, aug1, aug2, aug3,
                    jun0, jun1, jun2, jun3]
