from astropy.io import fits
from glob import glob
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
plt.close()
sns.set_style('ticks')

def uvplot(uvf, bins=20, binning='lin', uv_bounds=None, label=None):
    
    dfits = fits.open(uvf)
    
    #Reading in polarized real & imaginary visibilities & weights (assuming ALMA data)
    rlimwt = dfits[0].data['data']
    rlimwt = np.mean(rlimwt, axis=3)
    
    
    rlxx = rlimwt[:,0,0,:,0,0]
    print(rlxx.shape)
    rlyy = rlimwt[:,0,0,:,1,0]
    imxx = rlimwt[:,0,0,:,0,1]
    imyy = rlimwt[:,0,0,:,1,1]
    wtxx = rlimwt[:,0,0,:,0,2]
    wtyy = rlimwt[:,0,0,:,1,2]
    
    u_meters = dfits[0].data['UU']
    v_meters = dfits[0].data['VV']
    
    pa = 0.0
    inc = np.radians(85.0)
    u1 = (u_meters*np.cos(pa)-v_meters*np.sin(pa))*np.cos(inc)
    v1 = (u_meters*np.sin(pa)+v_meters*np.cos(pa))
    #Rotation of (u,v) coordinates back to sky PA
    u2 = u1*np.cos(-1.0*pa)-v1*np.sin(-1.0*pa)
    v2 = u1*np.sin(-1.0*pa)+v1*np.cos(-1.0*pa)

    

    #Taking weighted average of polarized visibilities & weights
    rlcor = np.array(0.25*(1.0/wtxx + 1.0/wtyy)*(rlxx*wtxx + rlyy*wtyy))
    imcor = np.array(0.25*(1.0/wtxx + 1.0/wtyy)*(imxx*wtxx + imyy*wtyy))
    wtcor = 4.0/(1.0/wtxx + 1.0/wtyy)

    #Converting to kilolambda and Jy*km/s
    freq = dfits[0].header['CRVAL4']
    u_klam = u2*freq*1e-03
    v_klam = v2*freq*1e-03
    
    #Radial Averaging
    uv_dist = (u_klam**2 + v_klam**2)**0.5
    print(len(np.where((uv_dist > 23) & (uv_dist < 312))[0]))
    
    if not uv_bounds:
        uv_bounds = (np.min(uv_dist), np.max(uv_dist))
        print(np.min(uv_dist), np.max(uv_dist))
        
    if binning == 'lin':
        binbounds = np.linspace(uv_bounds[0], uv_bounds[1], bins+1)
    elif binning == 'log':
        binbounds = np.logspace(np.log10(uv_bounds[0]), np.log10(uv_bounds[1]), bins+1)
            
    binctr = np.zeros(bins)
    rlavg = np.zeros(bins)
    imavg = np.zeros(bins)
    rlwtedSEM = np.zeros(bins)
    imwtedSEM = np.zeros(bins)
    N = np.zeros(bins)
    rlxbar = np.zeros(bins)
    rlterm1 = np.zeros(bins)
    rlterm2 = np.zeros(bins)
    rlterm3 = np.zeros(bins)
    imxbar = np.zeros(bins)
    imterm1 = np.zeros(bins)
    imterm2 = np.zeros(bins)
    imterm3 = np.zeros(bins)
    
    for i in range(bins):
        binctr[i] = (binbounds[i] + binbounds[i+1])/2.0
        
        rlavg[i] = (np.nansum((rlcor*wtcor)[np.where( (uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]) )]) / (np.nansum(wtcor[np.where( (uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]) )])))*1e03
        imavg[i] = (np.nansum((imcor*wtcor)[np.where( (uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]) )]) / (np.nansum(wtcor[np.where( (uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]) )])))*1e03
        

        #Standard Error of the Mean
        N[i] = np.size(np.where( (uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]) ))
        rlxbar[i] = np.nansum((rlcor*wtcor)[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])/np.nansum(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])
        imxbar[i] = np.nansum((imcor*wtcor)[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])/np.nansum(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])
        rlterm1[i] = np.nansum(((rlcor*wtcor)[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])*rlxbar[i])**2)
        rlterm2[i] = -2*rlxbar[i]*np.nansum((wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]))*((wtcor*rlcor)[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])*rlxbar[i]))
        rlterm3[i] = rlxbar[i]**2*np.nansum((wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]))**2)
        imterm1[i] = np.nansum(((imcor*wtcor)[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])*imxbar[i])**2)
        imterm2[i] = -2*imxbar[i]*np.nansum((wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]))*((wtcor*imcor)[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])*imxbar[i]))
        imterm3[i] = imxbar[i]**2*np.nansum((wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]))**2)
        
        rlwtedSEM[i] = np.sqrt((N[i]/((N[i]-1)*np.nansum(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])**2))*(rlterm1[i]+rlterm2[i]+rlterm3[i]))
        imwtedSEM[i] = np.sqrt((N[i]/((N[i]-1)*np.nansum(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])**2))*(imterm1[i]+imterm2[i]+imterm3[i]))
        
    rlwtedSEM = rlwtedSEM*1e03
    imstedSEM = imwtedSEM*1e03


    #Plotting: setting up gridspec, setting ticks and labels, etc.
    rlax.set_xlim(uv_bounds)
    imax.set_xlabel(r'$\mathcal{R}_{uv}$ ($k\lambda$)', fontsize=15, labelpad=10)
    
    rlax.axhline(c='k', ls='--', linewidth=1)
    rlax.errorbar(binctr, rlavg, yerr=rlwtedSEM, fmt=':o', markersize=5)
    rlax.set_ylabel('Re (mJy)', fontsize=15)


    imax.errorbar( binctr, imavg, yerr=imwtedSEM, fmt=':o', markersize=5, label=label)
    imax.axhline(c='k', ls='--', linewidth=1)
    imax.set_ylim(rlax.get_ylim())
    imax.set_ylabel('Im (mJy)', fontsize=15)
    imax.legend()
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
fig, (rlax, imax) = plt.subplots(2, 1, sharey=False, sharex=True)
files = glob("*FINAL.uvf")
labels = ['mar', 'aug', 'jun']
savefig = None

for uvf, label in zip(files, labels):
    uvplot(uvf, bins = 6, binning='lin', uv_bounds=(23, 100), label=label)

    
    
# if savefig:
#     plt.savefig(''+savefig+'.eps')
plt.show(block=False)
