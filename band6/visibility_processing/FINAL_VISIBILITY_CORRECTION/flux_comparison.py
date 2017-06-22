from astropy.io import fits
from glob import glob
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
plt.close()

#INPUT VARIABLES//Opening fits file

# fileinp = raw_input('If you would like to deproject observed data, please enter its filename here. Otherwise, please hit return: ')
def uvplot(uvf, bins=20, binning='lin', uv_bounds=None):
    
    dfits = fits.open(uvf)

    #Reading in polarized real & imaginary visibilities & weights (assuming ALMA data)
    rlimwt = dfits[0].data['data']
    
    rlxx = rlimwt[:,0,0,:,0,0]
    rlyy = rlimwt[:,0,0,:,1,0]
    imxx = rlimwt[:,0,0,:,0,1]
    imyy = rlimwt[:,0,0,:,1,1]
    wtxx = rlimwt[:,0,0,:,0,2]
    wtyy = rlimwt[:,0,0,:,1,2]
    
    u_meters = dfits[0].data['UU']
    v_meters = dfits[0].data['VV']
    

    #Taking weighted average of polarized visibilities & weights
    rlcor = np.array(0.25*(1.0/wtxx + 1.0/wtyy)*(rlxx*wtxx + rlyy*wtyy))
    imcor = np.array(0.25*(1.0/wtxx + 1.0/wtyy)*(imxx*wtxx + imyy*wtyy))
    wtcor = 4.0/(1.0/wtxx + 1.0/wtyy)

    #Converting to kilolambda and Jy*km/s
    freq = dfits[0].header['CRVAL4']
    u_klam = u_meters*freq*1e-03
    v_klam = v_meters*freq*1e-03
    
    #Radial Averaging
    uv_dist = (u_klam**2 + v_klam**2)**0.5
    
    if not uv_bounds:
        uv_bounds = (np.min(uv_dist), np.max(uv_dist))
        
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
        
        rlavg[i] = (np.sum((rlcor*wtcor)[np.where( (uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]) )]) / (np.sum(wtcor[np.where( (uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]) )])))*1e03
        imavg[i] = (np.sum((imcor*wtcor)[np.where( (uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]) )]) / (np.sum(wtcor[np.where( (uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]) )])))*1e03


        #Standard Error of the Mean
        N[i] = np.size(np.where( (uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]) ))
        rlxbar[i] = np.sum((rlcor*wtcor)[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])/np.sum(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])
        imxbar[i] = np.sum((imcor*wtcor)[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])/np.sum(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])
        rlterm1[i] = np.sum(((rlcor*wtcor)[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])*rlxbar[i])**2)
        rlterm2[i] = -2*rlxbar[i]*np.sum((wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]))*((wtcor*rlcor)[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])*rlxbar[i]))
        rlterm3[i] = rlxbar[i]**2*np.sum((wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]))**2)
        imterm1[i] = np.sum(((imcor*wtcor)[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])*imxbar[i])**2)
        imterm2[i] = -2*imxbar[i]*np.sum((wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]))*((wtcor*imcor)[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])*imxbar[i]))
        imterm3[i] = imxbar[i]**2*np.sum((wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]-np.mean(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))]))**2)
        rlwtedSEM[i] = np.sqrt((N[i]/((N[i]-1)*np.sum(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])**2))*(rlterm1[i]+rlterm2[i]+rlterm3[i]))
        imwtedSEM[i] = np.sqrt((N[i]/((N[i]-1)*np.sum(wtcor[np.where((uv_dist>=binbounds[i]) & (uv_dist<=binbounds[i+1]))])**2))*(imterm1[i]+imterm2[i]+imterm3[i]))
        
    rlwtedSEM = rlwtedSEM*1e03
    imstedSEM = imwtedSEM*1e03


    #Plotting: setting y scales to be the same for both plots
    rlmax = np.ceil(max(rlavg) + rlwtedSEM[np.where(rlavg==max(rlavg))])
    rlmin = np.floor(min(rlavg) - rlwtedSEM[np.where(rlavg==min(rlavg))])
    imrange = np.ceil(max(np.absolute(imavg)) + max(imwtedSEM))
    
    rlrange = rlmax - rlmin

    #Plotting: setting up gridspec, setting ticks and labels, etc.
    gs = gridspec.GridSpec(2, 1, height_ratios=[rlrange,2*imrange])

    
    fig.xlim(uv_bounds)
    fig.axhline(c='k', ls='--', linewidth=2)
    
    rlax.ylabel('Re (mJy)', fontsize=15)
    
    rlax.errorbar(binctr, rlavg, yerr=rlwtedSEM, marker='.', c='b', fmt=' ' )

    # plt.xticks(np.arange(, max(binbounds)+.0001, 20))

    imax.errorbar( binctr, imavg, yerr=imwtedSEM, marker='.', fmt=' ' )
    # plt.xlim(0, xmax)
    # imax.axhline(c='k', ls='--', linewidth=2)
    # plt.xticks(np.arange(0, max(binbounds)+.0001, 20))
    # plt.yticks(np.linspace(-1*imrange, imrange, 3))
    imax.xlabel(r'$\mathcal{R}_{uv}$ ($k\lambda$)', fontsize=15, labelpad=10)
    imax.ylabel('Im (mJy)', fontsize=15)
    
fig, (rlax, imax) = plt.subplots(1, 2, sharey=True)


files = glob("*FINAL.uvf")
savefig = None

uvplot(files[0])

    
    
    
    
# if savefig:
#     plt.savefig(''+savefig+'.eps')
plt.show()
