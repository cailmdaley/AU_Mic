                                    #Modeling Code
from radmc3dPy import *
from astropy.io import fits
import os
import numpy as np
import matplotlib.pylab as plb

#==============================================================================#

    #read in & print model parameters
analyze.writeDefaultParfile('ppdisk')
setup.problemSetupDust('ppdisk', crd_sys='sph', nx=[20], ny=[20], nz=[20], xbound='[0.1*au,100*au]', ybound='[0.0,pi]', zbound='[0.0,2*pi]', nw=[10], wbound=['0.1', '1500'], bgdens=1e-30, dusttogas=0.01, gasspec_vturb=1000, gsmax='10', gsmin='1', hpr_prim_rout=0.0, hrdisk=0.1, hrpivot='100*au', itempdecoup=0, lines_mode=0, mdisk=1e28, mstar='0.48*ms', nphot=1e6, plh=1/.7, plsig1=2.8, rin='8.8*au', rdisk='40.3*au', rstar='0.86*rs', sig0=0.0, srim_plsig=0.0, srim_rout=0.0, tstar=3600.)
par=analyze.readParams()
par.printPar()

#Create Model
os.system('radmc3d mctherm')
image.makeImage(npix=300., wav=1351.0, incl=89.5, posang=128.41+90, sizeau=150.)
im = image.readImage()

#Variables
modelname='aumic_model'
os.system('rm -rf '+modelname+'*')
vis = ['18aug2015_aumic_spw0.vis','18aug2015_aumic_spw1.vis','18aug2015_aumic_spw2.vis','18aug2015_aumic_spw3.vis','24jun2015_aumic1_spw0.vis','24jun2015_aumic1_spw1.vis','24jun2015_aumic1_spw2.vis','24jun2015_aumic1_spw3.vis','26mar2014_aumic_spw0.vis','26mar2014_aumic_spw1.vis','26mar2014_aumic_spw2.vis','26mar2014_aumic_spw3.vis']
uvf = ['18aug2015_aumic_spw0.uvf','18aug2015_aumic_spw1.uvf','18aug2015_aumic_spw2.uvf','18aug2015_aumic_spw3.uvf','24jun2015_aumic1_spw0.uvf','24jun2015_aumic1_spw1.uvf','24jun2015_aumic1_spw2.uvf','24jun2015_aumic1_spw3.uvf','26mar2014_aumic_spw0.uvf','26mar2014_aumic_spw1.uvf','26mar2014_aumic_spw2.uvf','26mar2014_aumic_spw3.uvf']
weightfile = ['18aug2015_aumic_spw0_weights.fits','18aug2015_aumic_spw1_weights.fits','18aug2015_aumic_spw2_weights.fits','18aug2015_aumic_spw3_weights.fits','24jun2015_aumic1_spw0_weights.fits','24jun2015_aumic1_spw1_weights.fits','24jun2015_aumic1_spw2_weights.fits','24jun2015_aumic1_spw3_weights.fits','26mar2014_aumic_spw0_weights.fits','26mar2014_aumic_spw1_weights.fits','26mar2014_aumic_spw2_weights.fits','26mar2014_aumic_spw3_weights.fits']
# 18aug2014 phasecenter='J2000 20h45m09.854710s -031d20m32.52034s'
# 24jun2015 phasecenter='J2000 20h45m09.867700s -31d20m32.89000s'
# 26mar2014 phasecenter='J2000 20h45m09.844300s -031d20m32.36000s'
coord=['20h45m09.854710s -031d20m32.52034s','20h45m09.854710s -031d20m32.52034s','20h45m09.854710s -031d20m32.52034s','20h45m09.854710s -031d20m32.52034s','20h45m09.867700s -031d20m32.89000s','20h45m09.867700s -031d20m32.89000s','20h45m09.867700s -031d20m32.89000s','20h45m09.867700s -031d20m32.89000s','20h45m09.844300s -031d20m31.36000s','20h45m09.844300s -031d20m31.36000s','20h45m09.844300s -031d20m31.36000s','20h45m09.844300s -031d20m31.36000s']
ra = [15.*(20.+45./60.+9.85471/3600.),15.*(20.+45./60.+9.85471/3600.),15.*(20.+45./60.+9.85471/3600.),15.*(20.+45./60.+9.85471/3600.),15.*(20.+45./60.+9.867700/3600.),15.*(20.+45./60.+9.867700/3600.),15.*(20.+45./60.+9.867700/3600.),15.*(20.+45./60.+9.867700/3600.),15.*(20.+45./60.+9.844300/3600.),15.*(20.+45./60.+9.844300/3600.),15.*(20.+45./60.+9.844300/3600.),15.*(20.+45./60.+9.844300/3600.)]
dec = [-31.-20./60.-32.52034/3600.,-31.-20./60.-32.52034/3600.,-31.-20./60.-32.52034/3600.,-31.-20./60.-32.52034/3600.,-31.-20./60.-32.89/3600.,-31.-20./60.-32.89/3600.,-31.-20./60.-32.89/3600.,-31.-20./60.-32.89/3600.,-31.-20./60.-32.36/3600.,-31.-20./60.-32.36/3600.,-31.-20./60.-32.36/3600.,-31.-20./60.-32.36/3600.]
starflux = 9.417385e-05
pixsize = '0.03arcsec'
model,data,datrlimwt,modrlimwt,openweights,weights,datrlxx,datrlyy,datimxx,datimyy,datrlI,datimI,modrlI,modimI,chi,redchi = [None]*len(vis),[None]*len(vis),[None]*len(vis),[None]*len(vis),[None]*len(vis),[None]*len(vis),[None]*len(vis),[None]*len(vis),[None]*len(vis),[None]*len(vis),[None]*len(vis),[None]*len(vis),[None]*len(vis),[None]*len(vis),[None]*len(vis),[None]*len(vis)


#Add star position & subtract star flux, create convolved model visibilities
for i in range(len(model)):
    im.writeFits(''+modelname+''+str(i)+'.fits', dpc=8.9, coord=''+coord[i]+'')
    model[i] = fits.open(''+modelname+''+str(i)+'.fits')
    model[i][0].data[0,150,150] = model[i][0].data[0,150,150] - starflux
    model[i][0].header['CRVAL1'] = ra[i]; model[i][0].header['CRVAL2'] = dec[i]
    os.system('rm -rf '+modelname+''+str(i)+'*')
    model[i].writeto(''+modelname+''+str(i)+'.fits')
    model[i].close
    os.system('fits in='+modelname+''+str(i)+'.fits op=xyin out='+modelname+''+str(i)+'.im')
    os.system('uvmodel model='+modelname+''+str(i)+'.im vis='+vis[i]+' options=replace out='+modelname+''+str(i)+'.vis')
    os.system('fits in='+modelname+''+str(i)+'.vis op=uvout out='+modelname+''+str(i)+'.uvf')

#Imaging convolved model & residuals
    os.system('uvmodel model='+modelname+''+str(i)+'.im vis='+vis[i]+' options=subtract out='+modelname+''+str(i)+'_residuals.vis')
    os.system('invert vis='+modelname+''+str(i)+'_residuals.vis map='+modelname+''+str(i)+'_residuals.mp beam='+modelname+''+str(i)+'_residuals.bm cell='+pixsize+' imsize=512 options=systemp,mfs robust=2')
    os.system('clean map='+modelname+''+str(i)+'_residuals.mp beam='+modelname+''+str(i)+'_residuals.bm out='+modelname+''+str(i)+'_residuals.cl  niters=10000') #cutoff=6.290E-05
    os.system('restor map='+modelname+''+str(i)+'_residuals.mp beam='+modelname+''+str(i)+'_residuals.bm model='+modelname+''+str(i)+'_residuals.cl out='+modelname+''+str(i)+'_residuals.cm')

    os.system('invert vis='+modelname+''+str(i)+'.vis map='+modelname+''+str(i)+'.mp beam='+modelname+''+str(i)+'.bm cell='+pixsize+' imsize=512 options=systemp,mfs robust=2')
    os.system('clean map='+modelname+''+str(i)+'.mp beam='+modelname+''+str(i)+'.bm out='+modelname+''+str(i)+'.cl  niters=10000') #cutoff=6.290E-05
    os.system('restor map='+modelname+''+str(i)+'.mp beam='+modelname+''+str(i)+'.bm model='+modelname+''+str(i)+'.cl out='+modelname+''+str(i)+'.cm')
    os.system('cgdisp in='+modelname+''+str(i)+'.cm device=/xs labtyp=arcsec beamtyp=b,l,3')

#Read in visibilities, convert data from polarized to Stokes I
for i in range(len(model)):
    data[i] = fits.open(''+uvf[i]+'')
    datrlimwt[i] = data[i][0].data['data']
    datrlxx[i] = datrlimwt[i][:,0,0,0,0,0,0]; datrlyy[i] = datrlimwt[i][:,0,0,0,0,1,0]
    datimxx[i] = datrlimwt[i][:,0,0,0,0,0,1]; datimyy[i] = datrlimwt[i][:,0,0,0,0,1,1]
    datrlI[i] = np.array((datrlxx[i]+datrlyy[i])/2.)
    datimI[i] = np.array((datimxx[i]+datimyy[i])/2.)

    model[i] = fits.open(''+modelname+''+str(i)+'.uvf')
    modrlimwt[i] = model[i][0].data['data']
    modrlI[i] = modrlimwt[i][::2,0,0,0,0,0]; modimI[i] = modrlimwt[i][::2,0,0,0,0,1]

    openweights[i] = fits.open(''+weightfile[i]+'')
    weights[i] = np.reshape(openweights[i][0].data, len(datrlxx[i]))

    data[i].close; model[i].close; openweights[i].close

#Calculate chi^2
    chi[i] = np.sum((datrlI[i] - modrlI[i])**2*weights[i] + (datimI[i] - modimI[i])**2*weights[i])
    redchi[i] = chi[i]/len(datrlI[i])
redchisum=np.sum(redchi)



