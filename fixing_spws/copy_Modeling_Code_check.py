#Modeling Code
from radmc3dPy import *
from astropy.io import fits
import os
import numpy as np
import matplotlib.pylab as plb
#==============================================================================#



#Variables and input files:
modelname='aumic_model'
filenames = ['18aug2015_aumic_spw0.corrected_weights','18aug2015_aumic_spw1.corrected_weights','18aug2015_aumic_spw2.corrected_weights','18aug2015_aumic_spw3.corrected_weights','24jun2015_aumic1_spw0.corrected_weights','24jun2015_aumic1_spw1.corrected_weights','24jun2015_aumic1_spw2.corrected_weights','24jun2015_aumic1_spw3.corrected_weights','26mar2014_aumic_spw0.corrected_weights','26mar2014_aumic_spw1.corrected_weights','26mar2014_aumic_spw2.corrected_weights','26mar2014_aumic_spw3.corrected_weights']
coord=['20h45m09.854710s -031d20m32.52034s','20h45m09.854710s -031d20m32.52034s','20h45m09.854710s -031d20m32.52034s','20h45m09.854710s -031d20m32.52034s','20h45m09.867700s -031d20m32.89000s','20h45m09.867700s -031d20m32.89000s','20h45m09.867700s -031d20m32.89000s','20h45m09.867700s -031d20m32.89000s','20h45m09.844300s -031d20m31.36000s','20h45m09.844300s -031d20m31.36000s','20h45m09.844300s -031d20m31.36000s','20h45m09.844300s -031d20m31.36000s']
ra = [15.*(20.+45./60.+9.85471/3600.),15.*(20.+45./60.+9.85471/3600.),15.*(20.+45./60.+9.85471/3600.),15.*(20.+45./60.+9.85471/3600.),15.*(20.+45./60.+9.867700/3600.),15.*(20.+45./60.+9.867700/3600.),15.*(20.+45./60.+9.867700/3600.),15.*(20.+45./60.+9.867700/3600.),15.*(20.+45./60.+9.844300/3600.),15.*(20.+45./60.+9.844300/3600.),15.*(20.+45./60.+9.844300/3600.),15.*(20.+45./60.+9.844300/3600.)]
dec = [-31.-20./60.-32.52034/3600.,-31.-20./60.-32.52034/3600.,-31.-20./60.-32.52034/3600.,-31.-20./60.-32.52034/3600.,-31.-20./60.-32.89/3600.,-31.-20./60.-32.89/3600.,-31.-20./60.-32.89/3600.,-31.-20./60.-32.89/3600.,-31.-20./60.-32.36/3600.,-31.-20./60.-32.36/3600.,-31.-20./60.-32.36/3600.,-31.-20./60.-32.36/3600.]
# filenames = ["18aug2015_aumic_spw0.corrected_weights", "24jun2015_aumic1_spw0.corrected_weights", "24jun2015_aumic1_spw1.corrected_weights", "24jun2015_aumic1_spw2.corrected_weights", "24jun2015_aumic1_spw3.corrected_weights"]
# # 18aug2014 phasecenter='J2000 20h45m09.854710s -031d20m32.52034s'
# # 24jun2015 phasecenter='J2000 20h45m09.867700s -31d20m32.89000s'
# # 26mar2014 phasecenter='J2000 20h45m09.844300s -031d20m32.36000s'
# coord =['20h45m09.854710s -031d20m32.52034s', '20h45m09.867700s -031d20m32.89000s', '20h45m09.867700s -031d20m32.89000s', '20h45m09.867700s -031d20m32.89000s', '20h45m09.867700s -031d20m32.89000s']
# ra = [311.29106129166667, 311.2911154166666,  311.2911154166666, 311.2911154166666, 311.2911154166666]
# dec =[-31.34236676111111, -31.342469444444443, -31.342469444444443, -31.342469444444443, -31.342469444444443]
starflux = 9.417385e-05
pixsize = '0.03arcsec'

#Create Model
def create_model():
    """
    Returns a model image of the disk
    """

    analyze.writeDefaultParfile('ppdisk')
    setup.problemSetupDust('ppdisk', crd_sys='sph', nx=[20], ny=[20], nz=[20], xbound='[0.1*au,100*au]', ybound='[0.0,pi]', zbound='[0.0,2*pi]', nw=[10], wbound=['0.1', '1500'], bgdens=1e-30, dusttogas=0.01, gasspec_vturb=1000, gsmax='10', gsmin='1', hpr_prim_rout=0.0, hrdisk=0.1, hrpivot='100*au', itempdecoup=0, lines_mode=0, mdisk=1e28, mstar='0.31*ms', nphot=1e6, plh=1/.7, plsig1=2.8, rin='8.8*au', rdisk='40.3*au', rstar='0.86*rs', sig0=0.0, srim_plsig=0.0, srim_rout=0.0, tstar=3600.)

    par=analyze.readParams()
    par.printPar()

    os.system('radmc3d mctherm')
    image.makeImage(npix=300., wav=1351.0, incl=89.5, posang=218.41, sizeau=150.)
    im = image.readImage()

    return im

#Add star position & subtract star flux, create convolved model visibilities
def model_convolve(im, modelname, num, filename, coord, ra, dec):
    """
    Create model fits file with correct header information, remove stellar emission, convolve with ALMA visiblities to create model .vis and .uvf files.
    """

    os.system('rm -rf '+modelname+''+str(num)+'*')

    im.writeFits(''+modelname+''+str(num)+'.fits', dpc=8.9, coord=''+coord+'')
    model = fits.open(''+modelname+''+str(num)+'.fits')
    model[0].data[0,150,150] -= starflux
    model[0].header['CRVAL1'] = ra
    model[0].header['CRVAL2'] = dec

    #Clear pre-existing models before writing
    os.system('rm -rf '+modelname+''+str(num)+'.fits')
    model.writeto(''+modelname+''+str(num)+'.fits')
    model.close
    os.system('fits in='+modelname+''+str(num)+'.fits op=xyin out='+modelname+''+str(num)+'.im')
    os.system('uvmodel model='+modelname+''+str(num)+'.im vis='+filename+'.vis options=replace out='+modelname+''+str(num)+'.vis')
    os.system('fits in='+modelname+''+str(num)+'.vis op=uvout out='+modelname+''+str(num)+'.uvf')

def get_chi(filename, modelname, num):
    """
    Return chi^2 of model using seperate weightfile.
    """

    data = fits.open(''+filename+'.uvf')
    datrlimwt = data[0].data['data']

    #splitting visibilities by time removes one index from datrlimwt array: so, must check size of datrlimwt
    if len(datrlimwt.shape) == 7:
        datrlxx = datrlimwt[:,0,0,0,0,0,0]
        datrlyy = datrlimwt[:,0,0,0,0,1,0]
        datimxx = datrlimwt[:,0,0,0,0,0,1]
        datimyy = datrlimwt[:,0,0,0,0,1,1]
        weights = datrlimwt[:,0,0,0,0,0,2]
    elif len(datrlimwt.shape) == 6:
        datrlxx = datrlimwt[:,0,0,0,0,0]
        datrlyy = datrlimwt[:,0,0,0,1,0]
        datimxx = datrlimwt[:,0,0,0,0,1]
        datimyy = datrlimwt[:,0,0,0,1,1]
        weights = datrlimwt[:,0,0,0,0,2]
    datrlI = np.array((datrlxx+datrlyy)/2.)
    datimI = np.array((datimxx+datimyy)/2.)

    model = fits.open(''+modelname+''+str(num)+'.uvf')
    modrlimwt = model[0].data['data']
    modrlI = modrlimwt[::2,0,0,0,0,0]
    modimI = modrlimwt[::2,0,0,0,0,1]


    data.close; model.close


#Calculate chi^2
    chi = np.sum((datrlI - modrlI)**2*weights + (datimI - modimI)**2*weights)
    redchi = chi/len(datrlI)

    return chi, redchi




#Make this part do the file splitting
#os.system("uvlist vis=24jun2015_aumic1_spw3.vis options=baseline")
#os.system("prthd in=24jun2015_aumic1_spw3.vis")
def vis_split(vis, select, suff, filenames):
    os.system("rm -r {}.vis".format(vis+suff))
    os.system("uvaver vis={}.vis select={} out={}.vis".format(vis, select, vis+suff))
    filenames.append(vis+suff)

newfiles = []
for i in [4,5,6,7]:
    print(filenames[i])
    vis_split(filenames[i], "'time(15JUN24:03:45:36.0,15JUN24:04:20:00.0)'", ".timesplit", newfiles)

def create_uvf(filename):
    os.system("rm -r {}.uvf".format(filename))
    os.system('fits in={}.vis op=uvout out={}.uvf'.format(filename, filename))

#Make this image residuals? Doesn't seem that important
# def model_residuals(modelname, num, filename):
#         os.system('uvmodel model='+modelname+''+str(num)+'.im vis='+filename+'.vis options=subtract out='+modelname+''+str(num)+'_residuals.vis')
#         os.system('invert vis='+modelname+''+str(num)+'_residuals.vis map='+modelname+''+str(num)+'_residuals.mp beam='+modelname+''+str(num)+'_residuals.bm cell='+pixsize+' imsize=512 options=systemp,mfs robust=2')
#         os.system('clean map='+modelname+''+str(num)+'_residuals.mp beam='+modelname+''+str(num)+'_residuals.bm out='+modelname+''+str(num)+'_residuals.cl  niters=10000') #cutoff=6.290E-05
#         os.system('restor map='+modelname+''+str(num)+'_residuals.mp beam='+modelname+''+str(num)+'_residuals.bm model='+modelname+''+str(num)+'_residuals.cl out='+modelname+''+str(num)+'_residuals.cm')
#         os.system('invert vis='+modelname+''+str(num)+'.vis map='+modelname+''+str(num)+'.mp beam='+modelname+''+str(num)+'.bm cell='+pixsize+' imsize=512 options=systemp,mfs robust=2')
#         os.system('clean map='+modelname+''+str(num)+'.mp beam='+modelname+''+str(num)+'.bm out='+modelname+''+str(num)+'.cl  niters=10000') #cutoff=6.290E-05
#         os.system('restor map='+modelname+''+str(num)+'.mp beam='+modelname+''+str(num)+'.bm model='+modelname+''+str(num)+'.cl out='+modelname+''+str(num)+'.cm')
#         #os.system('cgdisp in='+modelname+''+str(num)+'.cm device=/xs labtyp=arcsec beamtyp=b,l,3')



recreate_uvf = True
if recreate_uvf:
    for filename in newfiles:
        create_uvf(filename)

filenames = ['18aug2015_aumic_spw0.corrected_weights','18aug2015_aumic_spw1.corrected_weights','18aug2015_aumic_spw2.corrected_weights','18aug2015_aumic_spw3.corrected_weights','24jun2015_aumic1_spw0.corrected_weights.timesplit','24jun2015_aumic1_spw1.corrected_weights.timesplit','24jun2015_aumic1_spw2.corrected_weights.timesplit','24jun2015_aumic1_spw3.corrected_weights.timesplit','26mar2014_aumic_spw0.corrected_weights','26mar2014_aumic_spw1.corrected_weights','26mar2014_aumic_spw2.corrected_weights','26mar2014_aumic_spw3.corrected_weights']
im = create_model()

chis = []
redchis = []
for i in range(len(filenames)):
    model_convolve(im, "aumic_model", i, filenames[i], coord[i], ra[i], dec[i])
    chi, redchi = get_chi(filenames[i], "aumic_model", i)

    chis.append(chi)
    redchis.append(redchi)

i=7
model_convolve(im, "aumic_model", i, filenames[i], coord[i], ra[i], dec[i])
chi, redchi = get_chi(filenames[i], "aumic_model", i)
