from astropy.io import fits
import numpy as np
from scipy.spatial.distance import pdist,squareform
import matplotlib.pyplot as plt
import subprocess
from glob import glob


def var_vis(infile, outfile):
    ''' Calculate the variance in a visibility map at each u,v point and each channel. Calculate the variance based on the 70 nearest neighbors'''
    plt.clf(); plt.close()

    im = fits.open('{}.uvf'.format(infile))
    u,v = im[0].data['UU'],im[0].data['VV']
    freq0 = im[0].header['crval4']
    klam = freq0/1e3
    # for now assume it is ALMA data
    vis = (im[0].data['data']).squeeze()
    if vis.shape[1] == 2:
        real = (vis[:,0,0]+vis[:,1,0])/2.
        imag = (vis[:,0,1]+vis[:,1,1])/2.
    else:
        real = vis[:,:,0]
        imag = vis[:,:,1]

    nuv = u.size
    print(infile)
    print("{} visibilities".format(nuv))
    print(" Max u baseline = {}".format(np.max(u)*klam))
    print(" Max v baseline = {}".format(np.max(v)*klam))
    nfreq = 1
    uv = u**2+v**2
    nclose = 1500 #number of nearby visibility points to use when measuring the dispersion
    uvwidth = 22 #area around a particular uv point to consider when searching for the nearest nclose neighbors (smaller numbers help make the good run faster, but could result in many points for which the weight cannot be calculated and is left at 0)
    max_dist = np.zeros(nuv)

    import time
    start=time.time()
    real_weight = np.zeros((nuv,nfreq), dtype=np.float32)
    real_nclose_arr = np.zeros(len(u))
    imag_weight = np.zeros((nuv,nfreq), dtype=np.float32)
    imag_nclose_arr = np.zeros(len(u))

    bad_points = [] #list of bad point indices used to calculate 0-weight point percentage
    for iuv in range(nuv):
        w = (np.abs(u-u[iuv])*klam < uvwidth) & (np.abs(v-v[iuv])*klam < uvwidth)
        s = np.argsort(np.sqrt((v[w]-v[iuv])**2+(u[w]-u[iuv])**2))

        #Reals
        real_wf = (real[w][s] !=0)
        real_nclose_arr[iuv] = real_wf.sum()
        if real_wf.sum()>nclose:
            real_weight[iuv] = 1/np.std(real[w][s][real_wf][:nclose])**2
        else:
            #print iuv,real_wf.sum(),np.sqrt(u[iuv]**2+v[iuv]**2)*klam
            bad_points.append(iuv)

        #Imaginaries
        imag_wf = (imag[w][s] !=0)
        imag_nclose_arr[iuv] = imag_wf.sum()
        if imag_wf.sum()>nclose:
            imag_weight[iuv] = 1/np.std(imag[w][s][imag_wf][:nclose])**2
        else:
            #print iuv,imag_wf.sum(),np.sqrt(u[iuv]**2+v[iuv]**2)*klam
            if bad_points[-1] != iuv:
                bad_points.append(iuv)

    # plt.plot(np.sqrt(u**2+v**2)*klam,nclose_arr,'.')
    #plt.scatter(real_weight, imag_weight)
    #plt.plot(imag_weight, imag_weight)
    #plt.xlabel("Real Weight")
    #plt.ylabel("Imaginary Weight")
    #plt.savefig("{}_weight_scatter.png".format(infile))

    #Calculating total weight from imaginary and real components
    total_weight = np.sqrt(real_weight*imag_weight)

    #Reshape weights to fit in weights column of uv fits file; place the weights into both xx and yy polarization columns.
    if len(im[0].data['data'].shape) == 7:
        total_weight = np.reshape(total_weight, len(im[0].data['data'][:,0,0,0,0,0,2]))
        im[0].data['data'][:,0,0,0,0,0,2], im[0].data['data'][:,0,0,0,0,1,2] = total_weight, total_weight
    else:
        total_weight = np.reshape(total_weight, len(im[0].data['data'][:,0,0,0,0,2]))
        im[0].data['data'][:,0,0,0,0,2], im[0].data['data'][:,0,0,0,1,2] = total_weight, total_weight


    subprocess.call('rm {}.uvf'.format(outfile), shell=True)
    im.writeto('{}.uvf'.format(outfile), overwrite=True)
    im.close()

    print("{}% visibilities have zero weight".format((len(bad_points)/float(nuv)) / 100.))
    print('Elapsed time (min): {}'.format((time.time()-start)/60.))

    return real_weight, imag_weight


def create_vis(filename):
    subprocess.call('rm -r {}.vis'.format(filename), shell=True)
    subprocess.call('fits in={}.uvf op=uvin out={}.vis'.format(
        filename, filename), shell=True)


# uvfs = ['18aug2015_aumic_spw0','18aug2015_aumic_spw1','18aug2015_aumic_spw2','18aug2015_aumic_spw3','24jun2015_aumic1_spw0','24jun2015_aumic1_spw1','24jun2015_aumic1_spw2','24jun2015_aumic1_spw3','26mar2014_aumic_spw0','26mar2014_aumic_spw1','26mar2014_aumic_spw2','26mar2014_aumic_spw3']

# files = glob("*.timing.uvf")
# files = [f[:-4] for f in files]
# print(files)
# input('ok?')
# start_dir = './'
# end_dir = '../../data_files/'
# 
# for f in files:
#     var_vis(f)
#     create_vis("{}.reweighted".format(f))

#var_vis("../24jun2015_aumic1_spw3.timeflag.corrected_weights")
#var_vis("24jun2015_aumic1_spw3")
#create_vis("24jun2015_aumic1_spw3.corrected_weights")
