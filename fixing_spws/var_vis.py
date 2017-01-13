from astropy.io import fits
import numpy as np
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import subprocess


def var_vis(filename):
    ''' Calculate the variance in a visibility map at each u,v point and each channel. Calculate the variance based on the 70 nearest neighbors'''
    plt.clf()
    plt.close()

    im = fits.open('{}.uvf'.format(filename))
    u, v = im[0].data['UU'], im[0].data['VV']
    freq0 = im[0].header['crval4']
    klam = freq0 / 1e3
    # for now assume it is ALMA data
    vis = (im[0].data['data']).squeeze()
    if vis.shape[1] == 2:
        real = (vis[:, 0, 0] + vis[:, 1, 0]) / 2.
        imag = (vis[:, 0, 1] + vis[:, 1, 1]) / 2.
    else:
        real = vis[:, :, 0]
        imag = vis[:, :, 1]

    nuv = u.size
    nfreq = 1
    uv = u**2 + v**2
    nclose = 50  # number of nearby visibility points to use when measuring the dispersion
    # area around a particular uv point to consider when searching for the
    # nearest nclose neighbors (smaller numbers help make the good run faster,
    # but could result in many points for which the weight cannot be
    # calculated and is left at 0)
    uvwidth = 100
    max_dist = np.zeros(nuv)

    import time
    start = time.time()
    real_weight = np.zeros((nuv, nfreq), dtype=np.float32)
    real_nclose_arr = np.zeros(len(u))
    imag_weight = np.zeros((nuv, nfreq), dtype=np.float32)
    imag_nclose_arr = np.zeros(len(u))

    for iuv in range(nuv):
        w = (np.abs(u - u[iuv]) * klam <
             uvwidth) & (np.abs(v - v[iuv]) * klam < uvwidth)
        s = np.argsort(np.sqrt((v[w] - v[iuv])**2 + (u[w] - u[iuv])**2))

        # Reals
        real_wf = (real[w][s] != 0)
        real_nclose_arr[iuv] = real_wf.sum()
        if real_wf.sum() > nclose:
            real_weight[iuv] = 1 / np.std(real[w][s][real_wf][:nclose])**2
        else:
            print iuv, real_wf.sum(), np.sqrt(u[iuv]**2 + v[iuv]**2) * klam

        # Imaginaries
        imag_wf = (imag[w][s] != 0)
        imag_nclose_arr[iuv] = imag_wf.sum()
        if imag_wf.sum() > nclose:
            imag_weight[iuv] = 1 / np.std(imag[w][s][imag_wf][:nclose])**2
        else:
            print iuv, imag_wf.sum(), np.sqrt(u[iuv]**2 + v[iuv]**2) * klam

    print 'Elapsed time (hrs): ', (time.time() - start) / 3600.

    # plt.plot(np.sqrt(u**2+v**2)*klam,nclose_arr,'.')
    plt.scatter(real_weight, imag_weight)
    plt.plot(imag_weight, imag_weight)
    plt.xlabel("Real Weight")
    plt.ylabel("Imaginary Weight")
    plt.savefig("{}_weight_scatter.png".format(filename))

    return real_weight, imag_weight

    # diff = real_weight-imag_weight
    # plt.hist(diff, 20)
    # plt.xlabel("Real - Imaginary")
    # plt.ylabel("Number of Visibilities")
    # plt.savefig("{}_weight_difference_hist.png".format(filename))
    #
    #
    # Calculating total weight from imaginary and real components
    total_weight = np.sqrt(real_weight * imag_weight)

    # Reshape weights to fit in weights column of uv fits file; place the
    # weights into both xx and yy polarization columns.
    total_weight = np.reshape(total_weight, len(
        im[0].data['data'][:, 0, 0, 0, 0, 0, 2]))
    im[0].data['data'][:, 0, 0, 0, 0, 0, 2] = total_weight
    im[0].data['data'][:, 0, 0, 0, 0, 1, 2] = total_weight

    subprocess.call('rm {}.corrected_weights.uvf'.format(filename))
    im.writeto('{}.corrected_weights.uvf'.format(filename), clobber=True)
    im.close()
    return


# hdu=fits.PrimaryHDU(weight)
# hdu.writeto('HIP73145_weights2.fits')

uvf = ['18aug2015_aumic_spw0', '18aug2015_aumic_spw1', '18aug2015_aumic_spw2', '18aug2015_aumic_spw3',
       '24jun2015_aumic1_spw0', '24jun2015_aumic1_spw1', '24jun2015_aumic1_spw2', '24jun2015_aumic1_spw3', '26mar2014_aumic_spw0', '26mar2014_aumic_spw1', '26mar2014_aumic_spw2', '26mar2014_aumic_spw3']

uvfs = ['24jun2015_aumic1_spw0.timesplit', '24jun2015_aumic1_spw1.timesplit',
        '24jun2015_aumic1_spw2.timesplit', '24jun2015_aumic1_spw3.timesplit']
for uvf in uvfs:
    var_vis(uvf)

# real_weight, imag_weight = var_vis(uvf[0])
# real_weight0, imag_weight0 = var_vis(uvf[4])
# real_weight1, imag_weight1 = var_vis(uvf[5])
# real_weight2, imag_weight2 = var_vis(uvf[6])
# real_weight3, imag_weight3 = var_vis(uvf[7])
# real_weight4, imag_weight4 = var_vis(uvf[10])
