from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import matplotlib
# matplotlib.use('Qt4Agg')
matplotlib.use('GTK3Agg')

import astropy.units as u
import logging as log
import numpy as np
import matplotlib.pyplot as plt
import re
import glob

from astropy.modeling import (models, fitting, Model)
from ccdproc import CCDData
from goodman_ccd.core import cosmicray_rejection
from goodman_ccd.core import identify_targets
from goodman_ccd.core import trace_targets
from goodman_ccd.core import get_extraction_zone
from goodman_ccd.core import remove_background
from numpy import ma
from scipy import signal





def extract(ccd, spatial_profile, sampling_step=1, plots=False):

    assert isinstance(ccd, CCDData)

    spatial_length, dispersion_length = ccd.data.shape

    # create variance model
    rdnoise = float(ccd.header['RDNOISE'])
    gain = float(ccd.header['GAIN'])

    variance2D = (rdnoise + np.absolute(ccd.data) * gain) / gain
    if ccd.header['OBSTYPE'] == 'OBJECT':
        cr_mask = cosmicray_rejection(ccd=ccd, mask_only=True)
        cr_mask = np.logical_not(cr_mask).astype(int)
    else:
        log.info('only objects get cosmic ray rej')
        cr_mask = np.ones(ccd.data.shape, dtype=int)
    # print(cr_mask)
    # if False in cr_mask:
    #     print('YEs! there is one')
    # plt.imshow(cr_mask)
    # plt.show()

    model_fitter = fitting.LevMarLSQFitter()

    if isinstance(spatial_profile, models.Gaussian1D):
        amplitude = spatial_profile.amplitude.value
        mean = spatial_profile.mean.value
        stddev = spatial_profile.stddev.value
        new_model = models.Gaussian1D(amplitude=amplitude,
                                      mean=mean,
                                      stddev=stddev)
        # print('Fixed ', new_model.mean.fixed)
    else:
        raise NotImplementedError
    simple_sum = np.sum(ccd.data, axis=0)

    # plt.show()
    out_spectrum = np.empty(dispersion_length)
    for i in range(0, dispersion_length, sampling_step):
        sample = np.median(ccd.data[:, i:i + sampling_step], axis=1)
        fitted_profile = model_fitter(model=new_model,
                                      x=range(len(sample)),
                                      y=sample)
        profile = fitted_profile(range(sample.size))
        # enforce positivity
        pos_profile = np.array([np.max([0, x]) for x in profile])
        # enforce normalization
        nor_profile = np.array([x / pos_profile.sum() for x in pos_profile])

        # print(i)

        if sampling_step > 1:
            # TODO (simon): Simplify to Pythonic way
            right = min((i + sampling_step), dispersion_length)

            for e in range(i, right, 1):
                mask = cr_mask[:,e]
                data  = ma.masked_invalid(ccd.data[:, e])
                V = variance2D[:, e]
                P = nor_profile
                a = [(P[z] / V[z]) / np.sum(P ** 2 / V) for z in range(P.size)]
                weights = (nor_profile / variance2D[:, e]) / np.sum(nor_profile ** 2 / variance2D[:, e])
                # print('SUMN ', np.sum(a), np.sum(weights), np.sum(nor_profile), np.sum(P * weights))
                # if 1500 < e < 1510:
                #     plt.ion()
                #     plt.title('Variance')
                #     plt.plot(variance2D[:, e], label='Variance')
                #     # plt.show()
                #     plt.title('Data')
                #     plt.plot(data, label='Data')
                #     plt.legend(loc='best')
                #     # plt.show()
                #     plt.draw()
                #     plt.pause(1)
                #     plt.clf()
                #     plt.title('Norm Profile')
                #     plt.plot(nor_profile, label='Normalized Profile')
                #     # plt.show()
                #     plt.draw()
                #     plt.pause(1)
                #     plt.clf()
                #     plt.title('Weights')
                #     plt.plot(weights, label='Old')
                #     plt.plot(a, label='New')
                #     plt.legend(loc='best')
                #     # plt.show()
                #     plt.draw()
                #     plt.pause(1)
                #     plt.clf()
                #     plt.ioff()

                out_spectrum[e] = np.sum(data * mask * nor_profile)
                # print(np.sum(ccd.data[:, e] * mask * nor_profile))
                # if 0 in cr_mask[:,e]:
                #     print(cr_mask[:,e])
                # # print(e)

        # print(cr_mask[:,i])
        # plt.plot(cr_mask[:,i], color='k')
        # plt.plot(nor_profile)
        # plt.show()
    if plots:
        manager = plt.get_current_fig_manager()
        manager.window.maximize()
        ratio = np.mean(simple_sum / out_spectrum)
        # print(ratio)
        plt.title(ccd.header['OBJECT'])
        plt.plot(simple_sum, label='Simple Sum', color='k', alpha=0.5)
        plt.plot(out_spectrum * ratio, color='r',
                 label='Optimal * {:.1f}'.format(ratio), alpha=0.5)
        plt.legend(loc='best')
        plt.show()

        #     plt.ion()
        #     plt.title(i)
        #     plt.plot(fitted_profile(range(len(sample))), color='r', alpha=0.1)
        #     plt.plot(sample, color='k', alpha=0.1)
        #     plt.draw()
        #     plt.pause(.1)
        # plt.ioff()
        # plt.close()


def optimal_extraction(ccd, n_sigma_extract=10, plots=False):

    assert isinstance(ccd, CCDData)

    iccd = remove_background(ccd=ccd, plots=True)
    profile_model = identify_targets(ccd=iccd, plots=True)
    del(iccd)
    if isinstance(profile_model, Model):
        trace_targets(ccd=ccd, profile=profile_model, plots=True)
        # extract(ccd=ccd,
        #         spatial_profile=profile_model,
        #         n_sigma_extract=10,
        #         sampling_step=5)
        if 'CompoundModel' in profile_model.__class__.name:
            for submodel_name in profile_model.submodel_names:

                nccd, model, zone = get_extraction_zone(
                    ccd=ccd,
                    model=profile_model[submodel_name],
                    n_sigma_extract=n_sigma_extract,
                    plots=True)

                nccd = remove_background(ccd=nccd,
                                         plots=True)

                extract(ccd=nccd,
                        spatial_profile=model,
                        sampling_step=10,
                        plots=True)

                if plots:
                    plt.imshow(nccd)
                    plt.show()

        else:

            nccd, model, zone = get_extraction_zone(
                ccd=ccd,
                model=profile_model,
                n_sigma_extract=n_sigma_extract,
                plots=True)

            nccd = remove_background(ccd=nccd,
                                         plots=True)

            extract(ccd=nccd,
                    spatial_profile=model,
                    sampling_step=10,
                    plots=True)

            if plots:
                plt.imshow(nccd)
                plt.show()

    elif profile_model is None:
        log.warning("Didn't receive identified targets.")
    else:
        log.error('Got wrong input')

if __name__ == '__main__':
    # files = ['/data/simon/data/soar/work/GOO_BLU_SPE_2017-05-09/RED/cfzsto_1086_HD104237_400M2_GG455.fits',
    #          '/data/simon/data/soar/work/20161114_eng_3/RED4/cfzsto_0216_EG21_1200M5_GG455.fits',
    #          '/user/simon/Downloads/goodman_target_simulator-master/dummy2.fits',
    #          '/data/simon/data/soar/work/GOO_BLU_SPE_2017-05-09/RED/cfzsto_1119_HgAr_400M2_GG455.fits']
    # files = glob.glob('/data/simon/data/soar/work/GOO_BLU_SPE_2017-05-09/RED/cfzsto*')
    files = ['/user/simon/data/soar/observer/SciCase001/2014-11-10/SYZY_400M2/RED/cfzsto_0074.cvso114_400M2_GG455.fits',
             '/user/simon/data/soar/observer/SciCase001/2014-11-10/SYZY_400M2/RED/cfzsto_0075.cvso114_400M2_GG455.fits',
             '/user/simon/data/soar/observer/SciCase001/2014-11-10/SYZY_400M2/RED/cfzsto_0076.cvso114_400M2_GG455.fits']




    i=10
    for sfile in files:
        ccd = CCDData.read(sfile, unit=u.adu)

        if ccd.data.ndim == 3:
            ccd.data = ccd.data[0]
        # for z in range(ccd.data.shape[1]):
        #     plt.plot(ccd.data[:,z])
        # plt.plot()

        optimal_extraction(ccd=ccd)
        # remove_background(ccd=ccd, g=i)
        i += 1

