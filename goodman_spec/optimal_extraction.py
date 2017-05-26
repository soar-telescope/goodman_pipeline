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
from numpy import ma
from scipy import signal


def identify_targets(ccd, plots=False):
    """Identify targets cross correlating spatial profile with a gaussian model

    The method of cross-correlating a gaussian model to the spatial profile was
    mentioned in Marsh 1989, then I created my own implementation. The spatial
    profile is obtained by finding the median across the full dispersion axis.
    For goodman the spectrum and ccd are very well aligned, there is a slight
    variation from one configuration to another but in general is acceptable.

    Args:
        ccd (object): a ccdproc.CCDData instance

    Returns:
        profile_model (object): an astropy.modeling.Model instance, it could be
        a Gaussian1D or CompoundModel (several Gaussian1D). Each of them
        represent a point source spectrum found.

    """
    if isinstance(ccd, CCDData):
        slit_size = re.sub('[a-zA-Z"]', '', ccd.header['SLIT'])
        serial_binning = int(ccd.header['CCDSUM'].split()[0])
        # order will be used for finding the peaks later but also as an initial
        # estimate for stddev of gaussian
        order = int(round(float(slit_size) / (0.15 * serial_binning)))
        # averaging overall spectral dimension because in goodman the spectra is
        # deviated very little
        profile_median = np.median(ccd.data, axis=1)

        # Gaussian has to be defined at the middle for it to work
        gaussian = models.Gaussian1D(amplitude=profile_median.max(),
                                     mean=profile_median.size // 2,
                                     stddev=order).rename('Gaussian')

        # do cross correlation of data with gaussian
        # this is useful for cases with low signal to noise ratio
        cross_corr = signal.correlate(in1=profile_median,
                                      in2=gaussian(range(profile_median.size)),
                                      mode='same')
        # filter cross correlation data
        filt_cross_corr = np.where(np.abs(cross_corr > cross_corr.min()
                                          + 0.03 * cross_corr.max()),
                                   cross_corr,
                                   None)
        peaks = signal.argrelmax(filt_cross_corr, axis=0, order=order)[0]

        profile_model = None
        for i in range(len(peaks)):
            low_lim = np.max([0, peaks[i] - 5])
            hig_lim = np.min([peaks[i] + 5, profile_median.size])
            amplitude = profile_median[low_lim: hig_lim].max()
            gaussian = models.Gaussian1D(amplitude=amplitude,
                                         mean=peaks[i],
                                         stddev=order).rename('Gaussian_{:d}'.format(i))
            if profile_model is not None:
                profile_model += gaussian
            else:
                profile_model = gaussian
        #     plt.axvline(peaks[i])
        # print(profile_model)

        # fit model to profile
        fit_gaussian = fitting.LevMarLSQFitter()
        profile_model = fit_gaussian(model=profile_model,
                                     x=range(profile_median.size),
                                     y=profile_median)
        # print(fit_gaussian.fit_info)
        if plots:
            # plt.plot(cross_corr, label='Cross Correlation', linestyle='--')
            plt.plot(profile_model(range(profile_median.size)),
                     label='Fitted Model')
            plt.plot(profile_median, label='Profile (median)', linestyle='--')
            plt.legend(loc='best')
            plt.show()

        # print(profile_model)
        if fit_gaussian.fit_info['ierr'] not in [1, 2, 3, 4]:
            log.warning('There is some problem with the fitting. Returning None.')
            return None
        else:
            return profile_model

    else:
        log.error('Not a CCDData instance')
        return None


def trace_targets(ccd, profile, sampling_step=5, pol_deg=2, plots=True):
    """Find the trace of the target's spectrum on the image

    Args:
        ccd (object): Instance of ccdproc.CCDData
        profile (object): Instance of astropy.modeling.Model, contains the
        spatial profile of the 2D spectrum.
        sampling_step (int): Frequency of sampling in pixels
        pol_deg (int): Polynomial degree for fitting the trace
        plots (bool): If True will show plots (debugging)

    Returns:

    """
    # added two assert for debugging purposes
    assert isinstance(ccd, CCDData)
    assert isinstance(profile, Model)

    # Get image dimensions
    spatial_length, dispersion_length = ccd.data.shape

    # Initialize model fitter
    model_fitter = fitting.LevMarLSQFitter()

    # Initialize the model to fit the traces
    trace_model = models.Polynomial1D(degree=pol_deg)

    # Will store the arrays for the fitted location of the target obtained
    # in the fitting
    trace_points = None

    # List that will contain all the Model instances corresponding to traced
    # targets
    all_traces = None

    # Array defining the points to be sampled
    sampling_axis = range(0,
                          dispersion_length // sampling_step * sampling_step,
                          sampling_step)

    # Loop to go through all the sampling points
    for i in sampling_axis:
        # Fit the inital model to the data
        fitted_profile = model_fitter(model=profile,
                                      x=range(ccd.data[:, i].size),
                                      y=ccd.data[:, i])
        if model_fitter.fit_info['ierr'] not in [1, 2, 3, 4]:
            log.error(
                "Fitting did not work fit_info['ierr'] = \
                {:d}".format(model_fitter.fit_info['ierr']))

        # alternatively could use fitted_profile.param_names
        # the mean_keys is another way to tell how many submodels are in the
        # model that was parsed.
        mean_keys = [key for key in dir(fitted_profile) if 'mean' in key]

        if trace_points is None:
            trace_points = np.ndarray((len(mean_keys),
                                       dispersion_length // sampling_step))

        # store the corresponding value in the proper array for later fitting
        # a low order polinomial
        for e in range(trace_points.shape[0]):
            trace_points[e][i // sampling_step] =\
               fitted_profile.__getattribute__(mean_keys[e]).value

    # fit a low order polynomial for the trace
    for trace_num in range(trace_points.shape[0]):
        fitted_trace = model_fitter(model=trace_model,
                                    x=sampling_axis,
                                    y=trace_points[trace_num])

        fitted_trace.rename('Trace_{:d}'.format(trace_num))

        if model_fitter.fit_info['ierr'] not in [1, 2, 3, 4]:
            log.error(model_fitter.fit_info['ierr'])
        else:
            # RMS Error calculation
            RMSError = np.sqrt(
                np.sum(np.array([fitted_trace(sampling_axis) -
                                 trace_points[trace_num]]) ** 2))
            log.info('Trace Fit RMSE: {:.3f}'.format(RMSError))

            if all_traces is None:
                all_traces = [fitted_trace]
            else:
                all_traces.append(fitted_trace)

        if plots:
            plt.plot(sampling_axis,
                     trace_points[trace_num],
                     marker='o',
                     label='Data Points')
            plt.plot(fitted_trace(range(dispersion_length)),
                     label='Model RMSE: {:.2f}'.format(RMSError))

    if plots:
        plt.title('Targets Trace')
        plt.xlabel('Dispersion Axis')
        plt.ylabel('Spatial Axis')
        plt.imshow(ccd.data, cmap='YlGnBu', clim=(0, 300))

        for trace in trace_points:
            plt.plot(sampling_axis, trace, marker='.', linestyle='--')
            # print(trace)
        plt.legend(loc='best')
        plt.show()
    return all_traces



def remove_background(ccd, profile_model=None, plots=False):
    ccd_copia = ccd.copy()
    data = ma.masked_invalid(ccd.data)
    x, y = ccd.data.shape
    median = ma.median(data, axis=0)

    data -= median
    data.set_fill_value(-np.inf)
    ccd.data = data.filled()

    # ccd.write('/user/simon/dummy_{:d}.fits'.format(g), clobber=True)
    return ccd


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


def get_extraction_zone(ccd, model, n_sigma_extract, plots=False):

    spatial_length, dispersion_length = ccd.data.shape

    mean = model.mean.value
    stddev = model.stddev.value
    extract_width = n_sigma_extract // 2 * stddev

    low_lim = np.max([0, int(mean - extract_width)])
    hig_lim = np.min([int(mean + extract_width), spatial_length])

    nccd = ccd.copy()
    nccd.data = ccd.data[low_lim:hig_lim, :]
    nccd.header['HISTORY'] = 'Subsection of CCD [{:d}:{:d}, :]'.format(low_lim,
                                                                       hig_lim)
    model.mean.value = extract_width

    if plots:
        plt.imshow(ccd.data)
        plt.axhspan(low_lim, hig_lim, color='r', alpha=0.2)
        plt.show()

    return nccd, model


def optimal_extraction(ccd, n_sigma_extract=10, plots=False):

    assert isinstance(ccd, CCDData)

    ccd = remove_background(ccd=ccd)
    profile_model = identify_targets(ccd=ccd)
    if isinstance(profile_model, Model):
        trace_targets(ccd=ccd, profile=profile_model)
        # extract(ccd=ccd,
        #         spatial_profile=profile_model,
        #         n_sigma_extract=10,
        #         sampling_step=5)
        if 'CompoundModel' in profile_model.__class__.name:
            for submodel_name in profile_model.submodel_names:
                nccd, model = get_extraction_zone(ccd=ccd,
                                                  model=profile_model[submodel_name],
                                                  n_sigma_extract=n_sigma_extract)
                extract(ccd=nccd, spatial_profile=model, sampling_step=10)
                if plots:
                    plt.imshow(nccd)
                    plt.show()

        else:
            nccd, model = get_extraction_zone(ccd=ccd,
                                              model=profile_model,
                                              n_sigma_extract=n_sigma_extract)
            extract(ccd=nccd, spatial_profile=model, sampling_step=10)
            if plots:
                plt.imshow(nccd)
                plt.show()

    elif profile_model is None:
        log.warning("Didn't receive identified targets.")
    else:
        log.error('Got wrong input')

if __name__ == '__main__':
    files = ['/data/simon/data/soar/work/GOO_BLU_SPE_2017-05-09/RED/cfzsto_1086_HD104237_400M2_GG455.fits',
             '/data/simon/data/soar/work/20161114_eng_3/RED4/cfzsto_0216_EG21_1200M5_GG455.fits',
             '/user/simon/Downloads/goodman_target_simulator-master/dummy2.fits',
             '/data/simon/data/soar/work/GOO_BLU_SPE_2017-05-09/RED/cfzsto_1119_HgAr_400M2_GG455.fits']
    # files = glob.glob('/data/simon/data/soar/work/GOO_BLU_SPE_2017-05-09/RED/cfzsto*')
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

