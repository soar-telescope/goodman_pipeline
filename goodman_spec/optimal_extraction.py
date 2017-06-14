from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os

import matplotlib
# matplotlib.use('Qt4Agg')
matplotlib.use('GTK3Agg')

import astropy.units as u
import logging
import numpy as np
import matplotlib.pyplot as plt
import re
import glob
import pandas

from astropy.modeling import (models, fitting, Model)
from ccdproc import CCDData
from ccdproc import ImageFileCollection
from goodman_ccd.core import NightDataContainer
from goodman_ccd.core import cosmicray_rejection
from goodman_ccd.core import identify_targets
from goodman_ccd.core import trace_targets
from goodman_ccd.core import get_extraction_zone
from goodman_ccd.core import remove_background
from goodman_ccd.core import ra_dec_to_deg
from numpy import ma
from scipy import signal

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('goodmanccd')

# log.setLevel(level=logging.DEBUG)


def extract(ccd, trace, spatial_profile, sampling_step=1, plots=False):

    assert isinstance(ccd, CCDData)
    assert isinstance(trace, Model)

    spatial_length, dispersion_length = ccd.data.shape

    # create variance model
    rdnoise = float(ccd.header['RDNOISE'])
    gain = float(ccd.header['GAIN'])
    print('Original Name {:s}'.format(ccd.header['OFNAME']))

    variance_2d = (rdnoise + np.absolute(ccd.data) * gain) / gain
    if ccd.header['OBSTYPE'] == 'OBJECT':
        cr_mask = cosmicray_rejection(ccd=ccd, mask_only=True)
        cr_mask = np.logical_not(cr_mask).astype(int)
    else:
        log.info('Only OBSTYPE == OBJECT get cosmic ray rejection.')
        cr_mask = np.ones(ccd.data.shape, dtype=int)

    model_fitter = fitting.LevMarLSQFitter()

    # print(spatial_profile.mean.value)
    # print(trace.c0.value)

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
    log.debug('ccd.data is a masked array: '
              '{:s}'.format(str(np.ma.isMaskedArray(ccd.data))))

    ccd.data = np.ma.masked_invalid(ccd.data)
    # print(np.ma.isMaskedArray(ccd.data))
    np.ma.set_fill_value(ccd.data, 0)

    # print(np.ma.isMaskedArray(ccd.data))
    to_sum = ccd.data * cr_mask
    simple_sum = np.ma.sum(to_sum, axis=0)

    out_spectrum = np.empty(dispersion_length)
    for i in range(0, dispersion_length, sampling_step):
        # force the model to follow the trace
        new_model.mean.value = trace(i)

        # warn if the difference of the spectrum position in the trace at the
        # extremes of the sampling range is larger than 1 pixel.
        if np.abs(trace(i) - trace(i+sampling_step)) > 1:
            log.warning('Sampling step might be too large')

        sample = np.median(ccd.data[:, i:i + sampling_step], axis=1)
        fitted_profile = model_fitter(model=new_model,
                                      x=range(len(sample)),
                                      y=sample)

        profile = fitted_profile(range(sample.size))

        # enforce positivity
        pos_profile = np.array([np.max([0, x]) for x in profile])

        # enforce normalization
        nor_profile = np.array([x / pos_profile.sum() for x in pos_profile])

        if sampling_step > 1:
            # TODO (simon): Simplify to Pythonic way
            right = min((i + sampling_step), dispersion_length)

            for e in range(i, right, 1):
                mask = cr_mask[:,e]
                data = ma.masked_invalid(ccd.data[:, e])
                # print(ma.isMaskedArray(data))
                V = variance_2d[:, e]
                P = nor_profile
                a = [(P[z] / V[z]) / np.sum(P ** 2 / V) for z in range(P.size)]
                weights = (nor_profile / variance_2d[:, e]) / np.sum(nor_profile ** 2 / variance_2d[:, e])
                # print('SUMN ', np.sum(a), np.sum(weights), np.sum(nor_profile), np.sum(P * weights))



                # if e in range(5, 4001, 500):
                #     plt.plot(nor_profile * data.max()/ nor_profile.max(), label=str(e))
                #     plt.plot(data, label='Data')
                #     plt.legend(loc='best')
                #     plt.show()

                out_spectrum[e] = np.sum(data * mask * nor_profile)


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

    ccd.data = simple_sum
    # ccd.data = out_spectrum
    return ccd


def optimal_extraction(ccd, comp_list=None, n_sigma_extract=10, plots=False):

    assert isinstance(ccd, CCDData)

    comp_zones = []
    extracted = []

    if comp_list is None:
        comp_list = []

    iccd = remove_background(ccd=ccd,)

    profile_model = identify_targets(ccd=iccd, plots=True)
    del(iccd)

    if isinstance(profile_model, Model):
        traces = trace_targets(ccd=ccd, profile=profile_model, plots=True)
        # extract(ccd=ccd,
        #         spatial_profile=profile_model,
        #         n_sigma_extract=10,
        #         sampling_step=5)
        if 'CompoundModel' in profile_model.__class__.name:
            log.debug(profile_model.submodel_names)
            for m in range(len(profile_model.submodel_names)):

                submodel_name = profile_model.submodel_names[m]

                nccd, trace, model, zone = get_extraction_zone(
                    ccd=ccd,
                    trace=traces[m],
                    model=profile_model[submodel_name],
                    n_sigma_extract=n_sigma_extract,
                    plots=True)

                for comp in comp_list:
                    comp_zone = get_extraction_zone(ccd=comp,
                                                    trace=trace,
                                                    zone=zone)
                    # since a comparison lamp only needs only the relative line
                    # center in the dispersion direction, therefore the flux is
                    # not important we are only calculating the median along the
                    # spatial direction
                    comp_zone.data = np.median(comp_zone.data, axis=0)
                    comp_zones.append(comp_zone)

                original = nccd.copy()
                original.write('/data/simon/original.fits', clobber=True)

                nccd = remove_background(ccd=nccd,
                                         plots=True)

                nccd.write('/data/simon/bkg-removed.fits', clobber=True)

                original.data -= nccd.data

                original.write('/data/simon/background.fits', clobber=True)
                print('Written difference')

                extracted_ccd = extract(ccd=nccd,
                                        trace=trace,
                                        spatial_profile=model,
                                        sampling_step=10,
                                        plots=True)
                extracted.append(extracted_ccd)

                if plots:
                    plt.imshow(nccd)
                    plt.show()

        else:
            nccd, trace, model, zone = get_extraction_zone(
                ccd=ccd,
                trace=traces[0],
                model=profile_model,
                n_sigma_extract=n_sigma_extract,
                plots=True)

            for comp in comp_list:
                comp_zone = get_extraction_zone(ccd=comp,
                                                trace=trace,
                                                zone=zone)

                # since a comparison lamp only needs only the relative line
                # center in the dispersion direction, therefore the flux is not
                # important we are only calculating the median along the spatial
                # direction
                comp_zone.data = np.median(comp_zone.data, axis=0)
                comp_zones.append(comp_zone)

            nccd = remove_background(ccd=nccd,
                                     plots=True)

            extracted_ccd = extract(ccd=nccd,
                                    trace=trace,
                                    spatial_profile=model,
                                    sampling_step=10,
                                    plots=True)

            extracted.append(extracted_ccd)

            if plots:
                plt.imshow(nccd)
                plt.show()

        print(extracted)
        print(comp_zones)
        return extracted, comp_zones

    elif profile_model is None:
        log.warning("Didn't receive identified targets.")
    else:
        log.error('Got wrong input')


def process_spectroscopy_data(data_container):

    assert data_container.is_empty is False

    full_path = data_container.full_path

    for spec_group in data_container.spec_groups:
        comp_group = None
        object_group = None
        comp_ccd_list = []
        obstypes = spec_group.obstype.unique()
        if 'COMP' in obstypes:
            comp_group = spec_group[spec_group.obstype == 'COMP']
            # print(comp_group)
        if 'OBJECT' in obstypes:
            object_group = spec_group[spec_group.obstype == 'OBJECT']
            # print(object_group)
        for spec_file in object_group.file.tolist():
            file_path = os.path.join(full_path, spec_file)
            ccd = CCDData.read(file_path, unit=u.adu)
            ccd.header['OFNAME'] = (spec_file, 'Original File Name')
            if comp_group is not None:
                for comp_file in comp_group.file.tolist():
                    comp_path = os.path.join(full_path, comp_file)
                    comp_ccd = CCDData.read(comp_path, unit=u.adu)
                    comp_ccd.header['OFNAME'] = (comp_file,
                                                 'Original File Name')
                    comp_ccd_list.append(comp_ccd)

            extracted, comps = optimal_extraction(ccd=ccd, comp_list=comp_ccd_list)
            if True:
                manager = plt.get_current_fig_manager()
                manager.window.maximize()
                for edata in extracted:
                    plt.plot(edata.data, label=edata.header['OBJECT'])
                    if comps != []:
                        for comp in comps:
                            plt.plot(comp.data, label=comp.header['OBJECT'])
                plt.legend(loc='best')
                plt.show()

    print('\nEND')

if __name__ == '__main__':
    # # files = ['/data/simon/data/soar/work/GOO_BLU_SPE_2017-05-09/RED/cfzsto_1086_HD104237_400M2_GG455.fits',
    # #          '/data/simon/data/soar/work/20161114_eng_3/RED4/cfzsto_0216_EG21_1200M5_GG455.fits',
    # #          '/user/simon/Downloads/goodman_target_simulator-master/dummy2.fits',
    # #          '/data/simon/data/soar/work/GOO_BLU_SPE_2017-05-09/RED/cfzsto_1119_HgAr_400M2_GG455.fits']
    # # files = glob.glob('/data/simon/data/soar/work/GOO_BLU_SPE_2017-05-09/RED/cfzsto*')
    # files = ['/user/simon/data/soar/observer/SciCase001/2014-11-10/SYZY_400M2/RED/cfzsto_0074.cvso114_400M2_GG455.fits',
    #          '/user/simon/data/soar/observer/SciCase001/2014-11-10/SYZY_400M2/RED/cfzsto_0075.cvso114_400M2_GG455.fits',
    #          '/user/simon/data/soar/observer/SciCase001/2014-11-10/SYZY_400M2/RED/cfzsto_0076.cvso114_400M2_GG455.fits']

    pandas.set_option('display.expand_frame_repr', False)

    prefix = 'cfzsto'
    path = '/user/simon/data/soar/work/20161114_eng_3/RED'
    search_path = os.path.join(path, prefix + '*.fits')

    data_container = NightDataContainer(path=path,
                                        instrument='Red',
                                        technique='Spectroscopy')

    file_list = glob.glob(search_path)

    # file_list = ['cfzsto_0276_HILT600_1200M5_GG455.fits',
    #              'cfzsto_0278_HILT600_1200M5_GG455.fits',
    #              'cfzsto_0279_HILT600_2100.fits']

    keywords = ['date',
                'slit',
                'date-obs',
                'obstype',
                'object',
                'exptime',
                'obsra',
                'obsdec',
                'grating',
                'cam_targ',
                'grt_targ',
                'filter',
                'filter2',
                'gain',
                'rdnoise']

    ifc = ImageFileCollection(path, keywords=keywords, filenames=file_list)

    pifc = ifc.summary.to_pandas()

    pifc['radeg'] = ''
    pifc['decdeg'] = ''
    for i in pifc.index.tolist():
        radeg, decdeg = ra_dec_to_deg(pifc.obsra.iloc[i], pifc.obsdec.iloc[i])

        pifc.iloc[i, pifc.columns.get_loc('radeg')] = '{:.2f}'.format(radeg)

        pifc.iloc[i, pifc.columns.get_loc('decdeg')] = '{:.2f}'.format(decdeg)
        # now we can compare using degrees

    confs = pifc.groupby(['slit',
                          'radeg',
                          'decdeg',
                          'grating',
                          'cam_targ',
                          'grt_targ',
                          'filter',
                          'filter2',
                          'gain',
                          'rdnoise']).size().reset_index().rename(columns={0: 'count'})

    for i in confs.index:
        spec_group = pifc[((pifc['slit'] == confs.iloc[i]['slit']) &
                           (pifc['radeg'] == confs.iloc[i]['radeg']) &
                           (pifc['decdeg'] == confs.iloc[i]['decdeg']) &
                           (pifc['grating'] == confs.iloc[i]['grating']) &
                           (pifc['cam_targ'] == confs.iloc[i]['cam_targ']) &
                           (pifc['grt_targ'] == confs.iloc[i]['grt_targ']) &
                           (pifc['filter'] == confs.iloc[i]['filter']) &
                           (pifc['filter2'] == confs.iloc[i]['filter2']) &
                           (pifc['gain'] == confs.iloc[i]['gain']) &
                           (pifc['rdnoise'] == confs.iloc[i]['rdnoise']))]
        data_container.add_spec_group(spec_group=spec_group)

    # print(confs)


    # print(pifc)
    process_spectroscopy_data(data_container=data_container)


    # i=10
    # for sfile in files:
    #     ccd = CCDData.read(sfile, unit=u.adu)
    #
    #     if ccd.data.ndim == 3:
    #         ccd.data = ccd.data[0]
    #     # for z in range(ccd.data.shape[1]):
    #     #     plt.plot(ccd.data[:,z])
    #     # plt.plot()
    #
    #     optimal_extraction(ccd=ccd)
    #     # remove_background(ccd=ccd, g=i)
    #     i += 1

