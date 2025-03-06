from __future__ import absolute_import

from unittest import TestCase, skip
from ccdproc import CCDData
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from astropy.io import fits
from astropy.modeling import Model
from astropy.modeling import (models,
                              fitting)
import astropy.units as u
import collections
import mock
import numpy as np
import os
import pandas
import random
import re
import shutil
import logging

logging.disable(logging.CRITICAL)

# import all classes in core.py
from ..core import (GenerateDcrParFile,
                    NightDataContainer,
                    NoMatchFound,
                    NotEnoughLinesDetected,
                    NoTargetException,
                    ReferenceData,
                    SaturationValues,
                    SpectroscopicMode)


# import of functions in core.py
from ..core import (astroscrappy_lacosmic,
                    add_linear_wavelength_solution,
                    add_wcs_keys,
                    bias_subtract,
                    bin_reference_data,
                    call_cosmic_rejection,
                    classify_spectroscopic_data,
                    combine_data,
                    convert_time,
                    create_master_bias,
                    create_master_flats,
                    cross_correlation,
                    dcr_cosmicray_rejection,
                    define_trim_section,
                    extraction,
                    extract_fractional_pixel,
                    extract_optimal,
                    evaluate_wavelength_solution,
                    fix_keywords,
                    fractional_sum,
                    get_best_flat,
                    get_central_wavelength,
                    get_lines_in_lamp,
                    get_overscan_region,
                    get_spectral_characteristics,
                    get_slit_trim_section,
                    get_twilight_time,
                    identify_targets,
                    image_overscan,
                    image_trim,
                    interpolate_spectrum,
                    is_file_saturated,
                    linearize_spectrum,
                    name_master_flats,
                    normalize_master_flat,
                    ra_dec_to_deg,
                    read_fits,
                    record_trace_information,
                    save_extracted,
                    search_comp_group,
                    setup_logging,
                    trace,
                    trace_targets,
                    validate_ccd_region,
                    write_fits)


def fake_subprocess_popen(*args, stdout, stderr):
    raise OSError


def fake_dcr_communicate_stderr_no_dcr():
    return b'some message', b'dcr: not found'


class AddLinearWavelengthSolutionTest(TestCase):

    def setUp(self):
        self.ccd = CCDData(data=np.random.random_sample(200),
                           meta=fits.Header(),
                           unit='adu')
        self.ccd = add_wcs_keys(ccd=self.ccd)
        self.ccd.header.set('SLIT',
                            value='1.0_LONG_SLIT',
                            comment="slit [arcsec]")

    def test_add_wavelength_solution(self):

        calibration_lamp = 'non-existent.fits'

        crval1 = 3977.948
        npix = 4060
        cdelt = 0.9910068

        x_axis = np.linspace(crval1,
                             crval1 + cdelt * npix,
                             npix)

        self.ccd = add_linear_wavelength_solution(
            ccd=self.ccd,
            x_axis=x_axis,
            reference_lamp=calibration_lamp)

        self.assertEqual(self.ccd.header['CTYPE1'], 'LINEAR')
        self.assertEqual(self.ccd.header['CRVAL1'], crval1)
        self.assertEqual(self.ccd.header['CRPIX1'], 1)
        self.assertAlmostEqual(self.ccd.header['CDELT1'], cdelt, places=3)
        self.assertEqual(self.ccd.header['DCLOG1'],
                         'REFSPEC1 = {:s}'.format(calibration_lamp))


class AddWCSKeywordsTest(TestCase):

    def setUp(self):
        self.test_ccd = CCDData(data=np.ones((100, 100)),
                                meta=fits.Header(),
                                unit='adu')

    def test_add_wcs_keys(self):
        wcs_keys = ['BANDID1',
                    'APNUM1',
                    'WCSDIM',
                    'CTYPE1',
                    'CRVAL1',
                    'CRPIX1',
                    'CDELT1',
                    'CD1_1',
                    'LTM1_1',
                    'WAT0_001',
                    'WAT1_001',
                    'DC-FLAG',
                    'DCLOG1']



        self.test_ccd = add_wcs_keys(ccd=self.test_ccd)
        for key in wcs_keys:
            self.assertIn(key, self.test_ccd.header)

    @skip
    def test_add_wcs_keys_error(self):
        wcs_keys = ['BANDID1',
                    'APNUM1',
                    'WCSDIM',
                    'CTYPE1',
                    'CRVAL1',
                    'CRPIX1',
                    'CDELT1',
                    'CD1_1',
                    'LTM1_1',
                    'WAT0_001',
                    'WAT1_001',
                    'DC-FLAG',
                    'DCLOG1']


class BiasSubtractTest(TestCase):

    def setUp(self):
        self.ccd = CCDData(data=np.ones((100, 100)) * 100,
                           meta=fits.Header(),
                           unit='adu')
        self.master_bias = CCDData(data=np.ones((100, 100)) * 50,
                           meta=fits.Header(),
                           unit='adu')

        self.master_bias_name = os.path.join(os.getcwd(),
                                             'master_bias_file.fits')

    def test_bias_subtract(self):

        ccd = bias_subtract(ccd=self.ccd,
                            master_bias=self.master_bias,
                            master_bias_name=self.master_bias_name)
        np.testing.assert_array_equal(ccd.data, np.ones((100, 100)) * 50.)
        self.assertEqual(ccd.header['GSP_BIAS'],
                         os.path.basename(self.master_bias_name))


class BinningTest(TestCase):

    def test__bin_reference_data(self):
        wavelength = np.linspace(3000, 7000, 4000)
        intensity = np.random.random_sample(4000)

        for i in range(1, 4):
            new_wavelength, new_intensity = bin_reference_data(
                wavelength=wavelength,
                intensity=intensity,
                serial_binning=i)

            self.assertEqual(len(wavelength), len(intensity))
            self.assertEqual(len(new_wavelength), len(new_intensity))
            self.assertEqual(len(new_wavelength), np.floor(len(wavelength) / i))




class CentralWavelength(TestCase):

    def setUp(self):
        # 400m2
        self.grating = '400'
        self.grating_angle = 7.5
        self.camera_angle = 16.1
        self.reference_central_wavelength = 7001.54 * u.angstrom

    def test_get_central_wavelength(self):
        central_wavelength = get_central_wavelength(grating=self.grating,
                                                    grt_ang=self.grating_angle,
                                                    cam_ang=self.camera_angle)
        self.assertAlmostEqual(central_wavelength.value,
                               self.reference_central_wavelength.value,
                               places=2)


class ClassifySpectroscopicData(TestCase):

    def setUp(self):
        self.path = os.path.join(
            os.path.dirname(__file__),
            '../../data/test_data/test_classify_spectroscopic');
        if not os.path.isdir(self.path):
            os.mkdir(self.path)

    def tearDown(self):
        if os.path.isdir(self.path):
            shutil.rmtree(self.path)

    def create_fake_spectroscopic_data(self):
        if os.path.isdir(self.path):
            card_values = [
                {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:23:24.285', 'obsdec': '-39:12:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:24:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:24:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:24:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:25:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:25:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:25:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:26:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:27:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:28:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:23:34.285', 'obsdec': ' 39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:23:34.285', 'obsdec': ' 39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:23:34.285', 'obsdec': ' 39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:23:34.285', 'obsdec': ' 39:13:53.954'},
                {'obstype': 'OBJECT', 'object': 'NGC2070', 'obsra': '16:23:34.285', 'obsdec': ' 39:13:53.954'},
                {'obstype': 'COMP',   'object': 'HgArNe',  'obsra': '16:22:34.285', 'obsdec': '-39:23:53.954'},
                {'obstype': 'COMP',   'object': 'HgArNe',  'obsra': '16:24:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP',   'object': 'HgArNe',  'obsra': '16:24:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP',   'object': 'HgArNe',  'obsra': '16:24:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP',   'object': 'HgArNe',  'obsra': '16:25:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP',   'object': 'HgArNe',  'obsra': '16:25:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP',   'object': 'HgArNe',  'obsra': '16:25:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP',   'object': 'HgArNe',  'obsra': '16:26:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP',   'object': 'HgArNe',  'obsra': '16:27:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP',   'object': 'HgArNe',  'obsra': '16:28:34.285', 'obsdec': '-39:13:53.954'},
                {'obstype': 'COMP',   'object': 'HgArNe',  'obsra': '16:23:34.285', 'obsdec': '39:13:53.954'},
                {'obstype': 'COMP',   'object': 'HgArNe',  'obsra': '16:23:34.285', 'obsdec': '39:13:53.954'},
                {'obstype': 'COMP',   'object': 'HgArNe',  'obsra': '16:23:34.285', 'obsdec': '39:13:53.954'},
                {'obstype': 'COMP',   'object': 'HgArNe',  'obsra': '16:23:34.285', 'obsdec': '39:13:53.954'},
                {'obstype': 'COMP',   'object': 'HgArNe',  'obsra': '16:23:34.285', 'obsdec': '39:13:53.954'}]
            for i in range(len(card_values)):
                ccd = CCDData(data=np.ones((3, 3)),
                              meta=fits.Header(),
                              unit='adu')

                ccd.header.set('DATE', value='2019-03-22', comment='nc')
                ccd.header.set('SLIT', value='1.0" long slit', comment='nc')
                ccd.header.set('DATE-OBS', value='2019-03-22T09:59:33.654', comment='nc')
                ccd.header.set('OBSTYPE', value=card_values[i]['obstype'], comment='nc')
                ccd.header.set('OBJECT', value=card_values[i]['object'], comment='nc')
                ccd.header.set('EXPTIME', value='10', comment='nc')
                ccd.header.set('OBSRA', value=card_values[i]['obsra'], comment='nc')
                ccd.header.set('OBSDEC', value=card_values[i]['obsdec'], comment='nc')
                ccd.header.set('GRATING', value='SYZY_400', comment='nc')
                ccd.header.set('CAM_TARG', value='16.1', comment='nc')
                ccd.header.set('GRT_TARG', value='7.5', comment='nc')
                ccd.header.set('FILTER', value='<NO FILTER>', comment='nc')
                ccd.header.set('FILTER2', value='GG455', comment='nc')
                ccd.header.set('GAIN', value='1.48', comment='nc')
                ccd.header.set('RDNOISE', value='3.89', comment='nc')

                ccd.write(os.path.join(self.path, 'test_file_{:03d}.fits'.format(i)))


    def test_classify_spectroscopic_data__no_data(self):
        self.assertRaises(SystemExit, classify_spectroscopic_data, self.path, '*fits')

    def test_classify_spectroscopic_data(self):
        self.create_fake_spectroscopic_data()
        result = classify_spectroscopic_data(path=self.path, search_pattern='test_file')
        self.assertIsInstance(result, NightDataContainer)
        self.assertFalse(result.is_empty)


class CombineDataTest(TestCase):

    def setUp(self):
        self.ccd1 = CCDData(data=np.ones((100, 100)),
                            meta=fits.Header(),
                            unit='adu')
        self.ccd1.header.set('OBJECT', value='TestObject')
        self.ccd1.header.set('GRATING', value='Grating')
        self.ccd1.header.set('SLIT', value='1.05SlitSize')

        self.ccd2 = self.ccd1.copy()
        self.ccd2.data *= 2
        self.ccd3 = self.ccd1.copy()
        self.ccd3.data *= 5
        self.ccd1.header.set('GSP_FNAM', value='prefix_0001_image.fits')
        self.ccd2.header.set('GSP_FNAM', value='prefix_0002_image.fits')
        self.ccd3.header.set('GSP_FNAM', value='prefix_0003_image.fits')
        self.image_list = [self.ccd1, self.ccd2, self.ccd3]
        self.dest_path = os.getcwd()
        self.prefix = 'testing_'
        self.output_name = 'combine_data.fits'

    def tearDown(self):
        using_output_name_file_name = os.path.join(
            self.dest_path,
            self.output_name)
        if os.path.isfile(using_output_name_file_name):
            os.unlink(using_output_name_file_name)

        not_using_output_name_file_name_list = os.listdir(self.dest_path)
        if not_using_output_name_file_name_list:
            for _file in not_using_output_name_file_name_list:
                if '{:s}comb_'.format(self.prefix) in _file:
                    os.unlink(_file)

    def test_combine_data_median_prefix_ignored(self):
        combined = combine_data(
            image_list=self.image_list,
            dest_path=self.dest_path,
            prefix=self.prefix,
            output_name=self.output_name,
            method='median',
            save=True)

        np.testing.assert_array_equal(combined.data, np.ones((100, 100)) * 1.5)
        self.assertEqual(len(combined.header['GSP_IC*']), 3)
        self.assertTrue(self.prefix not in combined.header['GSP_FNAM'])

    def test_combine_data_median_prefix_used(self):
        combined = combine_data(
            image_list=self.image_list,
            dest_path=self.dest_path,
            prefix=self.prefix,
            method='median',
            save=True)

        np.testing.assert_array_equal(combined.data, np.ones((100, 100)) * 1.5)
        self.assertEqual(len(combined.header['GSP_IC*']), 3)
        self.assertTrue(self.prefix in combined.header['GSP_FNAM'])


class CosmicRayRejectionTest(TestCase):

    def setUp(self):
        self.ccd = CCDData(data=np.ones((100, 100)),
                           meta=fits.Header(),
                           unit='adu')
        self.file_name = 'cr_test.fits'

        self.ccd.header.set('CCDSUM', value='1 1')
        self.ccd.header.set('OBSTYPE', value='OBJECT')
        self.ccd.header.set('INSTCONF', value='Red')
        self.ccd.header.set('GSP_FNAM', value=self.file_name)
        self.ccd.header.set('GSP_COSM', value='none')

        self.red_path = os.getcwd()
        self.out_prefix = 'prefix'

    @skip
    def test_dcr_cosmicray_rejection(self):
        pass

    def test_call_cosmic_rejection_default_1x1(self):
        prefix = 'new_'
        initial_value = self.ccd.data[50, 50]
        self.ccd.data[50, 50] = 50000

        ccd, out_prefix = call_cosmic_rejection(ccd=self.ccd,
                                                image_name=self.file_name,
                                                out_prefix=self.out_prefix,
                                                red_path=self.red_path,
                                                keep_files=True,
                                                prefix=prefix,
                                                method='default',
                                                save=False)
        self.assertAlmostEqual(initial_value, ccd.data[50, 50])
        self.assertEqual(prefix + self.out_prefix, out_prefix)
        self.assertEqual(f"{prefix}{self.out_prefix}_{self.file_name}",
                         ccd.header['GSP_FNAM'])
        self.assertEqual('DCR', ccd.header['GSP_COSM'])

        self.assertTrue(os.path.isfile('dcr.par'))

    @mock.patch('subprocess.Popen', side_effect=fake_subprocess_popen)
    def test_call_cosmic_rejection_default_1x1_no_dcr_par(self,
                                                          subprocess_Popen_function):

        self.assertRaises(SystemExit,
                          call_cosmic_rejection,
                          self.ccd,
                          self.file_name,
                          self.out_prefix,
                          self.red_path,
                          os.getcwd())

    @mock.patch('subprocess.Popen.communicate',
                side_effect=fake_dcr_communicate_stderr_no_dcr)
    def test_dcr_cosmicray_rejection_no_dcr_executable(
            self, subprocess_Popen_communicate):

        self.assertRaises(SystemExit,
                          dcr_cosmicray_rejection,
                          self.red_path,
                          self.file_name,
                          'c',
                          os.getcwd())


    def test_call_cosmic_rejection_default_2x2(self):
        self.ccd.header.set('CCDSUM', value='2 2')
        prefix = 'new_'
        initial_value = self.ccd.data[50, 50]
        self.ccd.data[50, 50] = 50000

        ccd, out_prefix = call_cosmic_rejection(ccd=self.ccd,
                                                image_name=self.file_name,
                                                out_prefix=self.out_prefix,
                                                red_path=self.red_path,
                                                keep_files=True,
                                                prefix=prefix,
                                                method='default',
                                                save=True)
        self.assertAlmostEqual(initial_value, ccd.data[50, 50])
        self.assertEqual(prefix + self.out_prefix, out_prefix)
        self.assertEqual(f"{prefix}{self.out_prefix}_{self.file_name}",
                         ccd.header['GSP_FNAM'])
        self.assertEqual('LACosmic', ccd.header['GSP_COSM'])
        self.assertTrue(os.path.isfile('new_prefix_cr_test.fits'))

    def test_call_cosmic_rejection_default_3x3(self):
        self.ccd.header.set('CCDSUM', value='3 3')
        prefix = 'new_'
        initial_value = self.ccd.data[50, 50]
        self.ccd.data[50, 50] = 50000

        ccd, out_prefix = call_cosmic_rejection(ccd=self.ccd,
                                                image_name=self.file_name,
                                                out_prefix=self.out_prefix,
                                                red_path=self.red_path,
                                                keep_files=True,
                                                prefix=prefix,
                                                method='default',
                                                save=True)
        self.assertAlmostEqual(initial_value, ccd.data[50, 50])
        self.assertEqual(prefix + self.out_prefix, out_prefix)
        self.assertEqual(f"{prefix}{self.out_prefix}_{self.file_name}",
                         ccd.header['GSP_FNAM'])
        self.assertEqual('LACosmic', ccd.header['GSP_COSM'])

        self.assertTrue(os.path.isfile('new_prefix_cr_test.fits'))

    def test_call_cosmic_rejection_none(self):
        prefix = 'new_'
        ccd, out_prefix = call_cosmic_rejection(ccd=self.ccd,
                                                image_name=self.file_name,
                                                out_prefix=self.out_prefix,
                                                red_path=self.red_path,
                                                keep_files=True,
                                                prefix=prefix,
                                                method='none',
                                                save=True)
        self.assertEqual(self.out_prefix, out_prefix)
        self.assertEqual(f"{self.out_prefix}_{self.file_name}",
                         ccd.header['GSP_FNAM'])
        self.assertEqual('none', ccd.header['GSP_COSM'])
        self.assertTrue(os.path.isfile('prefix_cr_test.fits'))

    def test_call_cosmic_rejection_comp_lamp(self):
        self.ccd.header.set('OBSTYPE', value='COMP')
        prefix = 'new_'
        ccd, out_prefix = call_cosmic_rejection(ccd=self.ccd,
                                                image_name=self.file_name,
                                                out_prefix=self.out_prefix,
                                                red_path=self.red_path,
                                                keep_files=True,
                                                prefix=prefix,
                                                method='lacosmic',
                                                save=True)
        self.assertEqual(prefix + self.out_prefix, out_prefix)
        self.assertEqual(f"{prefix}{self.out_prefix}_{self.file_name}",
                         ccd.header['GSP_FNAM'])
        self.assertEqual('none', ccd.header['GSP_COSM'])

    def test_call_cosmic_rejection_arc_lamp(self):
        self.ccd.header.set('OBSTYPE', value='ARC')
        prefix = 'new_'
        ccd, out_prefix = call_cosmic_rejection(ccd=self.ccd,
                                                image_name=self.file_name,
                                                out_prefix=self.out_prefix,
                                                red_path=self.red_path,
                                                keep_files=True,
                                                prefix=prefix,
                                                method='dcr',
                                                save=True)
        self.assertEqual(prefix + self.out_prefix, out_prefix)
        self.assertEqual(f"{prefix}{self.out_prefix}_{self.file_name}",
                         ccd.header['GSP_FNAM'])
        self.assertEqual('none', ccd.header['GSP_COSM'])

    def test_call_cosmic_rejection_not_implemented_error(self):
        prefix = 'new_'
        self.assertRaises(NotImplementedError,
                          call_cosmic_rejection,
                          self.ccd,
                          self.file_name,
                          self.out_prefix,
                          self.red_path,
                          True,
                          prefix,
                          'not_implemented_method',
                          True)

    def tearDown(self):
        files_to_delete = ['dcr.par',
                           'goodman_log.txt',
                           'cosmic_cr_test.fits',
                           'new_prefix_cr_test.fits',
                           'prefix_cr_test.fits',
                           'crmask_cr_test.fits']

        for _file in files_to_delete:
            if os.path.isfile(_file):
                os.unlink(_file)


class CreateMasterBias(TestCase):

    def setUp(self):
        self.name = ''
        self.bias_files = ['bias_{}.fits'.format(i) for i in range(1, 12)]
        self.raw_data = os.getcwd()
        self.reduced_data = os.getcwd()

        self.overscan_region = '[1:10,1:100]'
        self.trim_section = '[11:90,11:90]'

        self.ccd = CCDData(data=np.ones((100, 100)),
                           meta=fits.Header(),
                           unit='adu')

        self.ccd.data[11:90, 10:90] += 1
        self.ccd.header.set('CCDSUM', value='1 1')
        self.ccd.header.set('INSTCONF', value='Red')
        self.ccd.header.set('RDNOISE', value=3.89)
        self.ccd.header.set('GAIN', value=1.48)

        for _file in self.bias_files:
            self.ccd.write(os.path.join(self.reduced_data, _file))

    def tearDown(self):
        for _file in self.bias_files:
            os.unlink(os.path.join(self.reduced_data, _file))

        if self.name != '':
            for _file in [self.name]:
                os.unlink(_file)

    def test_create_master_bias(self):
        master, self.name = create_master_bias(
            bias_files=self.bias_files,
            raw_data=self.raw_data,
            reduced_data=self.reduced_data,
            technique='Spectroscopy')

        self.assertTrue('master_bias' in self.name)
        self.assertEqual('master_bias_RED_SP_1x1_R03.89_G01.48.fits', self.name)

        self.assertTrue(all(
            [master.header[key] in self.bias_files for key in master.header['GSP_IC*'].keys()]))


class CreateMasterFlatsTest(TestCase):

    def setUp(self):
        self.flat_files = ['flat_{}.fits'.format(i) for i in range(1, 12)]
        self.raw_data = os.getcwd()
        self.reduced_data = os.getcwd()
        self.master_flat_name = 'master_flat.fits'
        self.saturation_level = 69257

        self.overscan_region = '[1:10,1:100]'
        self.trim_section = '[11:90,11:90]'

        self.ccd = CCDData(data=np.ones((100, 100)),
                           meta=fits.Header(),
                           unit='adu')

        self.ccd.data[11:90, 10:90] += 5
        self.ccd.header.set('CCDSUM', value='1 1')
        self.ccd.header.set('INSTCONF', value='Red')
        self.ccd.header.set('GAIN', value=1.48)
        self.ccd.header.set('RDNOISE', value=3.89)

        for _file in self.flat_files:
            self.ccd.write(os.path.join(self.reduced_data, _file))

        self.bias = CCDData(data=np.ones((100, 100)),
                            meta=fits.Header(),
                            unit='adu')
        self.bias.header.set('CCDSUM', value='1 1')
        self.bias.write(os.path.join(self.reduced_data, 'master_bias.fits'))

    def tearDown(self):
        for _file in self.flat_files:
            os.unlink(os.path.join(self.reduced_data, _file))
        for _file in ['master_bias.fits', self.master_flat_name]:
            if os.path.isfile(_file):
                os.unlink(_file)

    def test_create_master_flats_no_bias(self):
        master, name = create_master_flats(
            flat_files=self.flat_files,
            raw_data=self.raw_data,
            reduced_data=self.reduced_data,
            technique='Spectroscopy',
            overscan_region=self.overscan_region,
            trim_section=self.trim_section,
            master_bias_name='master_bias.fits',
            new_master_flat_name=self.master_flat_name,
            saturation_threshold=1,
            ignore_bias=True)

        self.assertEqual(self.master_flat_name, os.path.basename(name))
        self.assertTrue(os.path.isfile(name))
        self.assertEqual('none', master.header['GSP_BIAS'])

    def test_create_master_flats_with_bias(self):

        master, name = create_master_flats(
            flat_files=self.flat_files,
            raw_data=self.raw_data,
            reduced_data=self.reduced_data,
            technique='Spectroscopy',
            overscan_region=self.overscan_region,
            trim_section=self.trim_section,
            master_bias_name=os.path.join(self.reduced_data,
                                          'master_bias.fits'),
            new_master_flat_name=os.path.join(self.reduced_data,
                                              self.master_flat_name),
            saturation_threshold=1,
            ignore_bias=False)

        self.assertEqual(self.master_flat_name, os.path.basename(name))
        self.assertTrue(os.path.isfile(name))
        self.assertEqual('master_bias.fits', master.header['GSP_BIAS'])


    def test_create_master_flats_saturated_flats(self):

        file_to_replace = os.path.join(self.reduced_data, self.flat_files[0])
        os.unlink(file_to_replace)
        self.assertFalse(os.path.isfile(file_to_replace))
        self.ccd.data[11:90, 10:90] += 2 * self.saturation_level
        self.ccd.header.set('OBJECT', value='saturated')
        self.ccd.write(os.path.join(self.reduced_data, self.flat_files[0]),
                       overwrite=True)

        master, name = create_master_flats(
            flat_files=self.flat_files,
            raw_data=self.raw_data,
            reduced_data=self.reduced_data,
            technique='Spectroscopy',
            overscan_region=self.overscan_region,
            trim_section=self.trim_section,
            master_bias_name='master_bias.fits',
            new_master_flat_name=self.master_flat_name,
            saturation_threshold=1,
            ignore_bias=False)

        self.assertNotEqual(self.flat_files[0], master.header['GSP_IC01'])

    def test_create_master_flats_empty_list(self):
        master, name = create_master_flats(
            flat_files=[],
            raw_data=self.raw_data,
            reduced_data=self.reduced_data,
            technique='Spectroscopy',
            overscan_region=self.overscan_region,
            trim_section=self.trim_section,
            master_bias_name='master_bias.fits',
            new_master_flat_name=self.master_flat_name,
            saturation_threshold=1,
            ignore_bias=False)
        self.assertIsNone(master)
        self.assertIsNone(name)


class CrossCorrelationTest(TestCase):

    def setUp(self):
        self.binning = 1
        self.size = 5000
        self.reference_array = np.ones(self.size)
        self.compared_array = np.ones(self.size)
        self.x_axis = np.arange(0, self.size, 1)
        self.gaussian = models.Gaussian1D(stddev=5, amplitude=3000)
        self.gaussian.mean.value = int(self.size / 2.)
        self.reference_array += self.gaussian(self.x_axis)

    def test_cross_correlation_small_slit(self):

        offset = 500

        self.gaussian.mean.value -= offset

        self.compared_array += self.gaussian(self.x_axis)

        correlation_value = cross_correlation(reference=self.reference_array,
                                              compared=self.compared_array,
                                              slit_size=1,
                                              serial_binning=1,
                                              mode='full',
                                              plot=False)

        self.assertEqual(offset, correlation_value)

    def test_cross_correlation_large_slit(self):
        offset = 500

        self.gaussian.mean.value -= offset

        self.compared_array += self.gaussian(self.x_axis)

        box_kernel = Box1DKernel(width=5 / 0.15)

        self.compared_array = convolve(self.compared_array, box_kernel)

        correlation_value = cross_correlation(reference=self.reference_array,
                                              compared=self.compared_array,
                                              slit_size=5,
                                              serial_binning=1,
                                              mode='full',
                                              plot=False)

        self.assertEqual(offset, correlation_value)


class DefineTrimSectionTest(TestCase):

    def setUp(self):
        self.ccd = CCDData(data=np.ones((100, 100)),
                           meta=fits.Header(),
                           unit='adu')
        self.ccd.header.set('CCDSUM', value='1 1')
        self.ccd.header.set('TRIMSEC', value='[10:10,10:10]')

        self.full_path = os.path.join(os.getcwd(), 'testfile.fits')

        self.ccd.write(self.full_path)

    def test_define_trim_section_spectroscopy(self):

        expected_trim_section = '[51:4110,2:100]'
        trim_section = define_trim_section(sample_image=self.full_path,
                                           technique='Spectroscopy')
        self.assertEqual(expected_trim_section, trim_section)

    def test_define_trim_section_spectroscopy_2x2(self):
        self.ccd.header.set('CCDSUM', value='2 2')
        self.ccd.write(self.full_path, overwrite=True)

        expected_trim_section = '[26:2055,2:100]'
        trim_section = define_trim_section(sample_image=self.full_path,
                                           technique='Spectroscopy')
        self.assertEqual(expected_trim_section, trim_section)

    def test_define_trim_section_imaging(self):

        expected_trim_section = '[10:10,10:10]'
        trim_section = define_trim_section(sample_image=self.full_path,
                                           technique='Imaging')

        self.assertEqual(expected_trim_section, trim_section)

    def tearDown(self):
        if os.path.isfile(self.full_path):
            os.unlink(self.full_path)

class EvaluateWavelengthSolutionTest(TestCase):

    def test__evaluate_solution(self):
        differences = np.array([0.5] * 10)

        clipped_differences = np.ma.masked_array(differences,
                                                 mask=[0,
                                                       0,
                                                       1,
                                                       0,
                                                       0,
                                                       1,
                                                       0,
                                                       0,
                                                       1,
                                                       0])

        rms_error, n_points, n_rej = evaluate_wavelength_solution(
            clipped_differences=clipped_differences)

        self.assertEqual(rms_error, 0.5)
        self.assertEqual(n_points, 10)
        self.assertEqual(n_rej, 3)


class ExtractionTest(TestCase):

    def setUp(self):
        self.fake_image = CCDData(data=np.ones((100, 100)),
                                  meta=fits.Header(),
                                  unit='adu')
        self.fake_image.mask = np.zeros((100, 100), dtype=bool)
        self.fake_image.header.set('NAXIS', value=2)
        self.fake_image.header.set('NAXIS1', value=100)
        self.fake_image.header.set('NAXIS2', value=100)
        self.fake_image.header.set('OBSTYPE', value='COMP')
        self.fake_image.header['GSP_FNAM'] = 'fake-image.fits'

        # Create model aligned with pixels - represents the trace
        self.target_trace = models.Linear1D(slope=0, intercept=50.3)

        # Calculate the STDDEV
        self.stddev = 8.4

        # Calculate how many STDDEV will be extracted - N_STDDEV
        self.n_stddev = 2

        # Calculate how far the background is from the the center.
        self.distance = 1

        self.target_profile_gaussian = models.Gaussian1D(amplitude=1,
                                                         mean=50.3,
                                                         stddev=self.stddev)

        self.target_profile_moffat = models.Moffat1D(amplitude=1,
                                                     x_0=50.3,
                                                     gamma=self.stddev)

        self.reference_result_gaussian = np.ones(
            100) * self.target_profile_gaussian.fwhm * self.n_stddev

        self.reference_result_moffat = np.ones(
            100) * self.target_profile_moffat.fwhm * self.n_stddev

    def test_fractional_extraction(self):
        # Perform extraction
        extracted_array, background, info = extract_fractional_pixel(
            ccd=self.fake_image,
            target_trace=self.target_trace,
            target_fwhm=self.target_profile_gaussian.fwhm,
            extraction_width=self.n_stddev,
            background_spacing=self.distance)
        # assert isinstance(fake_image, CCDData)
        self.assertIsInstance(extracted_array, CCDData)

        np.testing.assert_array_almost_equal(extracted_array,
                                             self.reference_result_gaussian)

    def test_fractional_extraction_obstype_object(self):
        self.fake_image.header.set('OBSTYPE', value='OBJECT')
        # Perform extraction
        extracted_array, background, info = extract_fractional_pixel(
            ccd=self.fake_image,
            target_trace=self.target_trace,
            target_fwhm=self.stddev,
            extraction_width=self.n_stddev,
            background_spacing=self.distance)
        # assert isinstance(fake_image, CCDData)
        self.assertIsInstance(extracted_array, CCDData)

        np.testing.assert_array_almost_equal(extracted_array,
                                             np.zeros(extracted_array.shape))

    def test_fractional_sum(self):

        fake_image = np.ones((100, 100))
        low_limit = 50 + np.random.random()
        high_limit = 60 + np.random.random()

        frac_sum = fractional_sum(fake_image, 50, low_limit, high_limit)
        self.assertEqual(frac_sum, high_limit - low_limit)

    def test_extract_optimal(self):
        self.assertRaises(NotImplementedError, extract_optimal)

    def test_extract__optimal_not_implemented(self):
        self.assertRaises(NotImplementedError,
                          extraction,
                          self.fake_image,
                          self.target_trace,
                          self.target_profile_gaussian,
                          'optimal')

    def test_extraction_gaussian(self):
        extracted = extraction(ccd=self.fake_image,
                               target_trace=self.target_trace,
                               spatial_profile=self.target_profile_gaussian,
                               extraction_name='fractional')
        self.assertIsInstance(extracted, CCDData)
        np.testing.assert_array_almost_equal(extracted, self.reference_result_gaussian)

    def test_extraction_moffat(self):
        spatial_profile_moffat = models.Moffat1D(amplitude=self.target_profile_gaussian.amplitude.value,
                                                 x_0=self.target_profile_gaussian.mean.value,
                                                 gamma=self.target_profile_gaussian.stddev.value)
        extracted = extraction(ccd=self.fake_image,
                               target_trace=self.target_trace,
                               spatial_profile=spatial_profile_moffat,
                               extraction_name='fractional')
        self.assertIsInstance(extracted, CCDData)
        np.testing.assert_array_almost_equal(extracted, self.reference_result_moffat)

    def test_extraction_not_implemented_model(self):
        spatial_profile = models.BlackBody()
        self.assertRaises(NotImplementedError, extraction, ccd=self.fake_image,
                          target_trace=self.target_trace,
                          spatial_profile=spatial_profile,
                          extraction_name='fractional')

    def test_extraction_exception(self):
        self.assertRaises(NotImplementedError, extraction, ccd=self.fake_image,
                          target_trace=self.target_trace,
                          spatial_profile=self.target_profile_gaussian,
                          extraction_name='optimal')


class FitsFileIOAndOps(TestCase):

    def setUp(self):
        self.fake_image = CCDData(data=np.ones((100, 100)),
                                  meta=fits.Header(),
                                  unit='adu')

        self.file_name = 'sample_file.fits'
        self.target_non_zero = 4
        self.current_directory = os.getcwd()
        self.full_path = os.path.join(self.current_directory, self.file_name)
        self.parent_file = 'parent_file.fits'

        self.fake_image.header.set('CCDSUM',
                                   value='1 1',
                                   comment='Fake values')

        self.fake_image.header.set('OBSTYPE',
                                   value='OBJECT',
                                   comment='Fake values')

        self.fake_image.header.set('GSP_FNAM',
                                   value=self.file_name,
                                   comment='Fake values')

        self.fake_image.header.set('GSP_PNAM',
                                   value=self.parent_file,
                                   comment='Fake values')

        self.fake_image.write(self.full_path, overwrite=False)

    def test_write_fits(self):
        self.assertTrue(os.path.isfile(self.full_path))
        os.remove(self.full_path)
        write_fits(ccd=self.fake_image,
                   full_path=self.full_path,
                   parent_file=self.parent_file,
                   overwrite=False)
        self.assertTrue(os.path.isfile(self.full_path))

    def test_write_fits_directory_does_not_exist(self):
        os.unlink(self.full_path)

        base_path = os.path.dirname(self.full_path)
        self.full_path = os.path.join(base_path,
                                      'testdir',
                                      os.path.basename(self.full_path))
        self.assertFalse(os.path.isdir(os.path.dirname(self.full_path)))
        write_fits(ccd=self.fake_image,
                   full_path=self.full_path,
                   parent_file=self.parent_file,
                   overwrite=True)
        self.assertTrue(os.path.isdir(os.path.dirname(self.full_path)))


    def test_read_fits(self):
        self.fake_image.header.remove('GSP_PNAM')
        self.fake_image.write(self.full_path, overwrite=True)

        self.recovered_fake_image = read_fits(self.full_path)
        self.assertIsInstance(self.recovered_fake_image, CCDData)

    def test_image_overscan(self):
        data_value = 100.
        overscan_value = 0.1
        # alter overscan region to a lower number
        self.fake_image.data *= data_value
        self.fake_image.data[:, 0:5] = overscan_value

        overscan_region = '[1:6,:]'
        self.assertEqual(self.fake_image.data[:, 6:99].mean(), data_value)
        self.assertEqual(self.fake_image.data[:, 0:5].mean(), overscan_value)
        self.fake_image = image_overscan(ccd=self.fake_image,
                                         overscan_region=overscan_region)

        self.assertEqual(self.fake_image.data[:, 6:99].mean(),
                         data_value - overscan_value)
        self.assertEqual(self.fake_image.header['GSP_OVER'], overscan_region)

    def test_image_overscan_none(self):
        new_fake_image = image_overscan(ccd=self.fake_image,
                                        overscan_region=None)
        self.assertEqual(new_fake_image, self.fake_image)

    def test_image_trim(self):
        self.assertEqual(self.fake_image.data.shape, (100, 100))
        trim_section = '[1:50,:]'
        self.fake_image = image_trim(ccd=self.fake_image,
                                     trim_section=trim_section,
                                     trim_type='trimsec')

        self.assertEqual(self.fake_image.data.shape, (100, 50))
        self.assertEqual(self.fake_image.header['GSP_TRIM'], trim_section)

    def test_image_trim_invalid_type(self):
        self.assertEqual(self.fake_image.data.shape, (100, 100))
        trim_section = '[1:50,:]'
        self.fake_image = image_trim(ccd=self.fake_image,
                                     trim_section=trim_section,
                                     trim_type='invalid_type')
        self.assertEqual(self.fake_image.data.shape, (100, 50))
        self.assertEqual(self.fake_image.header['GSP_TRIM'], trim_section)

    def test_image_trim_trim_section_none(self):
        self.assertEqual(self.fake_image.data.shape, (100, 100))
        self.fake_image = image_trim(ccd=self.fake_image,
                                     trim_section=None,
                                     trim_type='trimsec')
        self.assertEqual(self.fake_image.data.shape, (100, 100))

    def test_save_extracted_target_zero(self):
        self.fake_image.header.set('GSP_FNAM', value=self.file_name)
        same_fake_image = save_extracted(ccd=self.fake_image,
                                         destination=self.current_directory,
                                         prefix='e',
                                         target_number=0)

        assert np.array_equal(same_fake_image.data, self.fake_image.data)
        self.assertEqual(same_fake_image.header, self.fake_image.header)
        self.assertEqual(same_fake_image.shape, self.fake_image.shape)
        self.assertTrue(os.path.isfile('e' + self.file_name))

    def test_save_extracted_target_non_zero(self):
        self.fake_image.header.set('GSP_FNAM', value=self.file_name)
        same_fake_image = save_extracted(ccd=self.fake_image,
                                         destination=self.current_directory,
                                         prefix='e',
                                         target_number=self.target_non_zero)
        assert np.array_equal(same_fake_image.data, self.fake_image.data)
        self.assertEqual(same_fake_image.header, self.fake_image.header)
        self.assertEqual(same_fake_image.shape, self.fake_image.shape)
        self.assertTrue(os.path.isfile('e' + re.sub('.fits',
                                       '_target_{:d}.fits'.format(
                                           self.target_non_zero),
                                       self.file_name)))

    def test_save_extracted_target_zero_comp(self):
        self.fake_image.header.set('GSP_FNAM', value=self.file_name)
        self.fake_image.header.set('OBSTYPE', value='ARC')
        self.fake_image.header.set('GSP_EXTR', value='100.00:101.00')
        same_fake_image = save_extracted(ccd=self.fake_image,
                                         destination=self.current_directory,
                                         prefix='e',
                                         target_number=0)
        expected_new_name = 'e' + re.sub(
            '.fits',
            '_' + re.sub(':', '-', same_fake_image.header['GSP_EXTR']) + '.fits',
            self.file_name)

        assert np.array_equal(same_fake_image.data, self.fake_image.data)
        self.assertEqual(same_fake_image.header, self.fake_image.header)
        self.assertEqual(same_fake_image.shape, self.fake_image.shape)
        self.assertEqual(same_fake_image.header['GSP_FNAM'], expected_new_name)
        self.assertTrue(os.path.isfile(self.fake_image.header['GSP_FNAM']))

    def tearDown(self):
        files_to_remove = [self.full_path, self.fake_image.header['GSP_FNAM']]

        for _file in files_to_remove:
            if os.path.isfile(_file):
                os.unlink(_file)

        if 'testdir' in self.full_path:
            os.rmdir(os.path.dirname(self.full_path))


class FixKeywordsTest(TestCase):

    def setUp(self):
        self.ccd = CCDData(data=np.ones((100, 100)),
                           meta=fits.Header(),
                           unit='adu')
        self.full_path = os.path.join(os.getcwd(), 'sample_file.fits')
        self.ccd.write(self.full_path)

    def tearDown(self):
        if os.path.isfile(self.full_path):
            os.unlink(self.full_path)

    def test_fix_keywords(self):
        # not really testing anything here
        fix_keywords(path=os.getcwd(), pattern="*.fits")


class GenerateDcrFile(TestCase):

    def setUp(self):
        self.create = GenerateDcrParFile()
        self.ccd = CCDData(data=np.ones((100, 100)),
                           meta=fits.Header(),
                           unit='adu')
        self.ccd.header.set('INSTCONF', value='Red')
        self.ccd.header.set('CCDSUM', value='1 1')

    def test_generate_dcr_par_file(self):
        serial, parallel = self.ccd.header['CCDSUM'].split()
        instrument = self.ccd.header['INSTCONF']

        self.assertEqual(serial, '1')
        self.assertEqual(instrument, 'Red')
        self.assertEqual(self.create._file_name, 'dcr.par')
        self.assertIsInstance(self.create._df, pandas.DataFrame)

        self.assertFalse(os.path.isfile(self.create._file_name))
        self.create()
        self.assertTrue(os.path.isfile(self.create._file_name))

        self.assertRaises(AssertionError, self.create, 'Green')

    def tearDown(self):
        if os.path.isfile(self.create._file_name):
            os.remove(self.create._file_name)


class GetLinesInLampTest(TestCase):

    def setUp(self):
        self.ccd = CCDData(data=np.ones(5000),
                           meta=fits.Header(),
                           unit='adu')
        self.ccd.header.set('SLIT', value='1.0_LONG_SLIT')
        self.ccd.header.set('CCDSUM', value='1 1')
        self.ccd.header.set('OBJECT', value='TestSubject')
        self.line_centers = np.arange(100, 4900, 200)

        self.gaussian = models.Gaussian1D(amplitude=10000, stddev=3)

        for line_center in self.line_centers:
            self.gaussian.mean.value = line_center
            self.ccd.data += self.gaussian(range(len(self.ccd.data)))

    def test_get_lines_in_lamp_wrong_input(self):

        expected_none = get_lines_in_lamp(ccd=list(range(100)))

        self.assertIsNone(expected_none)


    def test_get_lines_in_lamp_narrow_slit(self):
        recovered_lines = get_lines_in_lamp(ccd=self.ccd, plots=False)
        np.testing.assert_allclose(self.line_centers, recovered_lines)

    def test_get_lines_in_lamp_broad_slit(self):
        self.ccd.header.set('SLIT', value='5.0_LONG_SLIT')

        box_kernel = Box1DKernel(width=5.0 / 0.15)
        self.ccd.data = convolve(self.ccd.data, box_kernel)

        recovered_lines = get_lines_in_lamp(ccd=self.ccd, plots=False)
        np.testing.assert_allclose(self.line_centers, recovered_lines, atol=0.6)


class GetOverscanRegionTest(TestCase):

    def setUp(self):
        self.ccd = CCDData(data=np.ones((100, 100)),
                           meta=fits.Header(),
                           unit='adu')
        self.ccd.header.set('CCDSUM', value='1 1')
        self.ccd.header.set('TRIMSEC', value='[10:10,10:10]')

        self.full_path = os.path.join(os.getcwd(), 'testfile.fits')

        self.ccd.write(self.full_path)

    def tearDown(self):
        if os.path.isfile(self.full_path):
            os.unlink(self.full_path)

    def test_get_overscan_region_spectroscopy_blue(self):
        self.ccd.header.set('INSTCONF', 'Blue')
        self.ccd.write(self.full_path, overwrite=True)

        expected_overscan = '[1:16,1:100]'

        overscan_region = get_overscan_region(sample_image=self.full_path,
                                              technique='Spectroscopy')

        self.assertEqual(expected_overscan, overscan_region)

    def test_get_overscan_region_spectroscopy_red(self):
        self.ccd.header.set('INSTCONF', 'Red')
        self.ccd.write(self.full_path, overwrite=True)

        expected_overscan = '[6:49,1:100]'

        overscan_region = get_overscan_region(sample_image=self.full_path,
                                              technique='Spectroscopy')

        self.assertEqual(expected_overscan, overscan_region)

    def test_get_overscan_region_imaging(self):
        overscan_region = get_overscan_region(sample_image=self.full_path,
                                              technique='Imaging')

        self.assertIsNone(overscan_region)

    def test_get_overscan_region_other_technique(self):

        overscan_region = get_overscan_region(sample_image=self.full_path,
                                              technique='AnyOther')

        self.assertIsNone(overscan_region)


class GetSpectralCharacteristicsTest(TestCase):
    def setUp(self):
        self.ccd = CCDData(data=np.random.random_sample(200),
                           meta=fits.Header(),
                           unit='adu')
        self.pixel_size = 15 * u.micrometer
        self.goodman_focal_length = 377.3 * u.mm

    def test__get_spectral_characteristics(self):
        self.ccd.header.set('SLIT', '1.0_LONG_SLIT')
        self.ccd.header.set('GRATING', 'SYZY_400')
        self.ccd.header.set('GRT_ANG', 7.5)
        self.ccd.header.set('CAM_ANG', 16.1)
        self.ccd.header.set('CCDSUM', '1 1')

        spec_charact = get_spectral_characteristics(
            ccd=self.ccd,
            pixel_size=self.pixel_size,
            instrument_focal_length=self.goodman_focal_length)

        self.assertIsInstance(spec_charact, dict)
        self.assertEqual(len(spec_charact), 7)


class InterpolationTest(TestCase):

    def test_interpolate_spectrum(self):
        initial_array = np.sin(np.arange(0, 3 * np.pi))
        initial_length = len(initial_array)

        new_x_axis, new_array = interpolate_spectrum(spectrum=initial_array,
                                                     interpolation_size=100)

        self.assertEqual(len(new_x_axis), len(new_array))
        self.assertEqual(len(new_array), initial_length * 100)


class IsFileSaturatedTest(TestCase):

    def setUp(self):
        self.ccd = CCDData(data=np.ones((100, 100)),
                           meta=fits.Header(),
                           unit='adu')
        self.ccd.header.set('INSTCONF', value='Red')
        self.ccd.header.set('GAIN', value=1.48)
        self.ccd.header.set('RDNOISE', value=3.89)

        self.half_full_well = 69257

    def test_file_is_saturated(self):
        self.ccd.data[:10, :10] = self.half_full_well + 1
        self.assertTrue(is_file_saturated(ccd=self.ccd, threshold=1))

    def test_file_is_not_saturated(self):
        self.ccd.data[:10, :10] = self.half_full_well + 1
        self.ccd.data[0, 0] = 1
        self.assertFalse(is_file_saturated(ccd=self.ccd, threshold=1))


class LinearizeSpectrumTest(TestCase):

    def setUp(self):
        feature = models.Gaussian1D(amplitude=500, mean=3000, stddev=5)
        self.data = np.ones(5000) + feature(range(5000))
        self.solution_model = models.Polynomial1D(degree=3)
        self.solution_model.c0.value = 3500
        self.solution_model.c1.value = 1
        self.solution_model.c2.value = 1e-7

        self.non_linear_x_axis = self.solution_model(range(len(self.data)))

        self.feature_center = self.non_linear_x_axis[3000]

    def test_linearize_spectrum_nans_in_data(self):
        self.data[0:10] = np.nan
        self.assertRaises(SystemExit,
                          linearize_spectrum,
                          self.data,
                          self.solution_model)

    def test_linearize_spectrum_wrong_input(self):
        linear_data = linearize_spectrum(data=self.data,
                                         wavelength_solution=None)
        self.assertIsNone(linear_data)

    def test_linearize_spectrum(self):
        linear_x_axis, linear_data = linearize_spectrum(data=self.data,
                                         wavelength_solution=self.solution_model)

        new_gaussian = models.Gaussian1D(amplitude=100,
                                         mean=self.feature_center,
                                         stddev=5)
        fitter = fitting.LevMarLSQFitter()

        fitted_linear = fitter(new_gaussian, linear_x_axis, linear_data)

        np.testing.assert_array_almost_equal(self.feature_center,
                                             fitted_linear.mean.value,
                                             decimal=2)


class MasterFlatTest(TestCase):

    def setUp(self):
        # create a master flat
        self.master_flat = CCDData(data=np.ones((100, 100)),
                                   meta=fits.Header(),
                                   unit='adu')
        self.master_flat.header.set('GRATING', value='RALC_1200-BLUE')
        self.master_flat.header.set('SLIT', value='0.84" long slit')
        self.master_flat.header.set('FILTER2', value='<NO FILTER>')
        self.master_flat.header.set('WAVMODE', value='1200 m2')
        self.master_flat_name = 'master_flat_1200m2.fits'
        # expected master flat to be retrieved by get_best_flat
        self.reference_flat_name = 'master_flat_1200m2_0.84_dome.fits'
        # location of sample flats
        self.flat_path = os.path.join(os.path.dirname(__file__), '../../data/test_data/master_flat')
        slit = re.sub('[A-Za-z" ]',
                      '',
                      self.master_flat.header['SLIT'])
        self.flat_name_base = re.sub('.fits',
                                     '_' + slit + '*.fits',
                                     self.master_flat_name)

        # save a master flat with some random structure.

        self.master_flat_name_norm = 'flat_to_normalize.fits'
        # add a bias level
        self.master_flat.data += 300.
        # add noise
        self.master_flat.data += np.random.random_sample(
            self.master_flat.data.shape)

        self.master_flat.write(os.path.join(self.flat_path,
                                            self.master_flat_name_norm),
                               overwrite=False)

    def tearDown(self):
        full_path = os.path.join(self.flat_path,
                                 self.master_flat_name_norm)

        self.assertTrue(os.path.isfile(full_path))
        if os.path.isfile(full_path):
            os.unlink(full_path)
        self.assertFalse(os.path.isfile(full_path))

        # remove normalized flat
        norm_flat = re.sub('flat_to_', 'norm_flat_to_', full_path)
        if os.path.isfile(norm_flat):
            os.unlink(norm_flat)
        self.assertFalse(os.path.isfile(norm_flat))

    def test_get_best_flat(self):
        master_flat, master_flat_name = get_best_flat(
            flat_name=self.flat_name_base,
            path=self.flat_path)
        self.assertIsInstance(master_flat, CCDData)
        self.assertEqual(os.path.basename(master_flat_name),
                         self.reference_flat_name)

    def test_get_best_flat_fail(self):
        # Introduce an error that will never produce a result.
        wrong_flat_name = re.sub('1200m2', '1300m2', self.flat_name_base)
        master_flat, master_flat_name = get_best_flat(
            flat_name=wrong_flat_name,
            path=self.flat_path)
        self.assertIsNone(master_flat)
        self.assertIsNone(master_flat_name)

    def test_get_best_flat_relative_path(self):
        master_flat, master_flat_name = get_best_flat(
            flat_name=os.path.join(self.flat_path,self.flat_name_base),
            path=self.flat_path)
        self.assertIsInstance(master_flat, CCDData)
        self.assertEqual(os.path.basename(master_flat_name),
                         self.reference_flat_name)

    def test_normalize_master_flat(self):
        methods = ['mean', 'simple', 'full']
        for method in methods:
            self.assertNotAlmostEqual(self.master_flat.data.mean(), 1.)
            normalized_flat, normalized_flat_name = normalize_master_flat(
                master=self.master_flat,
                name=os.path.join(self.flat_path,
                                  self.master_flat_name_norm),
                method=method)

            self.assertAlmostEqual(normalized_flat.data.mean(), 1.,
                                   delta=0.001)
            self.assertEqual(normalized_flat.header['GSP_NORM'], method)
            self.assertIn('norm_', normalized_flat_name)


class NameMasterFlatsTest(TestCase):

    def setUp(self):
        self.reduced_data = os.getcwd()
        date = '2019-08-28'
        self.twilight_start_evening = '2019-08-27T23:45:00.022'
        self.twilight_end_morning = '2019-08-28T09:43:20.023'
        self.sun_set_time = '2019-08-27T22:21:00.437'
        self.sun_rise_time = '2019-08-28T11:07:11.851'

        self.header = fits.Header()
        self.header.set('GRATING', value='400_SYZY')
        self.header.set('GRT_TARG', value=7.5)
        self.header.set('CAM_TARG', value=16.1)
        self.header.set('FILTER', value='<NO FILTER>')
        self.header.set('FILTER2', value='<NO FILTER>')
        self.header.set('SLIT', value='1.0_LONG_SLIT')

    def test_name_master_flats_spectroscopy(self):

        expected_name = 'master_flat_TestSubject_400m2_GG455_1.0_dome.fits'

        self.header.set('FILTER2', value='GG455')
        self.header.set('DATE-OBS', value='2019-08-27T20:21:00.437')

        flat_name = name_master_flats(
            header=self.header,
            technique='Spectroscopy',
            reduced_data=self.reduced_data,
            sun_set=self.sun_set_time,
            sun_rise=self.sun_rise_time,
            evening_twilight=self.twilight_start_evening,
            morning_twilight=self.twilight_end_morning,
            target_name='TestSubject',
            get=False)
        self.assertEqual(expected_name, os.path.basename(flat_name))

    def test_name_master_flats_spectroscopy_no_grating_no_filter2(self):

        expected_name = 'master_flat_TestSubject_NO_GRATING_1.0_sky.fits'

        self.header.set('GRATING', value='<NO GRATING>')
        self.header.set('DATE-OBS', value='2019-08-27T23:40:00.022')

        flat_name = name_master_flats(
            header=self.header,
            technique='Spectroscopy',
            reduced_data=self.reduced_data,
            sun_set=self.sun_set_time,
            sun_rise=self.sun_rise_time,
            evening_twilight=self.twilight_start_evening,
            morning_twilight=self.twilight_end_morning,
            target_name='TestSubject',
            get=False)

        self.assertEqual(expected_name, os.path.basename(flat_name))

    def test_name_master_flats_imaging_no_filter(self):
        expected_name = 'master_flat_NO_FILTER_night.fits'
        self.header.set('DATE-OBS', value='2019-08-27T23:50:00.022')

        flat_name = name_master_flats(
            header=self.header,
            technique='Imaging',
            reduced_data=self.reduced_data,
            sun_set=self.sun_set_time,
            sun_rise=self.sun_rise_time,
            evening_twilight=self.twilight_start_evening,
            morning_twilight=self.twilight_end_morning,
            target_name='TestSubject',
            get=False)
        self.assertEqual(expected_name, os.path.basename(flat_name))

    def test_name_master_flats_imaging(self):
        expected_name = 'master_flat_u_BESSEL_night.fits'
        self.header.set('DATE-OBS', value='2019-08-27T23:50:00.022')
        self.header.set('FILTER', value='u-BESSEL')

        flat_name = name_master_flats(
            header=self.header,
            technique='Imaging',
            reduced_data=self.reduced_data,
            sun_set=self.sun_set_time,
            sun_rise=self.sun_rise_time,
            evening_twilight=self.twilight_start_evening,
            morning_twilight=self.twilight_end_morning,
            target_name='TestSubject',
            get=False)
        self.assertEqual(expected_name, os.path.basename(flat_name))

    def test_name_master_flats_get_true(self):
        expected_name = 'master_flat_TestSubject_400m2_GG455_1.0*.fits'

        self.header.set('FILTER2', value='GG455')
        self.header.set('DATE-OBS', value='2019-08-27T20:21:00.437')

        flat_name = name_master_flats(
            header=self.header,
            technique='Spectroscopy',
            reduced_data=self.reduced_data,
            sun_set=self.sun_set_time,
            sun_rise=self.sun_rise_time,
            evening_twilight=self.twilight_start_evening,
            morning_twilight=self.twilight_end_morning,
            target_name='TestSubject',
            get=True)
        self.assertEqual(expected_name, os.path.basename(flat_name))


class NightDataContainerTests(TestCase):

    def setUp(self):
        self.container = NightDataContainer(path=os.getcwd(),
                                            instrument='Red',
                                            technique='Spectroscopy')

        columns = ['file', 'obstype']
        sample_data_1 = [['file1.fits', 'OBJECT']]
        sample_data_2 = [['file1.fits', 'OBJECT'],
                         ['file2.fits', 'OBJECT']]

        self.sample_df_1 = pandas.DataFrame(sample_data_1, columns=columns)
        self.sample_df_2 = pandas.DataFrame(sample_data_2, columns=columns)

    def test___repr___method_empty(self):
        result = self.container.__repr__()
        self.assertEqual(result, 'Empty Data Container')

    def test___repr___method_not_empty(self):
        self.container.is_empty = False
        self.container.gain = 1
        self.container.rdnoise = 1
        self.container.roi = 'roi'
        result = self.container.__repr__()

        self.assertIn('Full Path: {:s}'.format(os.getcwd()), result)
        self.assertIn('Instrument: Red', result)
        self.assertIn('Technique: Spectroscopy', result)
        self.assertIn('Is Empty: False', result)

        _expected_content = ['Data Grouping Information',
                             'BIAS Group:',
                             'Group is Empty',
                             'Day FLATs Group:',
                             'Dome FLATs Group:',
                             'Sky FLATs Group:',
                             'COMP Group:',
                             'OBJECT Group',
                             'OBJECT + COMP Group:']

        for _line in _expected_content:
            self.assertIn(_line, result)

    def test___repr___method_imaging_not_empty(self):
        self.container.technique = 'Imaging'
        self.container.add_bias(bias_group=self.sample_df_2)
        self.container.add_day_flats(day_flats=self.sample_df_1)
        self.container.add_data_group(data_group=self.sample_df_2)

        #
        self.container.dome_flats = [self.sample_df_1]
        self.container.sky_flats = [self.sample_df_2]

        self.container.gain = 1
        self.container.rdnoise = 1
        self.container.roi = 'roi'
        result = self.container.__repr__()

        self.assertNotIn('Group is Empty', result)


    @skip
    def test__get_group_repr(self):
        pass

    def test_add_bias_imaging_insufficient_bias(self):
        self.container.technique = 'Imaging'
        self.container.add_bias(bias_group=self.sample_df_1)
        self.assertTrue(self.container.bias is None)
        self.assertTrue(self.container.is_empty)

    def test_add_bias_spectroscopy_insufficient_bias(self):
        self.container.add_bias(bias_group=self.sample_df_1)
        self.assertTrue(self.container.bias is None)
        self.assertTrue(self.container.is_empty)

    def test_add_bias(self):
        self.container.add_bias(bias_group=self.sample_df_2)
        self.container.add_bias(bias_group=self.sample_df_2)
        self.assertFalse(self.container.bias is None)
        self.assertFalse(self.container.is_empty)

    def test_add_day_flats(self):
        self.container.add_day_flats(day_flats=self.sample_df_1)
        self.assertIsInstance(self.container.day_flats[0], pandas.DataFrame)
        self.container.add_day_flats(day_flats=self.sample_df_2)
        self.assertFalse(self.container.day_flats is None)
        self.assertFalse(self.container.is_empty)

    def test_add_data_group(self):
        self.container.add_data_group(data_group=self.sample_df_1)
        self.assertIsInstance(self.container.data_groups[0], pandas.DataFrame)
        self.container.add_data_group(data_group=self.sample_df_2)
        self.assertFalse(self.container.data_groups is None)
        self.assertFalse(self.container.is_empty)

    def test_add_comp_group(self):
        self.container.add_comp_group(comp_group=self.sample_df_1)
        self.assertIsInstance(self.container.comp_groups[0], pandas.DataFrame)
        self.container.add_comp_group(comp_group=self.sample_df_2)
        self.assertFalse(self.container.comp_groups is None)
        self.assertFalse(self.container.is_empty)

    def test_add_object_group(self):
        self.container.add_object_group(object_group=self.sample_df_1)
        self.assertIsInstance(self.container.object_groups[0], pandas.DataFrame)
        self.container.add_object_group(object_group=self.sample_df_2)
        self.assertFalse(self.container.object_groups is None)
        self.assertFalse(self.container.is_empty)

    def test_add_spec_group(self):
        self.container.add_spec_group(spec_group=self.sample_df_1)
        self.assertIsInstance(self.container.spec_groups[0], pandas.DataFrame)
        self.container.add_spec_group(spec_group=self.sample_df_2)
        self.assertFalse(self.container.spec_groups is None)
        self.assertFalse(self.container.is_empty)

    def test_set_sun_times(self):
        _sun_set = '2019-01-01T18:00:00'
        _sun_rise = '2019-01-01T06:00:00'
        self.container.set_sun_times(sun_set=_sun_set, sun_rise=_sun_rise)

        self.assertEqual(self.container.sun_set_time, _sun_set)
        self.assertEqual(self.container.sun_rise_time, _sun_rise)

    def test_set_twilight_times(self):
        _evening = '2019-01-01T18:00:00'
        _morning = '2019-01-01T06:00:00'

        self.container.set_twilight_times(evening=_evening, morning=_morning)

        self.assertEqual(self.container.evening_twilight, _evening)
        self.assertEqual(self.container.morning_twilight, _morning)

    def test_set_readout(self):
        _gain = 1.48
        _rdnoise = 3.89
        _roi = 'Spectroscopic 2x2'

        self.container.set_readout(gain=_gain, rdnoise=_rdnoise, roi=_roi)

        self.assertEqual(self.container.gain, _gain)
        self.assertEqual(self.container.rdnoise, _rdnoise)
        self.assertEqual(self.container.roi, _roi)


class RaDecConversion(TestCase):

    def setUp(self):
        self.ra = '19:09:55.026'
        self.dec = '-68:18:01.901'
        self.reference_ra = 287.479275
        self.reference_dec = -68.3005281

    def test_ra_dec_to_deg_negative_dec(self):
        radeg, decdeg = ra_dec_to_deg(right_ascension=self.ra,
                                      declination=self.dec)
        self.assertAlmostEqual(radeg, self.reference_ra)
        self.assertAlmostEqual(decdeg, self.reference_dec)

    def test_ra_dec_to_deg_positive_dec(self):
        self.dec = '68:18:01.901'
        radeg, decdeg = ra_dec_to_deg(right_ascension=self.ra,
                                      declination=self.dec)
        self.assertAlmostEqual(radeg, self.reference_ra)
        self.assertAlmostEqual(decdeg, -1 * self.reference_dec)


class RecordTraceInformationTest(TestCase):

    def setUp(self):
        self.ccd = CCDData(data=np.ones((800, 2000)),
                           meta=fits.Header(),
                           unit='adu')

        self.all_keywords = ['GSP_TMOD',
                             'GSP_TORD',
                             'GSP_TC00',
                             'GSP_TC01',
                             'GSP_TC02',
                             'GSP_TERR']

        self.trace_info = collections.OrderedDict()

        self.trace_info['GSP_TMOD'] = ['Polinomial1D',
                                       'Model name used to fit trace']

        self.trace_info['GSP_TORD'] = [2, 'Degree of the model used to fit '
                                          'target trace']

        self.trace_info['GSP_TC00'] = [500, 'Parameter c0']
        self.trace_info['GSP_TC01'] = [1, 'Parameter c1']
        self.trace_info['GSP_TC02'] = [2, 'Parameter c2']
        self.trace_info['GSP_TERR'] = [0.5, 'RMS error of target trace']

    def test_record_trace_information(self):
        ccd = record_trace_information(ccd=self.ccd, trace_info=self.trace_info)
        new_keys = [key for key in ccd.header.keys()]

        self.assertTrue(all([key in new_keys for key in self.all_keywords]))
        self.assertEqual(ccd.header['GSP_TMOD'], 'Polinomial1D')
        self.assertEqual(ccd.header['GSP_TORD'], 2)


class ReferenceDataTest(TestCase):

    def setUp(self):
        self.rd = ReferenceData(
            reference_dir=os.path.join(os.path.dirname(__file__),
                                       '../../data/ref_comp'))
        self.ccd = CCDData(data=np.ones((800, 2000)),
                           meta=fits.Header(),
                           unit='adu')
        self.ccd.header.set('GRATING', value='400_SYZY')
        self.ccd.header.set('GRT_TARG', value=7.5)
        self.ccd.header.set('CAM_TARG', value=16.1)

        self.columns = ['object',
                        'grating',
                        'grt_targ',
                        'cam_targ',
                        'lamp_hga',
                        'lamp_ne',
                        'lamp_ar',
                        'lamp_fe',
                        'lamp_cu',]

        self.data_exist = [
            ['HgArNe',
             '400_SYZY',
             7.5,
             16.1,
             'TRUE',
             'TRUE',
             'FALSE',
             'FALSE',
             'FALSE'],
            ['HgAr',
             '400_SYZY',
             7.5,
             16.1,
             'TRUE',
             'FALSE',
             'FALSE',
             'FALSE',
             'FALSE']]

        self.data_does_not_exist = [
            ['HgArNe',
             'SYZY_800',
             7.5,
             16.1,
             'TRUE',
             'TRUE',
             'FALSE',
             'FALSE',
             'FALSE'],
            ['HgAr',
             'SYZY_800',
             7.5,
             16.1,
             'TRUE',
             'FALSE',
             'FALSE',
             'FALSE',
             'FALSE']]

    def test_get_reference_lamp_exist_with_lamps_status_key(self):
        self.ccd.header.set('LAMP_HGA', value='TRUE')
        self.ccd.header.set('LAMP_NE', value='TRUE')
        self.ccd.header.set('LAMP_AR', value='FALSE')
        self.ccd.header.set('LAMP_FE', value='FALSE')
        self.ccd.header.set('LAMP_CU', value='FALSE')
        self.ccd.header.set('LAMP_QUA', value='FALSE')
        self.ccd.header.set('LAMP_QPE', value=0)
        self.ccd.header.set('LAMP_BUL', value='FALSE')
        self.ccd.header.set('LAMP_DOM', value='FALSE')
        self.ccd.header.set('LAMP_DPE', value=0)


        self.ccd.header.set('WAVMODE', value='400_M2')

        ref_lamp = self.rd.get_reference_lamp(header=self.ccd.header)

        self.assertIsInstance(ref_lamp, CCDData)
        self.assertEqual(ref_lamp.header['LAMP_HGA'], self.ccd.header['LAMP_HGA'])
        self.assertEqual(ref_lamp.header['LAMP_NE'], self.ccd.header['LAMP_NE'])
        self.assertEqual(ref_lamp.header['WAVMODE'], self.ccd.header['WAVMODE'])

    def test_get_reference_lamp_no_match_status_keys(self):
        self.ccd.header.set('LAMP_HGA', value='TRUE')
        self.ccd.header.set('LAMP_NE', value='TRUE')
        self.ccd.header.set('LAMP_AR', value='TRUE')
        self.ccd.header.set('LAMP_FE', value='TRUE')
        self.ccd.header.set('LAMP_CU', value='TRUE')
        self.ccd.header.set('LAMP_QUA', value='TRUE')
        self.ccd.header.set('LAMP_QPE', value=0)
        self.ccd.header.set('LAMP_BUL', value='FALSE')
        self.ccd.header.set('LAMP_DOM', value='FALSE')
        self.ccd.header.set('LAMP_DPE', value=0)

        self.ccd.header.set('WAVMODE', value='400_M2')

        self.assertRaises(NoMatchFound,
                          self.rd.get_reference_lamp,
                          self.ccd.header)

    def test_get_reference_lamp_exist_with_object_key(self):
        self.ccd.header.set('OBJECT', value='HgArNe')
        self.ccd.header.set('WAVMODE', value='400_M2')

        self.ccd.header.set('GSP_P', value='1')
        self.ccd.header.set('GSP_P', value='2')
        self.ccd.header.set('GSP_P', value='3')
        self.ccd.header.set('GSP_P', value='4')

        self.ccd.header.set('GSP_A', value='0')
        self.ccd.header.set('GSP_A', value='100')
        self.ccd.header.set('GSP_A', value='0')
        self.ccd.header.set('GSP_A', value='200')

        ref_lamp = self.rd.get_reference_lamp(header=self.ccd.header)

        self.assertIsInstance(ref_lamp, CCDData)
        self.assertEqual(ref_lamp.header['OBJECT'], self.ccd.header['OBJECT'])
        self.assertEqual(ref_lamp.header['WAVMODE'], self.ccd.header['WAVMODE'])

    def test_get_reference_lamp_does_not_exist(self):
        self.ccd.header.set('OBJECT', value='HgArCu')
        self.ccd.header.set('WAVMODE', value='400_M5')

        self.assertRaises(NoMatchFound,
                          self.rd.get_reference_lamp,
                          self.ccd.header)

    def test_lamp_exist(self):
        self.ccd.header.set('LAMP_HGA', value='TRUE')
        self.ccd.header.set('LAMP_NE', value='TRUE')
        self.ccd.header.set('LAMP_AR', value='FALSE')
        self.ccd.header.set('LAMP_FE', value='FALSE')
        self.ccd.header.set('LAMP_CU', value='FALSE')
        self.ccd.header.set('LAMP_QUA', value='FALSE')
        self.ccd.header.set('LAMP_QPE', value=0)
        self.ccd.header.set('LAMP_BUL', value='FALSE')
        self.ccd.header.set('LAMP_DOM', value='FALSE')
        self.ccd.header.set('LAMP_DPE', value=0)
        self.ccd.header.set('WAVMODE', value='400_M2')
        self.assertTrue(self.rd.lamp_exists(header=self.ccd.header))

        # HgArNeCu is impossible
        self.ccd.header.set('LAMP_CU', value='TRUE')
        self.assertFalse(self.rd.lamp_exists(header=self.ccd.header))

    def test_check_comp_group__lamp_exists(self):

        comp_group = pandas.DataFrame(self.data_exist,
                                      columns=self.columns)

        new_group = self.rd.check_comp_group(comp_group=comp_group)

        self.assertIsInstance(new_group, pandas.DataFrame)
        self.assertFalse(comp_group.equals(new_group))
        self.assertEqual(len(new_group), 1)

    def test_check_comp_group__lamp_does_not_exist(self):
        comp_group = pandas.DataFrame(self.data_does_not_exist,
                                      columns=self.columns)

        new_group = self.rd.check_comp_group(comp_group=comp_group)

        self.assertIsInstance(new_group, pandas.DataFrame)
        self.assertTrue(comp_group.equals(new_group))

    def test__order_validation(self):
        should_be_true = self.rd._order_validation(range(10))

        self.assertTrue(should_be_true)

        should_be_false = self.rd._order_validation(range(10)[::-1])
        self.assertFalse(should_be_false)

    def test__load_nist_list(self):
        self.assertIsInstance(self.rd.nist, dict)
        self.assertEqual(0, len(self.rd.nist))

        self.rd._load_nist_list()
        self.assertIsInstance(self.rd.nist, dict)
        self.assertGreater(len(self.rd.nist), 0)


class SaturationValuesTest(TestCase):

    def setUp(self):
        self.ccd = CCDData(data=np.ones((100, 100)),
                           meta=fits.Header(),
                           unit='adu')
        self.ccd.header.set('INSTCONF', value='Red')
        self.ccd.header.set('GAIN', value=1.48)
        self.ccd.header.set('RDNOISE', value=3.89)

        self.half_full_well = 69257

        self.saturation_values = SaturationValues(ccd=self.ccd)

    def test_half_full_well_value(self):
        self.assertEqual(self.saturation_values.saturation_value,
                         self.half_full_well)

    def test_empty_result(self):
        self.ccd.header['GAIN'] = 2.3
        result = self.saturation_values.get_saturation_value(ccd=self.ccd)
        self.assertIsNone(result)
        self.assertIsNone(self.saturation_values.saturation_value)


class SearchCompGroupTest(TestCase):

    def setUp(self):
        columns = ['object',
                   'grating',
                   'cam_targ',
                   'grt_targ',
                   'filter',
                   'filter2',
                   'lamp_hga',
                   'lamp_ne',
                   'lamp_ar',
                   'lamp_fe',
                   'lamp_cu']

        self.object_group = pandas.DataFrame(
            data=[['NGC2070',
                   'SYZY_400',
                   16.1,
                   7.5,
                   '<NO FILTER>',
                   'GG455',
                   'TRUE',
                   'FALSE',
                   'FALSE',
                   'FALSE',
                   'FALSE'
                   ]],
            columns=columns)
        self.object_group_no_match = pandas.DataFrame(
            data=[['NGC2070',
                   'SYZY_600',
                   16.1,
                   7.5,
                   '<NO FILTER>',
                   'GG455',
                   'TRUE',
                   'FALSE',
                   'FALSE',
                   'FALSE',
                    'FALSE']],
            columns=columns)
        self.comp_groups = [
            pandas.DataFrame(
                data=[['HgArNe',
                       'SYZY_400',
                       16.1,
                       7.5,
                       '<NO FILTER>',
                       'GG455',
                       'TRUE',
                       'FALSE',
                       'FALSE',
                       'FALSE',
                       'FALSE']], columns=columns),
            pandas.DataFrame(
                data=[['CuArNe',
                       'SYZY_400',
                       11.6,
                       5.8,
                       '<NO FILTER>',
                       'GG455',
                       'TRUE',
                       'FALSE',
                       'FALSE',
                       'FALSE',
                       'FALSE']], columns=columns)]

        self.reference_data = ReferenceData(
            reference_dir=os.path.join(os.path.dirname(__file__),
                                       '../../data/ref_comp'))

    def test_search_comp_group(self):
        result = search_comp_group(
            object_group=self.object_group,
            comp_groups=self.comp_groups,
            reference_data=self.reference_data)
        self.assertIsInstance(result, pandas.DataFrame)
        self.assertFalse(result.empty)

    def test_search_comp_group_no_match(self):
        with self.assertRaises(NoMatchFound):
            search_comp_group(
                object_group=self.object_group_no_match,
                comp_groups=self.comp_groups,
                reference_data=self.reference_data)


class SlitTrimTest(TestCase):
    # TODO (simon): discuss with Bruno

    def setUp(self):
        # Create fake image
        self.fake_image = CCDData(data=np.ones((100, 100)),
                                  meta=fits.Header(),
                                  unit='adu')

        # define
        self.slit_low_limit = 5
        self.slit_high_limit = 95

        self.reference_slit_trim = '[1:100,{:d}:{:d}]'.format(
            self.slit_low_limit + 10 + 1,
            self.slit_high_limit - 10)

        # make a flat-like structure
        self.fake_image.data[self.slit_low_limit:self.slit_high_limit, :] = 100

    def test_get_slit_trim_section__slit_within_data(self):

        slit_trim = get_slit_trim_section(master_flat=self.fake_image)
        self.assertEqual(slit_trim, self.reference_slit_trim)

    def test_get_slit_trim_section__slit_full_data(self):
        self.fake_image.data[:, :] = 100

        slit_trim = get_slit_trim_section(master_flat=self.fake_image)
        self.assertEqual(slit_trim, '[1:100,1:100]')

    def test_image_trim_slit(self):
        self.fake_image = image_trim(ccd=self.fake_image,
                                     trim_section=self.reference_slit_trim,
                                     trim_type='slit')
        self.assertIsInstance(self.fake_image, CCDData)
        reference_size = (self.slit_high_limit - 10) - \
                         (self.slit_low_limit + 10)
        self.assertEqual(self.fake_image.data.shape, (reference_size, 100))

        self.assertEqual(self.fake_image.header['GSP_SLIT'],
                         self.reference_slit_trim)


class SpectroscopicModeTest(TestCase):

    def setUp(self):
        self.sm = SpectroscopicMode()
        self.ccd = CCDData(data=np.ones((800, 2000)),
                           meta=fits.Header(),
                           unit='adu')
        self.ccd.header.set('GRATING', value='SYZY_400')
        self.ccd.header.set('CAM_TARG', value='16.1')
        self.ccd.header.set('GRT_TARG', value='7.5')
        self.ccd.header.set('FILTER2', value='GG455')

    def test__call__(self):
        self.assertRaises(SyntaxError, self.sm)

        mode_m2_header = self.sm(header=self.ccd.header)

        self.assertEqual(mode_m2_header, 'm2')

        mode_m2_keywords = self.sm(grating=self.ccd.header['GRATING'],
                                   camera_targ=self.ccd.header['CAM_TARG'],
                                   grating_targ=self.ccd.header['GRT_TARG'],
                                   blocking_filter=self.ccd.header['FILTER2'])

        self.assertEqual(mode_m2_keywords, 'm2')

    def test_get_modes(self):
        mode_m2 = self.sm.get_mode(grating='400',
                                   camera_targ='16.1',
                                   grating_targ='7.5',
                                   blocking_filter='GG455')
        self.assertEqual(mode_m2, 'm2')

        mode_custom_400 = self.sm.get_mode(grating='400',
                                           camera_targ='16.1',
                                           grating_targ='6.6',
                                           blocking_filter='GG455')

        self.assertEqual(mode_custom_400, 'Custom_7000nm')

        mode_custom_2100 = self.sm.get_mode(grating='2100',
                                            camera_targ='16.1',
                                            grating_targ='7.5',
                                            blocking_filter='GG455')
        self.assertEqual(mode_custom_2100, 'Custom_1334nm')

    def test_get_mode_400_m1_from_header(self):
        header = fits.Header()
        header.set('GRATING', value='400_SYZY')
        header.set('CAM_TARG', value='11.6')
        header.set('GRT_TARG', value='5.8')
        header.set('FILTER2', value='NO_FILTER')

        spectroscopic_mode = self.sm(header=header)

        self.assertEqual(spectroscopic_mode, 'm1')

    def test_get_mode_400_m1_from_old_header(self):
        header = fits.Header()
        header.set('GRATING', value='SYZY_400')
        header.set('CAM_TARG', value='11.6')
        header.set('GRT_TARG', value='5.8')
        header.set('FILTER2', value='<NO FILTER>')

        spectroscopic_mode = self.sm(header=header)

        self.assertEqual(spectroscopic_mode, 'm1')

    def test_get_mode_1200_m3_from_header(self):
        header = fits.Header()
        header.set('GRATING', value='1200')
        header.set('CAM_TARG', value='39.4')
        header.set('GRT_TARG', value='20.2')
        header.set('FILTER2', value='NO_FILTER')

        spectroscopic_mode = self.sm(header=header)

        self.assertEqual(spectroscopic_mode, 'm3')

    def test_get_cam_grt_targ_angle(self):

        cam_targ, grt_targ = self.sm.get_cam_grt_targ_angle(1800, 'm10')
        self.assertIsNone(cam_targ)
        self.assertIsNone(grt_targ)

        cam_targ, grt_targ = self.sm.get_cam_grt_targ_angle(930, 'm5')
        self.assertEqual(cam_targ, '39.4')
        self.assertEqual(grt_targ, '19.7')

        cam_targ, grt_targ = self.sm.get_cam_grt_targ_angle(930, 'm7')
        self.assertIsNone(cam_targ)
        self.assertIsNone(grt_targ)


class TargetsTest(TestCase):

    def setUp(self):
        self.ccd = CCDData(data=np.ones((800, 2000)),
                           meta=fits.Header(),
                           unit='adu')

        self.ccd.header.set('GSP_FNAM',
                            value='fake-name.fits',
                            comment='Fake file name')

        self.profile_1 = models.Gaussian1D(amplitude=200,
                                           mean=100,
                                           stddev=10).rename('Profile_1')
        self.profile_2 = models.Gaussian1D(amplitude=200,
                                           mean=600,
                                           stddev=10).rename('Profile_2')

        self.profile_3 = models.Moffat1D(amplitude=200,
                                         x_0=600,
                                         gamma=3).rename('Profile_3')

        profile_sum = self.profile_1 + self.profile_2
        self.ccd2 = self.ccd.copy()
        self.no_targets_ccd = self.ccd.copy()
        for i in range(self.ccd.data.shape[1]):
            self.ccd.data[:, i] *= profile_sum(range(self.ccd.data.shape[0]))
            self.ccd2.data[:, i] *= self.profile_3(
                range(self.ccd2.data.shape[0]))

    def tearDown(self):
        del self.ccd
        del self.profile_1
        del self.profile_2
        del self.profile_3

    def test_identify_targets_moffat(self):
        self.ccd.header.set('OBSTYPE',
                            value='OBJECT',
                            comment='Fake values')
        self.ccd.header.set('SLIT',
                            value='1.03" long slit',
                            comment='Fake slit')
        self.ccd.header.set('CCDSUM',
                            value='1 1',
                            comment='Fake values')
        targets = identify_targets(ccd=self.ccd,
                                   fit_model='moffat',
                                   background_threshold=3,
                                   nfind=2,
                                   plots=False)
        self.assertEqual(len(targets), 2)
        for target in targets:
            self.assertIsInstance(target, Model)

    def test_identify_targets_gaussian(self):
        self.ccd.header.set('OBSTYPE',
                            value='OBJECT',
                            comment='Fake values')
        self.ccd.header.set('SLIT',
                            value='1.03" long slit',
                            comment='Fake slit')
        self.ccd.header.set('CCDSUM',
                            value='1 1',
                            comment='Fake values')
        targets = identify_targets(ccd=self.ccd,
                                   fit_model='gaussian',
                                   background_threshold=3,
                                   nfind=2,
                                   plots=False)
        self.assertEqual(len(targets), 2)
        for target in targets:
            self.assertIsInstance(target, Model)

    def test_identify_targets_empty_output(self):
        self.no_targets_ccd.header.set('OBSTYPE',
                                       value='OBJECT',
                                       comment='Fake values')
        self.no_targets_ccd.header.set('SLIT',
                                       value='1.03" long slit',
                                       comment='Fake slit')
        self.no_targets_ccd.header.set('CCDSUM',
                                       value='1 1',
                                       comment='Fake values')
        targets = identify_targets(ccd=self.no_targets_ccd,
                                   fit_model='gaussian',
                                   background_threshold=3,
                                   nfind=2,
                                   plots=False)
        self.assertEqual(len(targets), 0)

    def test_trace_gaussian(self):
        trace_model = models.Polynomial1D(degree=2)
        fitter = fitting.LevMarLSQFitter()
        test_trace, trace_rms = trace(ccd=self.ccd,
                                      model=self.profile_1,
                                      trace_model=trace_model,
                                      model_fitter=fitter,
                                      sampling_step=5)
        self.assertEqual(test_trace.c0.value, self.profile_1.mean.value)
        self.assertAlmostEqual(test_trace.c1.value, 0.)
        self.assertAlmostEqual(test_trace.c2.value, 0.)

    def test_trace_moffat(self):
        trace_model = models.Polynomial1D(degree=2)
        fitter = fitting.LevMarLSQFitter()
        test_trace, trace_rms = trace(ccd=self.ccd2,
                                      model=self.profile_3,
                                      trace_model=trace_model,
                                      model_fitter=fitter,
                                      sampling_step=5)
        self.assertEqual(test_trace.c0.value, self.profile_3.x_0.value)
        self.assertAlmostEqual(test_trace.c1.value, 0.)
        self.assertAlmostEqual(test_trace.c2.value, 0.)

    def test_trace_not_implemented(self):
        trace_model = models.Polynomial1D(degree=2)
        fitter = fitting.LevMarLSQFitter()
        self.assertRaises(NotImplementedError,
                          trace,
                          self.ccd2,
                          models.BlackBody(),
                          trace_model,
                          fitter,
                          5)

    def test_trace_targets(self):
        targets = [self.profile_1, self.profile_2]
        all_traces = trace_targets(ccd=self.ccd,
                                   target_list=targets,
                                   sampling_step=5,
                                   pol_deg=2,
                                   nfwhm=2,
                                   plots=False)
        for new_trace, profile, trace_info in all_traces:
            self.assertEqual(new_trace.c0.value, profile.mean.value)
            self.assertAlmostEqual(new_trace.c1.value, 0)
            self.assertAlmostEqual(new_trace.c2.value, 0)


class TimeConversionTest(TestCase):

    def setUp(self):
        self.test_time_str = '2018-01-17T12:05:44.250'
        self.test_time_sec = 1516190744.0

    def test_convert_time(self):
        self.assertEqual(convert_time(self.test_time_str), self.test_time_sec)

    def test_get_twilight_time(self):
        expected_evening_twilight = '2018-01-17T01:21:26.120'
        expected_morning_twilight = '2018-01-17T08:24:38.922'
        expected_sun_set_time = '2018-01-17T23:43:46.775'
        expected_sun_rise_time = '2018-01-17T10:02:04.510'
        evening_twilight, morning_twilight, sun_set, sun_rise\
            = get_twilight_time([self.test_time_str])

        self.assertEqual(expected_evening_twilight, evening_twilight)
        self.assertEqual(expected_morning_twilight, morning_twilight)
        self.assertEqual(expected_sun_set_time, sun_set)
        self.assertEqual(expected_sun_rise_time, sun_rise)


class ValidateCcdRegionTest(TestCase):

    def test_validate_ccd_region_valid(self):
        self.assertTrue(validate_ccd_region, '[1:1,10:10]')

    def test_validate_ccd_region_invalid(self):
        self.assertRaises(SyntaxError, validate_ccd_region, "10:10:10]")
