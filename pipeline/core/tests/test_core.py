from __future__ import absolute_import

from unittest import TestCase, skip
from ccdproc import CCDData
from astropy.io import fits
from astropy.modeling import Model
import numpy as np
import os

# import all classes in core.py
from ..core import (NightDataContainer,
                    NoMatchFound,
                    NotEnoughLinesDetected,
                    NoTargetException,
                    SpectroscopicMode)

from astropy.modeling import (models,
                              fitting)

# import of functions in core.py
from ..core import (add_wcs_keys,
                    call_cosmic_rejection,
                    classify_spectroscopic_data,
                    convert_time,
                    dcr_cosmicray_rejection,
                    extraction,
                    extract_fractional_pixel,
                    extract_optimal,
                    fractional_sum,
                    get_best_flat,
                    get_central_wavelength,
                    get_slit_trim_section,
                    get_twilight_time,
                    identify_targets,
                    image_overscan,
                    image_trim,
                    lacosmic_cosmicray_rejection,
                    normalize_master_flat,
                    ra_dec_to_deg,
                    read_fits,
                    search_comp_group,
                    trace,
                    trace_targets,
                    write_fits)

# class ExceptionHandling(TestCase):
#
#     def test_critical_error(self):
#         self.assertRaises(CriticalError)
#
#
#     def test_night_data_container(self):
#         pass
#
#
#     def test_no_match_found(self):
#         pass
#
#
#     def test_not_enough_lines_detected(self):
#         pass
#
#
#     def test_no_target_exception(self):
#         pass


def test_spectroscopic_mode():
    pass


class AddWCSKeywordsTest(TestCase):

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

        test_ccd = CCDData(data=np.ones((100, 100)),
                           meta=fits.Header(),
                           unit='adu')

        test_ccd = add_wcs_keys(ccd=test_ccd)
        for key in wcs_keys:
            self.assertIn(key, test_ccd.header)


class CosmicRayRejectionTest(TestCase):

    @skip
    def test_dcr_cosmicray_rejection(self):
        pass

    @skip
    def test_call_cosmic_rejection(self):
        pass


def test_classify_spectroscopic_data():
    pass


class TimeConversionTest(TestCase):

    def setUp(self):
        self.test_time_str = '2018-01-17T12:05:44.250'
        self.test_time_sec = 1516190744.0

    def test_convert_time(self):
        self.assertEqual(convert_time(self.test_time_str), self.test_time_sec)

    def test_get_twilight_time(self):
        expected_evening_twilight = '2018-01-17T01:21:26.113'
        expected_morning_twilight = '2018-01-17T08:24:38.919'
        expected_sun_set_time = '2018-01-17T23:43:46.782'
        expected_sun_rise_time = '2018-01-17T10:02:04.508'
        evening_twilight, morning_twilight, sun_set, sun_rise\
            = get_twilight_time([self.test_time_str])

        self.assertEqual(evening_twilight, expected_evening_twilight)
        self.assertEqual(morning_twilight, expected_morning_twilight)
        self.assertEqual(sun_set, expected_sun_set_time)
        self.assertEqual(sun_rise, expected_sun_rise_time)


class ExtractionTest(TestCase):

    def setUp(self):
        self.fake_image = CCDData(data=np.ones((100, 100)),
                             meta=fits.Header(),
                             unit='adu')
        self.fake_image.header['OBSTYPE'] = 'COMP'
        self.fake_image.header['GSP_FNAM'] = 'fake-image.fits'

        # Create model aligned with pixels - represents the trace
        self.target_trace = models.Linear1D(slope=0, intercept=50.3)

        # Calculate the STDDEV
        self.stddev = 8.4

        # Calculate how many STDDEV will be extracted - N_STDDEV
        self.n_stddev = 2

        # Calculate how far the background is from the the center.
        self.distance = 1

        self.target_profile = models.Gaussian1D(amplitude=1,
                                                mean=50.3,
                                                stddev=self.stddev)

        self.reference_result = np.ones(100) * self.stddev * self.n_stddev


    def test_fractional_extraction(self):
        # Perform extraction
        extracted_array, background = extract_fractional_pixel(
            ccd=self.fake_image,
            target_trace=self.target_trace,
            target_stddev=self.stddev,
            extraction_width=self.n_stddev,
            background_spacing=self.distance)
        # assert isinstance(fake_image, CCDData)
        self.assertIsInstance(extracted_array, CCDData)

        np.testing.assert_array_almost_equal(extracted_array,
                                             self.reference_result)

    def test_fractional_sum(self):

        fake_image = np.ones((100, 100))
        low_limit = 50 + np.random.random()
        high_limit = 60 + np.random.random()

        sum = fractional_sum(fake_image, 50, low_limit, high_limit)
        self.assertEqual(sum, high_limit - low_limit)

    def test_extract_optimal(self):
        self.assertRaises(NotImplementedError, extract_optimal)

    def test_extraction(self):
        extracted = extraction(ccd=self.fake_image,
                               target_trace=self.target_trace,
                               spatial_profile=self.target_profile,
                               extraction_name='fractional')
        self.assertIsInstance(extracted, CCDData)
        np.testing.assert_array_almost_equal(extracted, self.reference_result)

    def test_extraction_exception(self):
        self.assertRaises(NotImplementedError, extraction, ccd=self.fake_image,
                          target_trace=self.target_trace,
                          spatial_profile=self.target_profile,
                          extraction_name='optimal')


# class BackgroundValue(TestCase):


def test_get_best_flat():
    pass


def test_get_central_wavelength():
    pass


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

    def test_get_slit_trim_section(self):

        slit_trim = get_slit_trim_section(master_flat=self.fake_image)
        # print(fake_image.data[:,5])
        # print(slit_trim)
        self.assertEqual(slit_trim, self.reference_slit_trim)

    def test_image_trim_slit(self):
        # # define
        # slit_low_limit = 5
        # slit_high_limit = 95
        #
        # slit_trim = '[1:100,{:d}:{:d}]'.format(slit_low_limit + 10 + 1,
        #                                        slit_high_limit - 10)
        self.fake_image = image_trim(ccd=self.fake_image,
                                     trim_section=self.reference_slit_trim,
                                     trim_type='slit')
        self.assertIsInstance(self.fake_image, CCDData)
        reference_size = (self.slit_high_limit - 10) - (self.slit_low_limit + 10)
        self.assertEqual(self.fake_image.data.shape, (reference_size, 100))

        self.assertEqual(self.fake_image.header['GSP_SLIT'],
                         self.reference_slit_trim)





def test_lacosmic_cosmicray_rejection():
    pass


def test_normalize_master_flat():
    pass


def test_ra_dec_to_deg():
    pass


def test_search_comp_group():
    pass


class TargetsTest(TestCase):

    def setUp(self):
        self.ccd = CCDData(data=np.ones((800, 2000)),
                           meta=fits.Header(),
                           unit='adu')

        self.profile_1 = models.Gaussian1D(amplitude=200,
                                           mean=100,
                                           stddev=10).rename('Profile_1')
        self.profile_2 = models.Gaussian1D(amplitude=200,
                                           mean=600,
                                           stddev=10).rename('Profile_2')

        profile_sum = self.profile_1 + self.profile_2
        for i in range(self.ccd.data.shape[1]):
            self.ccd.data[:, i] *= profile_sum(range(self.ccd.data.shape[0]))

    def tearDown(self):
        del self.ccd
        del self.profile_1
        del self.profile_2

    def test_identify_targets(self):
        self.ccd.header.set('OBSTYPE',
                            value='OBJECT',
                            comment='Fake values')
        self.ccd.header.set('SLIT',
                            value='1.03" long slit',
                            comment='Fake slit')
        self.ccd.header.set('CCDSUM',
                            value='1 1',
                            comment='Fake values')
        targets = identify_targets(ccd=self.ccd, nfind=2, plots=False)
        self.assertEqual(len(targets), 2)
        for target in targets:
            self.assertIsInstance(target, Model)

    def test_trace(self):
        trace_model = models.Polynomial1D(degree=2)
        fitter = fitting.LevMarLSQFitter()
        test_trace = trace(ccd=self.ccd,
                           model=self.profile_1,
                           trace_model=trace_model,
                           model_fitter=fitter,
                           sampling_step=5)
        self.assertEqual(test_trace.c0.value, self.profile_1.mean.value)
        self.assertAlmostEqual(test_trace.c1.value, 0.)
        self.assertAlmostEqual(test_trace.c2.value, 0.)

    def test_trace_targets(self):
        targets = [self.profile_1, self.profile_2]
        all_traces = trace_targets(ccd=self.ccd,
                                   target_list=targets,
                                   sampling_step=5,
                                   pol_deg=2,
                                   nsigmas=2,
                                   plots=False)
        for new_trace, profile in all_traces:
            self.assertEqual(new_trace.c0.value, profile.mean.value)
            self.assertAlmostEqual(new_trace.c1.value, 0)
            self.assertAlmostEqual(new_trace.c2.value, 0)


class FitsFileIOAndOps(TestCase):

    def setUp(self):
        self.fake_image = CCDData(data=np.ones((100, 100)),
                                  meta=fits.Header(),
                                  unit='adu')
        self.fake_image.header.set('CCDSUM',
                                   value='1 1',
                                   comment='Fake values')

        self.file_name = 'sample_file.fits'
        self.current_directory = os.getcwd()
        self.full_path = os.path.join(self.current_directory, self.file_name)
        self.parent_file = 'parent_file.fits'
        self.fake_image.write(self.full_path, overwrite=False)

    def test_write_fits(self):
        self.assertTrue(os.path.isfile(self.full_path))
        os.remove(self.full_path)
        write_fits(ccd=self.fake_image,
                   full_path=self.full_path,
                   parent_file=self.parent_file,
                   overwrite=False)
        self.assertTrue(os.path.isfile(self.full_path))

    def test_read_fits(self):
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

    def test_image_trim(self):
        self.assertEqual(self.fake_image.data.shape, (100, 100))
        trim_section = '[1:50,:]'
        self.fake_image = image_trim(ccd=self.fake_image,
                                     trim_section=trim_section,
                                     trim_type='trimsec')

        self.assertEqual(self.fake_image.data.shape, (100, 50))
        self.assertEqual(self.fake_image.header['GSP_TRIM'], trim_section)

    def tearDown(self):
        self.assertTrue(os.path.isfile(self.full_path))
        os.remove(self.full_path)


