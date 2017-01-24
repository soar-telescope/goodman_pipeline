# -*- coding: utf8 -*-

"""
# PyGoodman CCD Reduction - CCD reductions for Goodman spectroscopic data.

This script performs CCD reduction for long-slit spectra taken with the
Goodman High Throughput Spectrograph at SOAR Telescope. The script will
make (in order):

- BIAS subtraction
- TRIMMING including slit edges identification (it does not work for MOS spectroscopy)
- FLAT correction
- COSMIC rays rejection (optional)

Users can add a flag in order to clean up science data from cosmic rays, which
are removed by using the LACosmic code (P. G. van Dokkum, 2001, PASP, 113, 1420).

## I/O Data Structure

This script was designed to make CCD reduction for any spectrograph configuration,
but the input directory should contain only an unique spectral configuration (binning,
grating, slit, gain, rdnoise, CCD ROI, etc). The input dir should contain only the
following frames:

- BIAS frames
- FLAT frames  (Flats taken between science exposures will be trimmed and bias subtracted.)
- ARC frames   (data from focus sequence will not be reduced)
- SCIENCE and/or STANDARD frames

Please, inspect you calibration and throw it out the bad ones. The output data has the same
filename of the input data, but with a prefix "fzh". It means data has its header updated (h),
bias subtracted (z) and flat corrected (f). The prefix "c_fzh" means that cosmic ray correction
was applied.

## How to use it...

It can be be executed in terminal running:

    $ python goodman_ccdreduction.py [options] raw_path red_path

More information abotu the options and how to use it can be otained by using...

    $ python goodman_ccdreduction.py --help
or
    $ python goodman_ccdreduction.py -h

Documentation for specific functions of the code can be found directly in the corresponding
function. (To be done...)

#ToDo

- Consider internal illumination correction (in fact this will not be done)
- Disable the flat correction if there is a no grating flat
- Automatically determine the best input parameters for LACOSMIC

David Sanmartim (dsanmartim at gemini.edu)
August 2016

Thanks to Bruno Quint for all comments and helping.

"""

import sys
import os
import glob
import argparse
import numpy as np
from astropy.io import fits
from astropy import log
from astropy import units as u
from scipy.interpolate import interp1d

from astropy.coordinates import EarthLocation
from astropy.time import Time, TimeDelta
from astroplan import Observer
from astroplan import get_IERS_A_or_workaround, download_IERS_A

import ccdproc
from ccdproc import ImageFileCollection
from ccdproc import CCDData
import warnings

__author__ = 'David Sanmartim'
__date__ = '2016-07-15'
__version__ = "0.1"
__email__ = "dsanmartim@ctio.noao.edu"
__maintainer__ = "Simon Torres"


class Main:
    def __init__(self):

        self.args = self.get_args()

        # Soar Geodetic Location and other definitions
        self.observatory = 'SOAR Telescope'
        self.Geodetic_Location = ['-70d44m01.11s', '-30d14m16.41s', 2748]
        self.longitude = self.Geodetic_Location[0]
        self.latitude = self.Geodetic_Location[1]
        self.elevation = self.Geodetic_Location[2]
        self.timezone = 'UTC'
        self.description = 'Soar Telescope on Cerro Pachon, Chile'

        # Set variables used globally
        self.master_bias = None
        self.slit1 = None
        self.slit2 = None
        self.master_flat = None
        self.master_flat_nogrt = None
        self.master_flat_name = None
        self.master_flat_nogrt_name = None

        # ToDo Check if the file already exist before download it
        # if get_IERS_A_or_workaround() is None:
        #     download_IERS_A(show_progress=True)

        # Memory Limit to be used
        self.memlim = 32E9
        # self.memlim = 1E7
        ## For computer with up to 8GB of RAM
        # self.memlim = 6E6

        # Taking some args from argparse method
        self.raw_path = str(os.path.join(self.args.raw_path[0], ''))
        self.red_path = str(os.path.join(self.args.red_path[0], ''))

        if self.raw_path == self.red_path:
            raise ValueError('raw_path may not be equal to red_path')
        else:
            pass

        # More args from argparse
        self.clean = self.args.clean
        self.slit = self.args.slit

        # checking the reduction directory
        if not os.path.isdir(self.red_path):
            os.mkdir(self.red_path)
        os.chdir(self.red_path)

        # About warnings
        warnings.filterwarnings('ignore')
        log.propagate = False

    def __call__(self, *args, **kwargs):

        # cleaning up the reduction dir
        self.clean_path(self.red_path)

        # Fixing header and shape of raw data
        self.fix_header_and_shape(self.raw_path, self.red_path, prefix='h_', overwrite=True)

        # Create image file collection for raw data
        ic = ImageFileCollection(self.red_path)

        # Getting twilight time
        twi_eve, twi_mor = self.get_twilight_time(ic, self.observatory, self.longitude, self.latitude,
                                                  self.elevation, self.timezone, self.description)
        # Create master_flats
        self.create_daymaster_flat(ic, twi_eve, twi_mor, self.slit, self.memlim)

        # Create master bias
        if len(ic.files_filtered(obstype='BIAS')) > 0:
            self.create_master_bias(ic, self.slit, self.memlim)
        else:
            log.info('No BIAS image detected')
            log.warning('The images will be processed but the results will not be optimal')

        # Reduce Night Flat frames (if they exist)
        self.reduce_nightflats(ic, twi_eve, twi_mor, self.slit, prefix='z')

        # Reduce Arc frames
        self.reduce_arc(ic, self.slit, prefix='fz')

        # Reduce Sci frames
        self.reduce_sci(ic, self.slit, self.clean, prefix='fz')

        return

    def get_args(self):
        # Parsing Arguments ---
        parser = argparse.ArgumentParser(description="PyGoodman CCD Reduction - CCD reductions for "
                                                     "Goodman spectroscopic data")

        parser.add_argument('-c', '--clean', action='store_true',
                            help="Clean cosmic rays from science data.")

        parser.add_argument('-s', '--slit', action='store_true',
                            help="Find slit edge to make an additional trimming (recommended).")

        parser.add_argument('raw_path', metavar='raw_path', type=str, nargs=1,
                            help="Full path to raw data (e.g. /home/jamesbond/soardata/).")

        parser.add_argument('red_path', metavar='red_path', type=str, nargs=1,
                            help="Full path to reduced data (e.g /home/jamesbond/soardata/RED/).")

        # parser.add_argument('--red-camera', action='store_true', default=False, dest='red_camera',
        #                    help='Enables Goodman Red Camera')

        args = parser.parse_args()
        return args

    @staticmethod
    def fit_spline3(y, order=3, nsum=5):
        """Fit a cubib spline to an 1D-array of N pixels.

        Args:
            y (1D array like): A 1-D array of monotonically increasing real values
            order (int) : order of t
            nsum (ins)  : number of array elements to be avareged

        Returns: It returns a function

        Examples:

        f = fit_spline3(y, order=5, nsum=2)
        x = np.arange(0,1500,1)
        ysmooth = f(x)

        """
        y_resampled = [np.median(y[i:i + nsum]) for i in range(0, len(y) - len(y) % nsum, nsum)]
        x_resampled = np.linspace(0, len(y), len(y_resampled))

        # Fitting
        f = interp1d(x_resampled, y_resampled, kind=order, bounds_error=True)

        # Return function to be constructed with any other x array
        return f

    # Local Minima and Maxima
    @staticmethod
    def local_minmax(data, nmin=1, nmax=1):
        """Find local minima-maxima points for a set of non-noisy data

        Args:
            data (1D array like): 1D array of non-noisy data
            nmin (int): number of local minina to be find
            nmax (int): number of local maxima to be find

        Returns: It returns a function

        """

        # Identifying indices of local minima-maxima points
        id_min = (np.gradient(np.sign(np.gradient(data))) > 0).nonzero()[0]  # index of local min
        id_max = (np.gradient(np.sign(np.gradient(data))) < 0).nonzero()[0]  # index of local max

        # Taking values at min/max points
        list_min, list_max = data[id_min], data[id_max]

        # Sorting minima-maxima values (bigger --> lower)
        list_min, id_min = (list(p) for p in zip(*sorted(zip(list_min, id_min), reverse=False)))
        list_max, id_max = (list(p) for p in zip(*sorted(zip(list_max, id_max), reverse=True)))

        # Taking the desired number of local minima-maxima points
        list_min, list_max, id_min, id_max = list_min[0:nmin], list_max[0:nmax], id_min[0:nmin], id_max[0:nmax]

        return list_min, list_max, id_min, id_max

    @staticmethod
    def clean_path(path):
        """
        Clean up FIST file in a directoy. It's not recursive
        """
        if os.path.exists(path):
            for _file in glob.glob(os.path.join(path, '*.fits')):
                os.remove(_file)

    @staticmethod
    def fix_header_and_shape(input_path, output_path, prefix, overwrite=False):

        """Remove/Update some  inconvenient parameters in the header of the Goodman FITS
        files. Some of these parameters contain non-printable ASCII characters. The ouptut
        files are created in the output_path. Also convert fits from 3D [1,X,Y] to 2D [X,Y].

        Args:
            input_path (str): Location of input data.
            output_path (str): Location of output data.
            prefix (str): Prefix to be added in the filename of output data
            overwrite (bool, optional): If true it it overwrite existing data

        Returns:
            Fits file with header and shape fixed.
        """

        for _file in sorted(glob.glob(os.path.join(input_path, '*.fits'))):

            ccddata, hdr = fits.getdata(_file, header=True, ignore_missing_end=True)

            # if not self.args.red_camera:

            # 3D to 2D
            if ccddata.ndim is 3:
                ccddata = ccddata[0]
                hdr['NAXIS'] = 2

            # keywords to remove
            key_list_to_remove = ['PARAM0', 'PARAM61', 'PARAM62', 'PARAM63', 'NAXIS3', 'INSTRUME']

            # Keyword to be changed (3 --> 2)
            try:
                hdr['N_PARAM'] -= len(key_list_to_remove)
                # Specific keywords to be removed
                for key in key_list_to_remove:
                    if (key in hdr) is True:
                        hdr.remove(keyword=key)
            except KeyError as key_error:
                log.debug(key_error)

            # Removing duplicated keywords
            key_list = []
            for key in hdr.iterkeys():
                if key in key_list:
                    hdr.remove(keyword=key)
                key_list.append(key)

            hdr.add_history('Header and Shape fixed.')
            fits.writeto(os.path.join(output_path, '') + prefix + os.path.basename(_file), ccddata, hdr,
                         clobber=overwrite)
            log.info('Header of ' + os.path.basename(_file) + ' has been updated --> ' + prefix
                     + os.path.basename(_file))
        log.info('Done: All headers have been updated.')
        print('\n')
        return

    def find_slitedge(self, ccddata):
        """Find slit edge by inspecting signal variation in the spatial direction
        of flat frames. The spatial direction is assumed to be axis=0 (or y axis
        in IRAF convention) and data are divided in two regions in which we are
        looking for slit edges.

        Args:
            ccddata (ccdproc.CCDData): The actual data contained in this ccdproc.CCDData object

        Returns (int):
            Pixel of slit edge 1 and slit edge 2 (bottom to top of the flat image).
        """
        # Reading and Collapsing flat in the dispersion direction
        flat_collapsed = np.sum(ccddata, axis=1) / ccddata.shape[1]

        # Excluding 3 first pixels in the spatial direction
        cut = 3
        c_flat = flat_collapsed[cut:-cut]
        c_lines = np.arange(0, c_flat.size, 1)

        # Fitting cubic spline. It's working very well with order=5, nsum=2
        func_splin3 = self.fit_spline3(c_flat, order=5, nsum=2)
        smooth_flat = func_splin3(c_lines)

        # Compute 1st and flat smoothed
        dy2 = np.gradient(np.gradient(smooth_flat))

        # Region one: it represent first 40 percent of all data. Region two: ... last 40%
        pixa, pixb = int(len(c_flat) * 0.4), int(len(c_flat) * 0.6)
        dy2_one, dy2_two = dy2[0:pixa], dy2[pixb:]

        # Reg. 1: Compute local min/max of the 2nd derivative
        list_min_1, _, id_min_1, _ = self.local_minmax(dy2_one, nmin=1, nmax=1)
        list_min_2, _, id_min_2, _ = self.local_minmax(dy2_two, nmin=1, nmax=1)

        # Slit edges are the local maxima/minima 1/2 [accounting the cutted pixels]
        slit_1, slit_2 = int(np.array(id_min_1) + cut), int(np.array(np.array(id_min_2) + pixb) + cut)

        return slit_1, slit_2

    @staticmethod
    def get_twilight_time(image_collection, observatory, longitude, latitude, elevation, timezone, description):
        """

        Args:
            image_collection:
            observatory:
            longitude:
            latitude:
            elevation:
            timezone:
            description:

        Returns:

        Old...

        image_collection: ccdproc object
        observatory: str, observatory name (e.g. 'Soar Telescope',
        long: str, dms or deg
        lat: str, dms or deg
        elevation: int, meters (define through ellipsoid WGS84)
        timezone: str, eg. 'UTC'
        description: str, short description of the observatory

        return: str, twilight evening and twilinght morning (format 'YYYY-MM-DDT00:00:00.00')
        """
        soar_loc = EarthLocation.from_geodetic(longitude, latitude, elevation * u.m, ellipsoid='WGS84')

        soar = Observer(name=observatory, location=soar_loc, timezone=timezone, description=description)

        dateobs_list = image_collection.values('date-obs')
        time_first_frame, time_last_frame = Time(min(dateobs_list)), Time(max(dateobs_list))

        twilight_evening = soar.twilight_evening_astronomical(Time(time_first_frame), which='nearest').isot
        twilight_morning = soar.twilight_morning_astronomical(Time(time_last_frame), which='nearest').isot

        return twilight_evening, twilight_morning

    @staticmethod
    def get_night_flats(image_collection, twilight_evening, twilight_morning):

        """

        Args:
            image_collection:
            twilight_evening:
            twilight_morning:

        Returns:

        """

        df = image_collection.summary.to_pandas()

        start_night = (Time(twilight_evening) - TimeDelta(1800.0, format='sec')).isot
        end_night = (Time(twilight_morning) + TimeDelta(1800.0, format='sec')).isot

        night_condition = ((Time(df['date-obs'].tolist()).jd < Time(start_night).jd) &
                           (Time(df['date-obs'].tolist()).jd > Time(end_night).jd))

        dfobj = df['file'][(df['obstype'] == 'FLAT') & night_condition]
        night_flat_list = dfobj.values.tolist()

        return night_flat_list

    @staticmethod
    def get_day_flats(image_collection, twilight_evening, twilight_morning):
        """
        image_collection: ccdproc object
        return: list of flats
        """
        df = image_collection.summary.to_pandas()

        start_night = (Time(twilight_evening) - TimeDelta(1800.0, format='sec')).isot
        end_night = (Time(twilight_morning) + TimeDelta(1800.0, format='sec')).isot

        day_condition = ((Time(df['date-obs'].tolist()).jd < Time(start_night).jd) |
                         (Time(df['date-obs'].tolist()).jd > Time(end_night).jd))

        dfobj = df['file'][(df['obstype'] == 'FLAT') & day_condition]
        dayflat_list = dfobj.values.tolist()

        return dayflat_list

    def create_daymaster_flat(self, image_collection, twilight_evening, twilight_morning, slit, memory_limit):
        """

        Args:
            image_collection:
            twilight_evening:
            twilight_morning:
            slit:
            memory_limit:

        Returns:

        """

        self.master_flat = {}
        self.master_flat_nogrt = {}

        # Creating dict. of flats. The key values are expected to be: GRATIN_ID and '<NO GRATING>'
        # if there is flat taken w/o grating
        df = image_collection.summary.to_pandas()
        grtobj = df['grating'][(df['obstype'] != 'BIAS')]
        grtobj = grtobj.unique()
        grt_list = grtobj.tolist()

        dic_all_flats = {}
        for grt in sorted(grt_list):
            start_night = (Time(twilight_evening) - TimeDelta(1800.0, format='sec')).isot
            end_night = (Time(twilight_morning) + TimeDelta(1800.0, format='sec')).isot

            day_condition = ((Time(df['date-obs'].tolist()).jd < Time(start_night).jd) |
                             (Time(df['date-obs'].tolist()).jd > Time(end_night).jd))

            dfobj = df['file'][(df['obstype'] == 'FLAT') & (df['grating'] == grt) & day_condition]
            dic_all_flats[str(grt)] = dfobj.tolist()

        # Dict. for flats with grating and without grating
        dic_flat = {grt: dic_all_flats[grt] for grt in dic_all_flats if grt != "<NO GRATING>"}
        dic_flatnogrt = {grt: dic_all_flats[grt] for grt in dic_all_flats if grt == "<NO GRATING>"}

        if np.size(dic_flat.values()) > 0:

            for grt in dic_flat.keys():

                flat_list = []
                log.info('Combining and trimming flat frames:')
                for filename in dic_flat[grt]:
                    log.info(filename)
                    ccd = CCDData.read(os.path.join(image_collection.location, '') + filename, unit=u.adu)
                    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])
                    flat_list.append(ccd)

                # combinning and trimming slit edges
                print('Flat list length: %s' % len(flat_list))
                if len(flat_list) >= 1:
                    master_flat = ccdproc.combine(flat_list, method='median', mem_limit=memory_limit,
                                                       sigma_clip=True,
                                                       sigma_clip_low_thresh=1.0, sigma_clip_high_thresh=1.0)
                    # self.master_flat.append(master_flat)
                else:
                    log.info('Flat list empty')
                    return
                if slit is True:
                    print('\n Finding slit edges... \n')
                    self.slit1, self.slit2 = self.find_slitedge(master_flat)
                    master_flat = ccdproc.trim_image(self.master_flat[self.slit1:self.slit2, :])
                # self.master_flat.append(master_flat)

                self.master_flat_name = self.get_flat_name(master_flat.header,get_name_only=True)
                self.master_flat[self.master_flat_name] = master_flat
                master_flat.write(self.master_flat_name, clobber=True)

                log.info('Done: master flat has been created --> ' + self.master_flat_name)
                print('\n')

        else:
            log.info('Flat files have not been found.')
            print('\n')

        if np.size(dic_flatnogrt.values()) > 0:
            # No grating flats
            for grt in dic_flatnogrt.keys():
                no_grating_ic = df[['file','filter', 'filter2']][(df['grating'] == grt) & (df['obstype'] == 'FLAT')]
                no_grating_filters_1 = no_grating_ic['filter'].unique()
                no_grating_filters_2 = no_grating_ic['filter2'].unique()
                for filter_1 in no_grating_filters_1:
                    for filter_2 in no_grating_filters_2:
                        # print("filter ", filter_1, filter_2)
                        no_grating_files = no_grating_ic.file[((df['filter'] == filter_1) & (df['filter2'] == filter_2))]
                        # print(no_grating_files)
                        flatnogrt_list = []
                        log.info('Combining and trimming flat frame taken without grating:')
                        for filename in no_grating_files:
                            log.info(filename)
                            ccd = CCDData.read(os.path.join(image_collection.location, '') + filename, unit=u.adu)
                            ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])

                            flatnogrt_list.append(ccd)

                        # combining and trimming slit edges
                        master_flat_nogrt = ccdproc.combine(flatnogrt_list, method='median',
                                                                 mem_limit=memory_limit,
                                                                 sigma_clip=True,
                                                                 sigma_clip_low_thresh=3.0, sigma_clip_high_thresh=3.0)

                        if slit is True:
                            master_flat_nogrt = ccdproc.trim_image(
                                self.master_flat_nogrt[self.slit1:self.slit2, :])

                        # self.master_flat.append(master_flat_nogrt)
                        # filter_string = ''
                        # if filter_1 != '<NO FILTER>':
                        #     filter_string += '_' + filter_1
                        # if filter_2 != '<NO FILTER>':
                        #     filter_string += '_' + filter_2

                        self.master_flat_nogrt_name = self.get_flat_name(master_flat_nogrt.header, get_name_only=True)
                        self.master_flat[self.master_flat_nogrt_name] = master_flat_nogrt
                        # print self.master_flat_nogrt_name
                        master_flat_nogrt.write(self.master_flat_nogrt_name, clobber=True)

                        log.info(
                            'Done: master flat have been created --> ' + self.master_flat_nogrt_name)
                        print('\n')


        else:
            log.info('Flat files taken without grating not found or not necessary')
            print('\n')

        return

    def create_master_bias(self, image_collection, slit, memory_limit):

        bias_list = []
        log.info('Combining and trimming bias frames:')
        for filename in image_collection.files_filtered(obstype='BIAS'):
            log.info(filename)
            ccd = CCDData.read(os.path.join(image_collection.location, '') + filename, unit=u.adu)
            # Finding overscan regions... getting from header and assuming it is at the right edge...
            # over_start = int((ccd.header['TRIMSEC'].split(':'))[1].split(',')[0]) - 1
            # over_start += 10 / int(ccd.header['CCDSUM'][0])
            # ccd = ccdproc.subtract_overscan(ccd, median=True, overscan_axis=1, overscan=ccd[:, over_start:])
            ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])

            bias_list.append(ccd)

        self.master_bias = ccdproc.combine(bias_list, method='median', mem_limit=memory_limit, sigma_clip=True,
                                           sigma_clip_low_thresh=3.0, sigma_clip_high_thresh=3.0)
        if slit is True:
            self.master_bias = ccdproc.trim_image(self.master_bias[self.slit1:self.slit2, :])
            # else:
            # self.master_bias = self.master_bias
        self.master_bias.header['HISTORY'] = "Trimmed."
        self.master_bias.write('master_bias.fits', clobber=True)

        # Now I obtained bias... subtracting bias from master flat
        # Testing if master_flats are not empty arrays
        if (not self.master_flat) is False:
            for master_flat_name in self.master_flat.keys():
                master_flat = self.master_flat[master_flat_name]
                fccd = ccdproc.subtract_bias(master_flat, self.master_bias)
                fccd.header['HISTORY'] = "Trimmed. Bias subtracted. Flat corrected."
                fccd.write(master_flat_name, clobber=True)

        if (not self.master_flat_nogrt) is False:
            ngccd = ccdproc.subtract_bias(self.master_flat_nogrt, self.master_bias)
            ngccd.header['HISTORY'] = "Trimmed. Bias subtracted. Flat corrected."
            ngccd.write(self.master_flat_nogrt_name, clobber=True)

        log.info('Done: a master bias have been created --> master_bias.fits')
        print('\n')
        return

    def reduce_nightflats(self, image_collection, twilight_evening, twilight_morning, slit, prefix):

        log.info('Reducing flat frames taken during the night...')

        df = image_collection.summary.to_pandas()

        # Night starts/ends 30min beforre/after twilight evening/morning
        start_night = (Time(twilight_evening) - TimeDelta(1800.0, format='sec')).isot
        end_night = (Time(twilight_morning) + TimeDelta(1800.0, format='sec')).isot

        night_condition = ((Time(df['date-obs'].tolist()).jd > Time(start_night).jd) &
                           (Time(df['date-obs'].tolist()).jd < Time(end_night).jd))

        dfobj = df['file'][(df['obstype'] == 'FLAT') & (df['grating'] != '<NO GRATING>') & night_condition]
        nightflat_list = dfobj.tolist()

        if len(nightflat_list) > 0:
            for filename in sorted(nightflat_list):
                log.info('Trimming and bias subtracting frame ' + filename + ' --> ' + prefix + filename)
                ccd = CCDData.read(os.path.join(image_collection.location, '') + filename, unit=u.adu)
                ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])
                ccd.header['HISTORY'] = "Trimmed"
                if slit is True:
                    ccd = ccdproc.trim_image(ccd[self.slit1:self.slit2, :])
                if self.master_bias is not None:
                    ccd = ccdproc.subtract_bias(ccd, self.master_bias)
                    ccd.header['HISTORY'] = "Bias subtracted."
                else:
                    ccd.header['HISTORY'] = "Bias NOT subtracted."
                    log.warning('No bias subtraction!')
                ccd.write(prefix + filename, clobber=True)
            log.info('Done --> Night flat frames have been reduced.')
            print('\n')
        return

    def reduce_arc(self, image_collection, slit, prefix):

        log.info('Reducing Arc frames...')

        arc_list = image_collection.files_filtered(obstype='COMP')

        if len(arc_list) > 0:
            for filename in sorted(arc_list):
                log.info('Reducing Arc frame ' + filename + ' --> ' + prefix + filename)
                ccd = CCDData.read(os.path.join(image_collection.location, '') + filename, unit=u.adu)
                ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])
                if slit is True:
                    ccd = ccdproc.trim_image(ccd[self.slit1:self.slit2, :])
                if self.master_bias is not None:
                    ccd = ccdproc.subtract_bias(ccd, self.master_bias)
                    ccd.header['HISTORY'] = "Bias subtracted."
                else:
                    ccd.header['HISTORY'] = "Bias NOT subtracted."
                    log.warning('No bias subtraction!')
                flat_name = self.get_flat_name(ccd.header)
                if flat_name != False:
                    ccd = ccdproc.flat_correct(ccd, self.master_flat[flat_name])
                    ccd.header['HISTORY'] = "Trimmed. Flat corrected."
                    ccd.write(prefix + filename, clobber=True)
                else:
                    log.info('No flat found to process ' + filename)
            log.info('Done --> Arc frames have been reduced.')
            print('\n')
        return

    def reduce_sci(self, image_collection, slit, clean, prefix):

        log.info('Reducing Sci/Std frames...')
        for filename in image_collection.files_filtered(obstype='OBJECT'):
            log.info('Reducing Sci/Std frame ' + filename + ' --> ' + prefix + filename)
            ccd = CCDData.read(os.path.join(image_collection.location, '') + filename, unit=u.adu)
            ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])
            if slit is True:
                ccd = ccdproc.trim_image(ccd[self.slit1:self.slit2, :])
            if self.master_bias is not None:
                try:
                    ccd = ccdproc.subtract_bias(ccd, self.master_bias)
                    ccd.header['HISTORY'] = "Bias subtracted."
                except ValueError as err:
                    log.error("Data must be of only one kind. Please check your source data.")
                    continue
            else:
                ccd.header['HISTORY'] = "Bias NOT subtracted."
                log.warning('No bias subtraction!')
            # print self.master_flat.header['filter'], self.master_flat.header['grating'], ccd.header['filter'], ccd.header['grating']
            print ccd.header['filter'], ccd.header['grating']

            flat_name = self.get_flat_name(ccd.header)
            if flat_name != False:
                # print flat_name, ccd.header['GRATING']
                ccd = ccdproc.flat_correct(ccd, self.master_flat[flat_name])
                # OBS: cosmic ray rejection is working pretty well by defining gain = 1. It's not working
                # when we use the real gain of the image. In this case the sky level changes by a factor
                # equal the gain.
                # Function to determine the sigfrac and objlim: y = 0.16 * exptime + 1.2
                value = 0.16 * float(ccd.header['EXPTIME']) + 1.2
                if clean is True:
                    log.info('Cleaning cosmic rays... ')
                    nccd, _ = ccdproc.cosmicray_lacosmic(ccd.data, sigclip=2.5, sigfrac=value, objlim=value,
                                                         gain=float(ccd.header['GAIN']),
                                                         readnoise=float(ccd.header['RDNOISE']),
                                                         satlevel=np.inf, sepmed=True, fsmode='median',
                                                         psfmodel='gaussy', verbose=True)
                    log.info('Cosmic rays have been cleaned ' + prefix + filename + ' --> ' + 'c' + prefix + filename)
                    print('\n')
                    nccd = np.array(nccd, dtype=np.double) / float(ccd.header['GAIN'])
                    ccd.header['HISTORY'] = "Trimmed. Flat corrected."
                    ccd.header['HISTORY'] = "Cosmic rays rejected."
                    fits.writeto('c' + prefix + filename, nccd, ccd.header, clobber=True)
                elif clean is False:
                    ccd.header['HISTORY'] = "Trimmed, Flat corrected."
                ccd.write(prefix + filename, clobber=True)
            else:
                log.info('No flat found to process ' + filename)
        log.info('Done: Sci/Std frames have been reduced.')
        print('\n')
        return

    def get_flat_name(self, header, get_name_only=False):

        name_text = ''
        if header['grating'] == '<NO GRATING>':
            name_text += '_nogrt'
        else:
            grating = header['grating'].split('_')[1]
            name_text += '_' + grating
        if header['filter'] != '<NO FILTER>':
            name_text += '_' + header['filter']
        if header['filter2'] != '<NO FILTER>':
            name_text += '_' + header['filter2']

        flat_name = 'master_flat' + name_text + '.fits'

        if flat_name in self.master_flat.keys() or get_name_only:
            log.info('Master flat Name: ' + flat_name)
            return flat_name
        else:
            log.error('There is no flat suitable for use')
            return False



        # def run(self):

        # cleaning up the reduction dir
        # self.clean_path(self.red_path)

        # Fixing header and shape of raw data
        # self.fix_header_and_shape(self.raw_path, self.red_path, prefix='h.', overwrite=True)

        # Create image file collection for raw data
        # ic = ImageFileCollection(self.red_path)

        # Getting twilight time
        # twi_eve, twi_mor = self.get_twilight_time(ic, self.observatory, self.longitude, self.latitude,
        #                                           self.elevation, self.timezone, self.description)
        # Create master_flats
        # self.create_daymaster_flat(ic, twi_eve, twi_mor, self.slit, self.memlim)

        # Create master bias
        # if len(ic.files_filtered(obstype='BIAS')) > 0:
        # self.create_master_bias(ic, self.slit, self.memlim)
        # else:
        # log.info('No bias detected')

        # Reduce Night Flat frames (if they exist)
        # self.reduce_nightflats(ic, twi_eve, twi_mor, self.slit, prefix='z')

        # Reduce Arc frames
        # self.reduce_arc(ic, self.slit, prefix='fz')

        # Reduce Sci frames
        # self.reduce_sci(ic, self.slit, self.clean, prefix='fz')

        # return


if __name__ == '__main__':

    # # Parsing Arguments ---
    # parser = argparse.ArgumentParser(description="PyGoodman CCD Reduction - CCD reductions for "
    #                                              "Goodman spectroscopic data")
    #
    # parser.add_argument('-c', '--clean', action='store_true',
    #                     help="Clean cosmic rays from science data.")
    #
    # parser.add_argument('-s', '--slit', action='store_true',
    #                     help="Find slit edge to make an additional trimming (recommended).")
    #
    # parser.add_argument('raw_path', metavar='raw_path', type=str, nargs=1,
    #                     help="Full path to raw data (e.g. /home/jamesbond/soardata/).")
    #
    # parser.add_argument('red_path', metavar='red_path', type=str, nargs=1,
    #                     help="Full path to reduced data (e.g /home/jamesbond/soardata/RED/).")
    #
    # # parser.add_argument('--red-camera', action='store_true', default=False, dest='red_camera',
    # #                    help='Enables Goodman Red Camera')
    #
    # args = parser.parse_args()

    main = Main()
    main()

# else:
#    print('goodman_ccdreduction.py is not being executed as main.')
