from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
from ccdproc import ImageFileCollection
import matplotlib.pyplot as plt
import time
import sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import re
from astropy.coordinates import EarthLocation
from astropy.time import Time, TimeDelta
from astroplan import Observer
from astropy import units as u
import logging
from .core import convert_time, get_twilight_time, ra_dec_to_deg

log = logging.getLogger('goodmanccd.nightorganizer')


class NightOrganizer(object):

    def __init__(self, full_path, instrument, technique, ignore_bias=False):
        """Initializes the NightOrganizer class

        This class contains methods to organize the data for processing. It will identify groups of OBJECTS, FLATS or
        COMPS (comparison lamps) whenever they exist. The product will be an object that will act as a data container.

        Args:
            args (object): Argparse object. Contains all the runtime arguments.
            night_dict (dict): A dictionary that contains full path, instrument and observational technique.

        """
        self.path = full_path
        self.instrument = instrument
        self.technique = technique
        self.ignore_bias = ignore_bias
        self.keywords = ['date',
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
        self.file_collection = None
        self.all_datatypes = None
        self.data_container = Night(path=self.path, instrument=self.instrument, technique=self.technique)
        self.day_time_data = None
        self.night_time_data = None

    def __call__(self):
        """Call method

        Creates a table with selected keywords that will allow to group the data in order to be classified according to
        the observational technique used, imaging or spectroscopy.

        Returns:
            data_container (object): Class used as storage unit for classified data.

        """

        ifc = ImageFileCollection(self.path, self.keywords)
        self.file_collection = ifc.summary.to_pandas()
        # add two columns that will contain the ra and dec in degrees
        # TODO (simon): This part creates a warning originated from Pandas. Fixit
        self.file_collection['radeg'] = ''
        self.file_collection['decdeg'] = ''
        for i in self.file_collection.index.tolist():
            radeg, decdeg = ra_dec_to_deg(self.file_collection.obsra.iloc[i], self.file_collection.obsdec.iloc[i])
            self.file_collection.radeg.iloc[i] = '{:.2f}'.format(radeg)
            self.file_collection.decdeg.iloc[i] = '{:.2f}'.format(decdeg)
            # now we can compare using degrees
        self.initial_checks()
        self.all_datatypes = self.file_collection.obstype.unique()
        if self.technique == 'Spectroscopy':
            self.day_time_data, self.night_time_data = self.separate_day_night()
            # print(len(self.day_time_data), len(self.night_time_data))
            if len(self.day_time_data) != 0:
                self.spectroscopy_day_time()
            else:
                log.warning('There is no day time data!')
                if not self.ignore_bias:
                    log.error('BIAS are needed for optimal results')
                    log.info('Check the argument --ignore-bias')
                    sys.exit('BIAS needed')
            if len(self.night_time_data) != 0:
                self.spectroscopy_night_time()
            else:
                log.warning('There is no night time data!')
                # return False
        elif self.technique == 'Imaging':
            self.imaging_night()

        return self.data_container

    def initial_checks(self):
        readout_confs = self.file_collection.groupby(['gain', 'rdnoise'])
        if len(readout_confs) > 1:
            log.warning('There are %s different readout modes in the data.', len(readout_confs))
            log.info('Sleeping 10 seconds')
            time.sleep(10)

    def spectroscopy_day_time(self):
        """Organizes day time calibration data for spectroscopy.

        This methods assumes that during the day only calibration data was taken, therefore it will search for
        bias and flat images only. Then each group will be divided in subgroups that correspond to different
        configurations using information in the header.
        Bias images are divided according to four keywords. GAIN, RDNOISE, OBSRA and OBSDEC. The last two define the
        intended position of the telescope. Although BIAS images are not affected by telescope position the separation
        is made for grouping purposes only, since a different set of bias will be most likely obtained at a different
        telescope position.
        Flat images are filtered by seven parameters, in header keywords words they are: GAIN, RDNOISE, GRATING,
        FILTER2, CAM_TARG, GRT_TARG and SLIT. For most of the keywords is easy to understand why they where considered,
        I would like to note though that only FILTER2 was considered because for spectroscopy, order blocking filters
        are located only on the second filter wheel, and the first filter wheel is dedicated for imaging filters.
        CAM_TARG and GRT_TARG are the Camera and Grating angles, those values are fixed and define the different that
        every grating can operate and that's why is a good classifier parameter.

        Notes:
            Bias or Flat will be only considered if there are three or more images.

        Notes:
            All bias images will be considered part of daytime data.

        """
        if len(self.day_time_data) != 0:
            # add bias
            bias_data = self.day_time_data[self.day_time_data.obstype == 'BIAS']

            if len(bias_data) > 2:
                bias_confs = bias_data.groupby(['gain',
                                                'rdnoise',
                                                'radeg',
                                                'decdeg']).size().reset_index().rename(columns={0: 'count'})
                for i in bias_confs.index:
                    bias_group = bias_data[((bias_data['gain'] == bias_confs.iloc[i]['gain']) &
                                            (bias_data['rdnoise'] == bias_confs.iloc[i]['rdnoise']) &
                                            (bias_data['radeg'] == bias_confs.iloc[i]['radeg']) &
                                            (bias_data['decdeg'] == bias_confs.iloc[i]['decdeg']))]
                    self.data_container.add_bias(bias_group)
            else:
                log.error('Not enough bias images.')
            # add flats
            flat_data = self.day_time_data[self.day_time_data.obstype == 'FLAT']
            if len(flat_data) > 2:
                # configurations
                confs = flat_data.groupby(['gain',
                                           'rdnoise',
                                           'grating',
                                           'filter2',
                                           'cam_targ',
                                           'grt_targ',
                                           'slit']).size().reset_index().rename(columns={0: 'count'})
                for i in confs.index:
                    # print(confs.iloc[i]['grating'])
                    flat_group = flat_data[((flat_data['gain'] == confs.iloc[i]['gain']) &
                                            (flat_data['rdnoise'] == confs.iloc[i]['rdnoise']) &
                                            (flat_data['grating'] == confs.iloc[i]['grating']) &
                                           (flat_data['filter2'] == confs.iloc[i]['filter2']) &
                                           (flat_data['cam_targ'] == confs.iloc[i]['cam_targ']) &
                                            (flat_data['grt_targ'] == confs.iloc[i]['grt_targ']) &
                                            (flat_data['slit'] == confs.iloc[i]['slit']))]
                    # print(flat_group.file)
                    self.data_container.add_day_flats(flat_group)
            else:
                log.error('Not enough flat images.')
            # if there are object data discard them
            # print(self.day_time_data)
            pass
        else:
            log.warning('There is no day time data.')

    def spectroscopy_night_time(self):
        """Organizes night time data for spectroscopy

        This method identifies all combinations of nine **key** keywords that can set appart different objects with
        their respective calibration data or not. The keywords used are: GAIN, RDNOISE, GRATING, FILTER2, CAM_TARG,
        GRT_TARG, SLIT, OBSRA and OBSDEC.

        This method populates the `data_container` class attribute which is an instance of the class Night.
        A data group is an instance of a Pandas DataFrame.

        """

        # confs stands for configurations
        confs = self.night_time_data.groupby(['gain',
                                              'rdnoise',
                                              'grating',
                                              'filter2',
                                              'cam_targ',
                                              'grt_targ',
                                              'slit',
                                              'radeg',
                                              'decdeg']).size().reset_index().rename(columns={0: 'count'})
        for i in confs.index:
            night_time_group = self.night_time_data[((self.night_time_data['gain'] == confs.iloc[i]['gain']) &
                                                     (self.night_time_data['rdnoise'] == confs.iloc[i]['rdnoise']) &
                                                     (self.night_time_data['grating'] == confs.iloc[i]['grating']) &
                                                     (self.night_time_data['filter2'] == confs.iloc[i]['filter2']) &
                                                     (self.night_time_data['cam_targ'] == confs.iloc[i]['cam_targ']) &
                                                     (self.night_time_data['grt_targ'] == confs.iloc[i]['grt_targ']) &
                                                     (self.night_time_data['slit'] == confs.iloc[i]['slit']) &
                                                     (self.night_time_data['radeg'] == confs.iloc[i]['radeg']) &
                                                     (self.night_time_data['decdeg'] == confs.iloc[i]['decdeg']))]
            self.data_container.add_data_group(night_time_group)
            # sss = night_time_group.groupby(['obstype']).size().reset_index().rename(columns={0: 'count'})
            # if 'OBJECT' not in sss.obstype.tolist() or len(sss) < 3:
            #     log.warning('Less than three obstype')
            #     # for i in sss.index:
            #     print(sss)
            #     print(' ')
            #     print(night_time_group)
            #     print(' ')
            #     time.sleep(3)

    def imaging_night(self):
        """Organizes data for imaging

        For imaging there is no discrimination regarding night data since the process is simpler. It is a three
        stage process classifying BIAS, FLAT and OBJECT datatype. The data is packed in groups that are pandas.DataFrame
        objects.

        """

        # bias data group
        date_obs_list = self.file_collection['date-obs'].tolist()
        afternoon_twilight, morning_twilight, sun_set, sun_rise = get_twilight_time(date_obs=date_obs_list)
        self.data_container.set_sun_times(sun_set, sun_rise)
        self.data_container.set_twilight_times(afternoon_twilight, morning_twilight)
        bias_group = self.file_collection[self.file_collection.obstype == 'BIAS'] # .tolist()
        if len(bias_group) > 2:
            bias_confs = bias_group.groupby(['gain',
                                            'rdnoise',
                                            'radeg',
                                            'decdeg']).size().reset_index().rename(columns={0: 'count'})
            for i in bias_confs.index:
                bias_group = bias_group[((bias_group['gain'] == bias_confs.iloc[i]['gain']) &
                                        (bias_group['rdnoise'] == bias_confs.iloc[i]['rdnoise']) &
                                        (bias_group['radeg'] == bias_confs.iloc[i]['radeg']) &
                                        (bias_group['decdeg'] == bias_confs.iloc[i]['decdeg']))]
                self.data_container.add_bias(bias_group)
        else:
            log.error('Not enough bias images.')

        # flats separation
        flat_data = self.file_collection[self.file_collection.obstype == 'FLAT']
        # confs stands for configurations
        confs = flat_data.groupby(['object', 'filter']).size().reset_index().rename(columns={0: 'count'})
        for i in confs.index:
            flat_group = flat_data[((flat_data['object'] == confs.iloc[i]['object']) &
                                    (flat_data['filter'] == confs.iloc[i]['filter']))]
            self.data_container.add_day_flats(flat_group)
        # science data separation
        science_data = self.file_collection[self.file_collection.obstype == 'OBJECT']
        # confs stands for configurations
        confs = science_data.groupby(['object', 'filter']).size().reset_index().rename(columns={0: 'count'})
        for i in confs.index:
            science_group = science_data[((science_data['object'] == confs.iloc[i]['object']) &
                                          (science_data['filter'] == confs.iloc[i]['filter']))]
            self.data_container.add_data_group(science_group)

    # def get_twilight_time(self):
    #     """Get end/start time of evening/morning twilight
    #
    #     Notes:
    #         Taken from David Sanmartim's development
    #
    #     Returns:
    #         twilight_evening (str): Evening twilight time in the format 'YYYY-MM-DDTHH:MM:SS.SS'
    #         twilight_morning (str): Morning twilight time in the format 'YYYY-MM-DDTHH:MM:SS.SS'
    #         sun_set_time (str): Sun set time in the format 'YYYY-MM-DDTHH:MM:SS.SS'
    #         sun_rise_time (str): Sun rise time in the format 'YYYY-MM-DDTHH:MM:SS.SS'
    #
    #     """
    #     # observatory(str): Observatory name.
    #     observatory = 'SOAR Telescope'
    #     geodetic_location = ['-70d44m01.11s', '-30d14m16.41s', 2748]
    #     # longitude (str): Geographic longitude in string format
    #     longitude = geodetic_location[0]
    #     # latitude (str): Geographic latitude in string format.
    #     latitude = geodetic_location[1]
    #     # elevation (int): Geographic elevation in meters above sea level
    #     elevation = geodetic_location[2]
    #     # timezone (str): Time zone.
    #     timezone = 'UTC'
    #     # description(str): Observatory description
    #     description = 'Soar Telescope on Cerro Pachon, Chile'
    #
    #     soar_loc = EarthLocation.from_geodetic(longitude, latitude, elevation * u.m, ellipsoid='WGS84')
    #
    #     soar = Observer(name=observatory, location=soar_loc, timezone=timezone, description=description)
    #
    #     dateobs_list = self.file_collection['date-obs'].tolist()
    #     time_first_frame, time_last_frame = Time(min(dateobs_list)), Time(max(dateobs_list))
    #
    #     twilight_evening = soar.twilight_evening_astronomical(Time(time_first_frame), which='nearest').isot
    #     twilight_morning = soar.twilight_morning_astronomical(Time(time_last_frame), which='nearest').isot
    #     sun_set_time = soar.sun_set_time(Time(time_first_frame), which='nearest').isot
    #     sun_rise_time = soar.sun_rise_time(Time(time_last_frame), which='nearest').isot
    #     log.debug('Sun Set ' + sun_set_time)
    #     log.debug('Sun Rise ' + sun_rise_time)
    #     return twilight_evening, twilight_morning, sun_set_time, sun_rise_time

    def separate_day_night(self):
        """Separates day and night time data

        Notes:
            For day time separation it only considers the afternoon twilight
            since some observers will get data even past the morning twilight.
            All bias data goes into day time data.

        Returns:
            day_time_data (object):
            night_time_data (object):
        """
        # print(self.file_collection)
        date_obs_list = self.file_collection['date-obs'].tolist()
        afternoon_twilight, morning_twilight, sun_set, sun_rise = get_twilight_time(date_obs=date_obs_list)
        self.data_container.set_sun_times(sun_set, sun_rise)
        self.data_container.set_twilight_times(afternoon_twilight, morning_twilight)
        # print(afternoon_twilight, morning_twilight)
        day_time_data = self.file_collection[((self.file_collection['date-obs'] < afternoon_twilight)
                                              | (self.file_collection['date-obs'] > morning_twilight)
                                              | (self.file_collection['obstype'] == 'BIAS'))]
        night_time_data = self.file_collection[((self.file_collection['date-obs'] > afternoon_twilight)
                                                & (self.file_collection['date-obs'] < morning_twilight)
                                                & (self.file_collection['obstype'] != 'BIAS'))]
        # print(night_time_data)
        # print(day_time_data)
        return day_time_data, night_time_data


class Night(object):
    """This class is designed to be the organized data container. It doesn't
    store image data but list of pandas.DataFrame objects. Also it stores
    critical variables such as sunrise and sunset times.

    """

    def __init__(self, path, instrument, technique):
        """Initializes all the variables for the class

        Args:
            path (str): Full path to the directory where raw data is located
            instrument (str): 'Red' or 'Blue' stating whether the data was taken
            using the Red or Blue Goodman Camera.
            technique (str): 'Spectroscopy' or 'Imaging' stating what kind of
            data was taken.
        """

        self.full_path = path
        self.instrument = instrument
        self.technique = technique
        self.bias = None
        self.day_flats = None
        self.dome_flats = None
        self.sky_flats = None
        self.data_groups = None
        self.sun_set_time = None
        self.sun_rise_time = None
        self.evening_twilight = None
        self.morning_twilight = None

    def add_bias(self, bias_group):
        """Adds a bias group

        Args:
            bias_group (pandas.DataFrame): Contains a set of keyword values of grouped image metadata

        """

        if len(bias_group) < 2:
            if self.technique == 'Imaging':
                log.error('Imaging mode needs BIAS to work properly. Go find some.')
            else:
                log.warning('BIAS are needed for optimal results.')
        else:
            if self.bias is None:
                self.bias = [bias_group]
            else:
                self.bias.append(bias_group)

    def add_day_flats(self, day_flats):
        """"Adds a daytime flat group

        Args:
            day_flats (pandas.DataFrame): Contains a set of keyword values of grouped image metadata

        """

        if self.day_flats is None:
            self.day_flats = [day_flats]
        else:
            self.day_flats.append(day_flats)

    def add_data_group(self, data_group):
        """Adds a data group

        Args:
            data_group (pandas.DataFrame): Contains a set of keyword values of grouped image metadata

        """

        if self.data_groups is None:
            self.data_groups = [data_group]
        else:
            self.data_groups.append(data_group)

    def set_sun_times(self, sun_set, sun_rise):
        """Sets values for sunset and sunrise

        Args:
            sun_set (str): Sun set time in the format 'YYYY-MM-DDTHH:MM:SS.SS'
            sun_rise (str):Sun rise time in the format 'YYYY-MM-DDTHH:MM:SS.SS'

        """

        self.sun_set_time = sun_set
        self.sun_rise_time = sun_rise

    def set_twilight_times(self, evening, morning):
        """Sets values for evening and morning twilight

        Args:
            evening (str): Evening twilight time in the format 'YYYY-MM-DDTHH:MM:SS.SS'
            morning (str): Morning twilight time in the format 'YYYY-MM-DDTHH:MM:SS.SS'

        """

        self.evening_twilight = evening
        self.morning_twilight = morning





