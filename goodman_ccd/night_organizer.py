from __future__ import print_function
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

log = logging.getLogger('goodmanccd.nightorganizer')


class NightOrganizer(object):
    def __init__(self, args, night_dict):
        """

        Args:
            args (object):
            night_dict (dict):
        """

        self.args = args
        self.path = night_dict['full_path']
        self.instrument = night_dict['instrument']
        self.technique = night_dict['technique']
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
                         'filter2']
        self.file_collection = None
        self.all_datatypes = None
        self.data_container = Night(path=self.path, instrument=self.instrument, technique=self.technique)
        self.day_time_data = None
        self.night_time_data = None

    def __call__(self):
        """

        Returns:

        """

        ifc = ImageFileCollection(self.path, self.keywords)
        self.file_collection = ifc.summary.to_pandas()
        self.all_datatypes = self.file_collection.obstype.unique()
        if self.technique == 'Spectroscopy':
            self.day_time_data, self.night_time_data = self.separate_day_night()
            # print(len(self.day_time_data), len(self.night_time_data))
            if len(self.day_time_data) != 0:
                self.spectroscopy_day_time()
            else:
                log.warning('There is no day time data!')
                if not self.args.ignore_bias:
                    log.error('BIAS are needed for optimal results')
                    log.info('Check the argument --ignore-bias')
                    sys.exit('BIAS needed')
            if len(self.night_time_data) != 0:
                self.spectroscopy_night_time()
            else:
                log.warning('There is no night time data!')
                return False
        elif self.technique == 'Imaging':
            self.imaging_night()

        return self.data_container

    def spectroscopy_day_time(self):
        """
        Notes:
            cheking that bias and flats are more than 2 because otherwise will not be enough
        :return:
        """
        if len(self.day_time_data) != 0:
            # add bias
            bias_group = self.file_collection[self.file_collection.obstype == 'BIAS']
            if len(bias_group) > 2:
                self.data_container.add_bias(bias_group)
            else:
                log.error('Not enough bias images.')
            # add flats
            # TODO (simon): Add multivariable filtering for spectroscopy flats
            flat_data = self.file_collection[self.file_collection.obstype == 'FLAT']
            if len(flat_data) > 2:
                # configurations
                confs = flat_data.groupby(['grating',
                                           'filter2',
                                           'cam_targ',
                                           'grt_targ',
                                           'slit']).size().reset_index().rename(columns={0: 'count'})
                for i in confs.index:
                    # print(confs.iloc[i]['grating'])
                    flat_group = flat_data[((flat_data['grating'] == confs.iloc[i]['grating']) &
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
        """

        Returns:

        """

        # confs stands for configurations
        confs = self.night_time_data.groupby(['grating',
                                              'filter2',
                                              'cam_targ',
                                              'grt_targ',
                                              'slit',
                                              'obsra',
                                              'obsdec']).size().reset_index().rename(columns={0: 'count'})
        for i in confs.index:
            night_time_group = self.night_time_data[((self.night_time_data['grating'] == confs.iloc[i]['grating']) &
                                                     (self.night_time_data['filter2'] == confs.iloc[i]['filter2']) &
                                                     (self.night_time_data['cam_targ'] == confs.iloc[i]['cam_targ']) &
                                                     (self.night_time_data['grt_targ'] == confs.iloc[i]['grt_targ']) &
                                                     (self.night_time_data['slit'] == confs.iloc[i]['slit']) &
                                                     (self.night_time_data['obsra'] == confs.iloc[i]['obsra']) &
                                                     (self.night_time_data['obsdec'] == confs.iloc[i]['obsdec']))]
            self.data_container.add_data_group(night_time_group)

    def imaging_night(self):
        """

        Returns:

        """
        # TODO (simon): modify it to work with the day time data and nigh time data separation
        # bias data group
        afternoon_twilight, morning_twilight, sun_set, sun_rise = self.get_twilight_time()
        self.data_container.set_sun_times(sun_set, sun_rise)
        self.data_container.set_twilight_times(afternoon_twilight, morning_twilight)
        bias_group = self.file_collection[self.file_collection.obstype == 'BIAS'] # .tolist()
        self.data_container.add_bias(bias_group)

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

    def get_twilight_time(self):
        """Get end/start time of evening/morning twilight

        Notes:
            Taken from David Sanmartim's development

        Returns:
            twilight_evening (str): Evening twilight time in the format 'YYYY-MM-DDTHH:MM:SS.SS'
            twilight_morning (str): Morning twilight time in the format 'YYYY-MM-DDTHH:MM:SS.SS'

        """
        # observatory(str): Observatory name.
        observatory = 'SOAR Telescope'
        geodetic_location = ['-70d44m01.11s', '-30d14m16.41s', 2748]
        # longitude (str): Geographic longitude in string format
        longitude = geodetic_location[0]
        # latitude (str): Geographic latitude in string format.
        latitude = geodetic_location[1]
        # elevation (int): Geographic elevation in meters above sea level
        elevation = geodetic_location[2]
        # timezone (str): Time zone.
        timezone = 'UTC'
        # description(str): Observatory description
        description = 'Soar Telescope on Cerro Pachon, Chile'

        soar_loc = EarthLocation.from_geodetic(longitude, latitude, elevation * u.m, ellipsoid='WGS84')

        soar = Observer(name=observatory, location=soar_loc, timezone=timezone, description=description)

        dateobs_list = self.file_collection['date-obs'].tolist()
        time_first_frame, time_last_frame = Time(min(dateobs_list)), Time(max(dateobs_list))

        twilight_evening = soar.twilight_evening_astronomical(Time(time_first_frame), which='nearest').isot
        twilight_morning = soar.twilight_morning_astronomical(Time(time_last_frame), which='nearest').isot
        sun_set_time = soar.sun_set_time(Time(time_first_frame), which='nearest').isot
        sun_rise_time = soar.sun_rise_time(Time(time_last_frame), which='nearest').isot
        log.debug('Sun Set ' + sun_set_time)
        log.debug('Sun Rise ' + sun_rise_time)
        return twilight_evening, twilight_morning, sun_set_time, sun_rise_time

    def separate_day_night(self):
        """Separates day and night time data

        Notes:
            For day time separation it only considers the afternoon twilight since some observers will get data even
            past the morning twilight
            All bias data goes into day time data.

        Returns:
            day_time_data (object):
            night_time_data (object):
        """
        afternoon_twilight, morning_twilight, sun_set, sun_rise = self.get_twilight_time()
        self.data_container.set_sun_times(sun_set, sun_rise)
        self.data_container.set_twilight_times(afternoon_twilight, morning_twilight)
        # print(afternoon_twilight, morning_twilight)
        day_time_data = self.file_collection[((self.file_collection['date-obs'] < afternoon_twilight)
                                             | (self.file_collection['obstype'] == 'BIAS'))]
        night_time_data = self.file_collection[((self.file_collection['date-obs'] > afternoon_twilight)
                                                & (self.file_collection['date-obs'] < morning_twilight)
                                                & (self.file_collection['obstype'] != 'BIAS'))]
        # print(night_time_data)
        # print(daytime_data)
        return day_time_data, night_time_data

    @staticmethod
    def convert_time(in_time):
        """Converts time to seconds since epoch

        Args:
            in_time (str): time obtained from header's keyword DATE-OBS

        Returns:
            time in seconds since epoch

        """
        return time.mktime(time.strptime(in_time, "%Y-%m-%dT%H:%M:%S.%f"))


class Night(object):

    def __init__(self, path, instrument, technique):
        """

        Args:
            path:
            instrument:
            technique:
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
        self.afternoon_twilight = None
        self.morning_twilight = None

    def add_bias(self, bias_group):
        """

        Args:
            bias_group:

        Returns:

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
        """

        Args:
            day_flats:

        Returns:

        """

        if self.day_flats is None:
            self.day_flats = [day_flats]
        else:
            self.day_flats.append(day_flats)

    def add_data_group(self, data_group):
        """

        Args:
            data_group:

        Returns:

        """

        if self.data_groups is None:
            self.data_groups = [data_group]
        else:
            self.data_groups.append(data_group)

    def set_sun_times(self, sun_set, sun_rise):
        """

        Args:
            sun_set:
            sun_rise:

        Returns:

        """

        self.sun_set_time = sun_set
        self.sun_rise_time = sun_rise

    def set_twilight_times(self, afternoon, morning):
        """

        Args:
            afternoon:
            morning:

        Returns:

        """

        self.afternoon_twilight = afternoon
        self.morning_twilight = morning





