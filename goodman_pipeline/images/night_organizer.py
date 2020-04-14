from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import sys
import pandas
import logging
from ccdproc import ImageFileCollection
from ..core import get_twilight_time, ra_dec_to_deg
from ..core import NightDataContainer

# log = logging.getLogger(__name__)


class NightOrganizer(object):

    def __init__(self, full_path, instrument, technique, ignore_bias=False,
                 ignore_flats=False):
        """Initializes the NightOrganizer class

        This class contains methods to organize the data for processing. It will
        identify groups of OBJECTS, FLATS or COMPS (comparison lamps) whenever
        they exist. The product will be an object that will act as a data
        container.

        Args:
            full_path (str): Full path to data location
            instrument (str): Instrument name, it refers to the camera used,
              for pipeline it can be Blue or Red
            technique (str): Technique used to obtain the data, Spectroscopy or
              Imaging.
            ignore_bias (bool): Flag that allows to bypass the processing of
              Bias.
            ignore_flats (bool): Flag that allows to bypass the processing of
              Flats.

        """
        self.log = logging.getLogger(__name__)
        self.path = full_path
        self.instrument = instrument
        self.technique = technique
        self.ignore_bias = ignore_bias
        self.ignore_flats = ignore_flats
        self.file_collection = None
        self.all_datatypes = None
        self.day_time_data = None
        self.night_time_data = None
        self.data_container = NightDataContainer(path=self.path,
                                                 instrument=self.instrument,
                                                 technique=self.technique)
        self.keywords = ['naxis',
                         'date',
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
                         'rdnoise',
                         'roi',
                         'wavmode']

    def __call__(self):
        """Call method

        Creates a table with selected keywords that will allow to group the data
        in order to be classified according to the observational technique used,
        imaging or spectroscopy.

        Returns:
            data_container_list (list): List of instances of NightDataContainer
            which is the class used as storage unit for classified data. The
            number of elements depends on the number of combinations of `GAIN`,
            `RDNOISE` and `ROI`.

        """

        ifc = ImageFileCollection(self.path, self.keywords)
        self.file_collection = ifc.summary.to_pandas()

        self.check_header_cards()
        if 3 in self.file_collection['naxis'].unique():
            raise IOError('One of the files of the night has a shape of a '
                          'three dimensional array. When images should have '
                          'two.')

        # if technique is Spectroscopy, ignore all data that has
        # WAVMODE = Imaging because assumes they are acquisition images
        if self.technique == 'Spectroscopy':
            _imaging_file = self.file_collection[
                self.file_collection.wavmode == 'IMAGING']
            if not _imaging_file.empty:
                self.log.warning("Ignoring all Imaging data. Assuming they are "
                                 "not science exposures.")
                self.log.info("If any of the files listed below is a science "
                              "file, please process it/them in a separate "
                              "folder.")

                for _file in _imaging_file['file']:
                    self.log.info("Discarding image: {:s}".format(_file))

            self.file_collection = self.file_collection[
                self.file_collection.wavmode != 'IMAGING'].reset_index(drop=True)

        elif self.technique == 'Imaging':
            self.log.warning("Ignoring all files where `wavmode` is not "
                             "Imaging.")
            self.file_collection = self.file_collection[
                self.file_collection.wavmode == 'IMAGING'].reset_index(drop=True)

        # add two columns that will contain the ra and dec in degrees

        self.file_collection['radeg'] = ''
        self.file_collection['decdeg'] = ''
        for i in self.file_collection.index.tolist():
            radeg, decdeg = ra_dec_to_deg(self.file_collection.obsra.iloc[i],
                                          self.file_collection.obsdec.iloc[i])

            self.file_collection.iloc[
                i, self.file_collection.columns.get_loc('radeg')] = \
                '{:.2f}'.format(radeg)

            self.file_collection.iloc[
                i, self.file_collection.columns.get_loc('decdeg')] = \
                '{:.2f}'.format(decdeg)
            # now we can compare using degrees

        readout_configurations = self.file_collection.groupby(
            ['gain',
             'rdnoise',
             'roi']).size().reset_index().rename(columns={0: 'count'})

        data_container_list = []
        for i in readout_configurations.index:
            self.log.info("Organizing data for this configuration: "
                          "Gain: {:.2f}, Noise: {:.2f}, ROI: {:s}"
                          "".format(readout_configurations.iloc[i]['gain'],
                                    readout_configurations.iloc[i]['rdnoise'],
                                    readout_configurations.iloc[i]['roi']))
            if not self.data_container.is_empty:
                self.log.debug("Reset data container")
                self.data_container = NightDataContainer(
                    path=self.path,
                    instrument=self.instrument,
                    technique=self.technique)

            self.data_container.set_readout(
                gain=readout_configurations.iloc[i]['gain'],
                rdnoise=readout_configurations.iloc[i]['rdnoise'],
                roi=readout_configurations.iloc[i]['roi'])

            sub_collection = self.file_collection[
                ((self.file_collection['gain'] ==
                  readout_configurations.iloc[i]['gain']) &
                 (self.file_collection['rdnoise'] ==
                  readout_configurations.iloc[i]['rdnoise']) &
                 (self.file_collection['roi'] ==
                  readout_configurations.iloc[i]['roi']))]

            self.all_datatypes = sub_collection.obstype.unique()
            if self.technique == 'Spectroscopy':
                self.spectroscopy_night(file_collection=sub_collection,
                                        data_container=self.data_container)
            elif self.technique == 'Imaging':
                self.imaging_night()

            if self.data_container.is_empty:
                self.log.warning("The following files will not be processed:")
                for _file in sub_collection['file'].tolist():
                    self.log.warning("{:s}".format(_file))
                self.log.debug('data_container is empty')
            else:
                self.log.info('Found valid data, appending to data container '
                              'list')
                data_container_list.append(self.data_container)

        # Warn the user in case the list of data_container element is empty or
        # all the elements are None
        # print(data_container_list)
        if len(data_container_list) == 0:
            return [None]
        elif not all(data_container_list):
            self.log.warning("It is possible that there is no valid data.")
            return [None]
        else:
            return data_container_list

    def check_header_cards(self):
        """Check if the header contains all the keywords (cards) expected.

        This is critical for old goodman data.

        Raises:
            ValueError: If any of the cards does not exists in the image's
              header
        """

        missing_cards = []
        for card in self.keywords:

            if self.file_collection[card].isnull().values.any():
                missing_cards.append(card.upper())

        if len(missing_cards) != 0:
            missing_cards = ', '.join(missing_cards)
            raise ValueError(
                "{:} ".format(missing_cards) +
                "card(s) is(are) not found in one of the headers. Please, "
                "check your data. A script is being developed to correct this "
                "automatically but, for now, you will have to add this "
                "keyword manually."
                )

    def spectroscopy_night(self, file_collection, data_container):
        """Organizes data for spectroscopy

        This method identifies all combinations of nine **key** keywords that
        can set apart different objects with their respective calibration data
        or not. The keywords used are:
          - ``GAIN``
          - ``RDNOISE``
          - ``GRATING``
          - ``FILTER2``
          - ``CAM_TARG``
          - ``GRT_TARG``
          - ``SLIT``
          - ``OBSRA``
          - ``OBSDEC``

        This method populates the `data_container` class attribute which is an
        instance of the :class:`goodman_pipeline.core.core.NightDataContainer`.
        A data group is an instance of a :class:`pandas.DataFrame`.

        """

        assert isinstance(file_collection, pandas.DataFrame)
        assert isinstance(data_container, NightDataContainer)

        # obtain a list of timestamps of observing time
        # this will only be used for naming flats
        dateobs_list = file_collection['date-obs'].tolist()

        # get times for twilights, sunset an sunrise
        afternoon_twilight, morning_twilight, sun_set, sun_rise = \
            get_twilight_time(date_obs=dateobs_list)

        # set times in data container
        data_container.set_sun_times(sun_set,
                                     sun_rise)
        data_container.set_twilight_times(afternoon_twilight,
                                          morning_twilight)

        # process bias
        bias_collection = file_collection[file_collection.obstype == 'BIAS']

        if not self.ignore_bias:
            if len(bias_collection) == 0:
                self.log.critical('There is no BIAS images for this '
                                  'configuration. Use --ignore-bias '
                                  'to proceed without BIAS.')
                # sys.exit('CRITICAL ERROR: BIAS not Found.')
                return False
            else:
                bias_conf = bias_collection.groupby(
                    ['gain',
                     'rdnoise',
                     'radeg',
                     'decdeg']).size().reset_index().rename(
                    columns={0: 'count'})

                # bias_conf
                for i in bias_conf.index:
                    bias_group = bias_collection[(
                        (bias_collection['gain'] == bias_conf.iloc[i]['gain']) &
                        (bias_collection['rdnoise'] == bias_conf.iloc[i][
                            'rdnoise']) &
                        (bias_collection['radeg'] == bias_conf.iloc[i][
                            'radeg']) &
                        (bias_collection['decdeg'] == bias_conf.iloc[i][
                            'decdeg']))]

                    data_container.add_bias(bias_group=bias_group)
        else:
            self.log.warning('Ignoring BIAS by request.')

        if not any([value in ['FLAT', 'LAMPFLAT'] for value in file_collection.obstype.unique()]) and \
                not self.ignore_flats:
            self.log.critical('There is no FLAT images. Use --ignore-flats to '
                              'continue without FLATs.')
            sys.exit('CRITICAL ERROR: FLAT not Found.')
        elif self.ignore_flats:
            self.log.warning('Ignoring FLAT images on request.')
            data_collection = file_collection[
                ((file_collection.obstype != 'BIAS') &
                 (file_collection.obstype != 'FLAT'))]
        else:
            # process non-bias i.e. flats and object ... and comp
            data_collection = file_collection[file_collection.obstype != 'BIAS']

        confs = data_collection.groupby(
            ['gain',
             'rdnoise',
             'grating',
             'filter2',
             'cam_targ',
             'grt_targ',
             'slit',
             'radeg',
             'decdeg']).size().reset_index().rename(columns={0: 'count'})

        for i in confs.index:

            data_group = data_collection[(
                (data_collection['gain'] == confs.iloc[i]['gain']) &
                (data_collection['rdnoise'] == confs.iloc[i]['rdnoise']) &
                (data_collection['grating'] == confs.iloc[i]['grating']) &
                (data_collection['filter2'] == confs.iloc[i]['filter2']) &
                (data_collection['cam_targ'] == confs.iloc[i]['cam_targ']) &
                (data_collection['grt_targ'] == confs.iloc[i]['grt_targ']) &
                (data_collection['slit'] == confs.iloc[i]['slit']) &
                (data_collection['radeg'] == confs.iloc[i]['radeg']) &
                (data_collection['decdeg'] == confs.iloc[i]['decdeg']))]

            group_obstype = data_group.obstype.unique()

            if any([value in ['FLAT', 'LAMPFLAT'] for value in group_obstype]) and len(group_obstype) == 1:

                data_container.add_day_flats(data_group)

            else:
                # Comparison lamps are processed as science data.
                data_container.add_spec_group(data_group)

                if any([value in ['FLAT', 'LAMPFLAT'] for value in group_obstype]):
                    # grab flats and put them in the flats group as well
                    object_flat_group = data_group[(
                        (data_group['obstype'] == 'FLAT') |
                        (data_group['obstype'] == 'LAMPFLAT'))]
                    data_container.add_day_flats(object_flat_group)

        return data_container

    def imaging_night(self):
        """Organizes data for imaging

        For imaging there is no discrimination regarding night data since the
        process is simpler. It is a three stage process classifying ``BIAS``,
        ``FLAT``  and ``OBJECT`` data type. The data is packed in groups that
        are :class:`pandas.DataFrame` objects.

        """

        # bias data group
        date_obs_list = self.file_collection['date-obs'].tolist()

        afternoon_twilight, morning_twilight, sun_set, sun_rise = \
            get_twilight_time(date_obs=date_obs_list)

        self.data_container.set_sun_times(sun_set=sun_set,
                                          sun_rise=sun_rise)

        self.data_container.set_twilight_times(evening=afternoon_twilight,
                                               morning=morning_twilight)

        bias_group = self.file_collection[
            self.file_collection.obstype == 'BIAS']  # .tolist()

        if len(bias_group) > 2:

            bias_confs = bias_group.groupby(
                ['gain',
                 'rdnoise',
                 'radeg',
                 'decdeg']).size().reset_index().rename(columns={0: 'count'})

            for i in bias_confs.index:

                bias_group = bias_group[(
                    (bias_group['gain'] == bias_confs.iloc[i]['gain']) &
                    (bias_group['rdnoise'] == bias_confs.iloc[i]['rdnoise']) &
                    (bias_group['radeg'] == bias_confs.iloc[i]['radeg']) &
                    (bias_group['decdeg'] == bias_confs.iloc[i]['decdeg']))]

                self.data_container.add_bias(bias_group)
        else:
            self.log.error('Not enough bias images.')
            sys.exit("Bias required for imaging.")

        # flats separation
        flat_data = self.file_collection[(
                (self.file_collection.obstype == 'FLAT') |
                (self.file_collection.obstype == 'LAMPFLAT'))]

        # confs stands for configurations
        confs = flat_data.groupby(
            ['object',
             'filter']).size().reset_index().rename(columns={0: 'count'})

        for i in confs.index:

            flat_group = flat_data[
                ((flat_data['object'] == confs.iloc[i]['object']) &
                 (flat_data['filter'] == confs.iloc[i]['filter']))]

            self.data_container.add_day_flats(flat_group)

        # science data separation
        science_data = self.file_collection[(
            (self.file_collection.obstype == 'OBJECT') |
            (self.file_collection.obstype == 'EXPOSE'))]

        # confs stands for configurations
        confs = science_data.groupby(
            ['object',
             'filter']).size().reset_index().rename(columns={0: 'count'})

        for i in confs.index:

            science_group = science_data[
                ((science_data['object'] == confs.iloc[i]['object']) &
                 (science_data['filter'] == confs.iloc[i]['filter']))]

            self.data_container.add_data_group(science_group)
