import argparse
import datetime
import glob
import logging

import os
import pymongo
import sys

from astropy.io import fits
from ccdproc import ImageFileCollection
from mongoengine import connect
from mongoengine import Document
from mongoengine import StringField, DateTimeField, FloatField, IntField
from mongoengine.errors import NotUniqueError

from ..core.core import NightDataContainer

from ..images.data_classifier import DataClassifier
from ..images.night_organizer import NightOrganizer
from ..images.image_processor import ImageProcessor

__version__ = __import__('goodman_pipeline').__version__


__arguments = ['--folder', os.getcwd()]


log = logging.getLogger(__name__)
log.setLevel(level=logging.DEBUG)


class FileRecord(Document):
    filename = StringField(required=True, unique=True, max_length=200)
    obstype = StringField(required=True)
    naxis = IntField(required=True)
    date = DateTimeField(required=True)
    slit = StringField(required=True)
    dateobs = DateTimeField(required=True)
    object = StringField(required=True)
    exptime = FloatField(required=True)
    obsra = StringField(required=True)
    obsdec = StringField(required=True)
    grating = StringField(required=True)
    cam_targ = FloatField(required=True)
    grt_targ = FloatField(required=True)
    filter = StringField(required=True)
    filter2 = StringField(required=True)
    gain = FloatField(required=True)
    rdnoise = FloatField(required=True)
    roi = StringField(required=True)
    wavmode = StringField(required=True)
    status = StringField(required=True, default='NEW')
    added = DateTimeField(default=datetime.datetime.now)


def get_args(arguments=None):
    parser = argparse.ArgumentParser(
        description='Reduce data live')

    parser.add_argument('--folder',
                        action='store',
                        dest='folder',
                        default=os.getcwd(),
                        help='Data location')

    parser.add_argument('--save-to',
                        action='store',
                        metavar='<save_to>',
                        type=str,
                        default=os.path.join(os.getcwd(), 'RED'),
                        help="Path to reduced data.")

    parser.add_argument('--reset',
                        action='store_true',
                        dest='reset',
                        help='Reset database')

    parser.add_argument('--ignore-bias',
                        action='store_true',
                        dest='ignore_bias',
                        help="Ignore bias correction")

    parser.add_argument('--ignore-flats',
                        action='store_true',
                        dest='ignore_flats',
                        help="Ignore flat field correction")

    parser.add_argument('--saturation',
                        action='store',
                        default=1.,
                        dest='saturation_threshold',
                        metavar='<value>',
                        help="Maximum percent of pixels above saturation "
                             "threshold. Default 1 percent.")

    parser.add_argument('--dbname',
                        action='store',
                        dest='dbname',
                        default='livereductiondb',
                        help='Database to store live reduction data')

    args = parser.parse_args(args=arguments)

    return args


class Reduce(object):

    def __init__(self, arguments=None):
        self.master = {
            'bias': None,
            'flat': None}
        self.data_containers_list = []
        self.args = get_args(arguments=arguments)

        self.dbclient = connect(self.args.dbname, host='localhost', port=27017)
        if self.args.reset:
            log.warning('Deleting {:d} database records'
                        ''.format(FileRecord.objects().count()))
            FileRecord.objects().delete()
            sys.exit(1)

        _file_list = []
        for _file_record in FileRecord.objects():
            _file_list.append(_file_record.filename)

        self.file_list = sorted(_file_list)
        self.data_classifier = DataClassifier()

        try:
            self.data_classifier(folder=self.args.folder)
            self.night_organizer = NightOrganizer(
                full_path=self.data_classifier.folder,
                instrument=self.data_classifier.instrument,
                technique=self.data_classifier.technique,
                ignore_bias=self.args.ignore_bias,
                ignore_flats=self.args.ignore_flats)
        except AttributeError:
            log.critical('Folder is empty')

        # self.db = self.dbclient[self.args.dbname]
        # self.files_in_path = self.db['files_in_path']

    def __update_args(self, new_args):
        pass

    def __insert_new_file_to_records(self, new_file_name, header):

        if FileRecord.objects(filename=new_file_name).count() == 0:
            if header['OBSTYPE'] == 'FOCUS':
                status = 'IGNORE'
            else:
                status = 'NEW'

            self.file_record = FileRecord(
                filename=new_file_name,
                obstype=header['OBSTYPE'],
                naxis=header['NAXIS'],
                date=header['DATE'],
                slit=header['SLIT'],
                dateobs=header['DATE-OBS'],
                object=header['OBJECT'],
                exptime=header['EXPTIME'],
                obsra=header['OBSRA'],
                obsdec=header['OBSDEC'],
                grating=header['GRATING'],
                cam_targ=header['CAM_TARG'],
                grt_targ=header['GRT_TARG'],
                filter=header['FILTER'],
                filter2=header['FILTER2'],
                gain=header['GAIN'],
                rdnoise=header['RDNOISE'],
                roi=header['ROI'],
                wavmode=header['WAVMODE'],
                status=status)

            self.file_record.save()

            log.debug('Saved new file {:s} with record id: {:s}'
                      ''.format(str(self.file_record.filename),
                                str(self.file_record.id)))
            return True
        else:
            log.debug('File {:s} already on database'.format(new_file_name))

        return False

    def __create_calibration_masters(self, calibration_type):
        if self.master[calibration_type] is None:
            if FileRecord.objects(type=calibration_type.upper()).count() > 3:
                log.info(
                    'Found {:d} {:s} files'.format(
                        FileRecord.objects(type=calibration_type.upper()).count(),
                        calibration_type))

                print(FileRecord.objects(type=calibration_type.upper()))

                process_data = ImageProcessor(
                    args=self.args,
                    data_container=self.data_containers_list)

                if calibration_type == 'bias':
                    for bias_group in self.data_containers_list.bias:
                        unique_obstype = bias_group.obstype.unique()
                        if len(unique_obstype) == 1 and \
                                unique_obstype[0] == 'BIAS' and \
                                not self.args.ignore_bias:
                            mb, mbn = process_data.create_master_bias(
                                bias_group=bias_group)

                            if self.master['bias'] is None:
                                self.master['bias'] = [mb]
                            else:
                                self.master['bias'].append(mb)

                if calibration_type == 'flat':
                    if self.master['bias'] is None:
                        log.error('Bias are required for making master {:s}'
                                  ''.format(calibration_type))

                print('Make {:s}'.format(calibration_type))
            else:
                log.warning('Not enough {:s} files'.format(calibration_type))
        else:
            log.debug('Master {:s} already exists'.format(calibration_type))

    def __call__(self, file_name=None):

        file_list = []

        if file_name is not None:
            file_list = [file_name]
        elif len(sorted(glob.glob(self.args.folder + '/*.fits'))) > 0:
            file_list = sorted(glob.glob(self.args.folder + '/*.fits'))

        new_files = [_file for _file in file_list if _file not in self.file_list]

        if new_files:
            for _file_name in new_files:
                log.debug(_file_name)
                header = fits.getheader(
                    os.path.join(self.args.folder, _file_name))

                self.__insert_new_file_to_records(new_file_name=_file_name,
                                                     header=header)
            print(FileRecord.objects(obstype='BIAS'))

            # if obstype in ['FLAT', 'OBJECT', 'COMP']:
            #     self.data_containers_list = self.night_organizer()[0]
            #     self.__create_calibration_masters(
            #         calibration_type='bias')
            #
            # if obstype in ['OBJECT', 'COMP']:
            #     self.__create_calibration_masters(
            #         calibration_type='flat')

        else:
            log.debug('Waiting for new files...')

        self.file_list = file_list








