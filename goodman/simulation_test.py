from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting
import logging as log
import argparse
import glob
import cPickle as pickle

from redspec import ScienceObject, MainApp
from process import Process, SciencePack



class TestSuite(object):

    def __init__(self):
        self.file_list = glob.glob('2016*fits')
        self.data = []
        self.headers = []
        self.args = self.get_args()

    def __call__(self, *args, **kwargs):
        for image in self.file_list:
            print(image)
            header = fits.getheader(image)
            # data = fits.getdata(image)
            name = header['OBJECT']
            obstime = header['DATE-OBS']
            ra = header['RA']
            dec = header['DEC']
            grating = header['GRATING']
            science_object = ScienceObject(name,
                                           image,
                                           obstime,
                                           ra,
                                           dec,
                                           grating)
            lamp_header = fits.getheader('fzh.0144_comp_400M2_GG455.fits')
            lamp_name = lamp_header['OBJECT']
            lamp_obstime = lamp_header['DATE-OBS']
            lamp_ra = lamp_header['RA']
            lamp_dec = lamp_header['DEC']
            lamp_grating = lamp_header['GRATING']
            science_object.add_lamp('fzh.0144_comp_400M2_GG455.fits',lamp_name, lamp_ra, lamp_dec)
            process = Process(science_object, self.args)
            extracted_data, rscience_object = process()
            # f = open('obj.save', 'rb')
            # loaded_obj = cPickle.load(f)
            # f.close()
            if isinstance(extracted_data, SciencePack):
                for data_index in range(len(extracted_data.data)):
                    plt.clf()
                    plt.title('Results %s' % name)
                    plt.plot(extracted_data.data[data_index], color='g', label='Extracted')
                    try:
                        print('Pickle File: ', image.replace('.fits', '_%s.pkl' % str(data_index + 1)))
                        f = open(image.replace('.fits', '_%s.pkl' % str(data_index + 1)), 'rb')
                        reference_spectrum = pickle.load(f)
                        f.close()
                        plt.plot(reference_spectrum, color='r', label='Reference')
                        diff = reference_spectrum - extracted_data.data[data_index]
                        plt.plot(diff, color='c', label='Difference (%.2f)' % np.mean(diff))
                    except:
                        print('Error')
                    plt.legend(loc='best')
                    plt.savefig('img/results-%s_%s.png' % (name, str(data_index + 1)), dpi=300)
                    # plt.show()
            else:
                print('No data extracted from this target %s' %name)




    @staticmethod
    def get_args():
        """Handles the argparse library and returns the arguments

        Returns:
            An object that contains all the variables parsed through the argument system
            The possible arguments to be returned are:

            -p or --data-path: has to be the source directory, where the data (images) is.
                    the location is **self.source**
                    default value is ./
            -d or --proc-path: is the destination where all new files/data will be placed.
                    the location is self.destiny
                    default value is ./
            -s or --search-pattern: the pattern that matches the reduced data that will be processed.
                    the location is self.pattern
                    default value is fc_
            -m or --obs-mode: is one of the predefined observing modes and the options are:
                    0: One or more lamps taken during the beginning or end of the night, i.e. single
                    calibration to all data in that night
                    1: One or more lamps around every science exposure.
                    2: An ASCII file will be read. This file contains a matching of science target files
                    to respective calibration lamp file that will be used for calibration.
                    3: No lamp used, use sky-lines instead.
                    the location is self.mode
                    default value is 0
            -l or --lamp-file: Name of the ASCII file that contains the relation between science files and lamp
                    files. An example is depicted below. Note that the lamps can be repeated.
                        #example of how the file should look
                        science_target_01.fits lamp_001.fits
                        science_target_02.fits lamp_001.fits
                        science_target_03.fits lamp_002.fits
                    the location is self.lamp_file
                    default value is lamps.txt

        Raises:
            In the case when -m or --obs-mode is set to 3 will requiere the name of file parsed with the -l or
            --lamp-file parameter an IOError is raised

        """
        leave = False
        parser = argparse.ArgumentParser('hola')

        parser.add_argument('-p', '--data-path',
                            action='store',
                            default='./',
                            type=str,
                            metavar='<Source Path>',
                            dest='source',
                            help='Path for location of raw data. Default <./>')

        parser.add_argument('-d', '--proc-path',
                            action='store',
                            default='./',
                            type=str,
                            metavar='<Destination Path>',
                            dest='destiny',
                            help='Path for destination of processed data. Default <./>')

        parser.add_argument('-s', '--search-pattern',
                            action='store',
                            default='fzh.',
                            type=str,
                            metavar='<Search Pattern>',
                            dest='pattern',
                            help="Pattern for matching the goodman's reduced data.")

        parser.add_argument('-m', '--obs-mode',
                            action='store',
                            default=0,
                            type=int,
                            metavar='<Observing Mode>',
                            dest='obsmode',
                            choices=[0, 1, 2, 3],
                            help='Defines the mode of matching lamps to science targets.')

        parser.add_argument('-r', '--reference-lamp',
                            action='store',
                            default='',
                            type=str,
                            metavar='<Lamp File>',
                            dest='lamp_all_night',
                            help="Name of reference lamp file for mode 0.\
                                 If not present, the first one in the list will be selected")

        parser.add_argument('-l', '--lamp-file',
                            action='store',
                            default='lamps.txt',
                            type=str,
                            metavar='<Lamp File>',
                            dest='lamp_file',
                            help="Name of an ASCII file describing which science target\
                                    uses which lamp. default <lamp.txt>")

        parser.add_argument('-t', '--telescope',
                            action='store_true',
                            default=False,
                            dest='telescope',
                            help="Enables the <Telescope> mode i.e. it run sequentially,\
                                    designed to use while observing at the telescope. Catches\
                                     new files arriving to the <source> folder. (NI!)")

        parser.add_argument('-i', '--interactive',
                            action='store_true',
                            default=False,
                            dest='interactive_ws',
                            help="Interactive wavelength solution. Disabled by default.")

        parser.add_argument('-o', '--output-prefix',
                            action='store',
                            default='g',
                            dest='output_prefix',
                            help="Prefix to add to calibrated spectrum.")

        parser.add_argument('-R', '--reference-files',
                            action='store',
                            default='refdata/',
                            dest='reference_dir',
                            help="Reference files location")

        parser.add_argument('--plots-enabled',
                            action='store_true',
                            default=False,
                            dest='plots_enabled',
                            help="Show plots of intermediate steps. For debugging only.")

        args = parser.parse_args()

        return args

if __name__ == '__main__':
    tester = TestSuite()
    tester()