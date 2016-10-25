"""Test for Process

Returns:
    True value

"""
from astropy.io import fits
from redspec import ScienceObject
from process import Process
import logging as log
import argparse
import textwrap
from sphinx import cmdline

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
                        1: One or more lamps right before OR right after every science exposure.
                        2: One or more lamps right before AND right after every science exposure.
                        3: An ASCII file will be read. This file contains a matching of sience target files
                        to respective calibration lamp file that will be used for calibration.
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
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                         description=textwrap.dedent(
                                             '''Extracts goodman spectra and does wavelength calibration.\n\n\
Supported Observing modes are:
    <0>: (Default) reads lamps taken at the begining or end of the night.\n\
    <1>: one or more lamps around science exposure.\n\
    <2>: ASCII file describing which science target uses which lamp.\n\
    <3>: No lamps. Uses the sky lines
    '''))

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
                        default='fc_',
                        type=str,
                        metavar='<Search Pattern>',
                        dest='pattern',
                        help="Pattern for matching the goodman's reduced data.")

    parser.add_argument('-m', '--obs-mode',
                        action='store',
                        default=0,
                        type=int,
                        metavar='<Observing Mode>',
                        dest='mode',
                        choices=[0, 1, 2, 3],
                        help='Defines the mode of matching lamps to science targets.')

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
                                         new files arriving to the <source> folder.")

    parser.add_argument('-i', '--interactive',
                        action='store_true',
                        default=False,
                        dest='interactive_ws',
                        help="Interactive wavelength solution")

    args = parser.parse_args()

    return args


single = True

if single:
    sourcepath = '/data/simon/data/soar/work/goodman/test/'
    obj = ScienceObject('R339704',
                        'fc_0045.SO2016A-019_0320.fits',
                        '2016-03-20T23:54:15.96',
                        65.69922083333333,
                        -76.57505888888889)
    # obj.add_lamp('fc_0046.SO2016A-019_0320.fits', 'HgAr', 65.69922083333333, -76.57505888888889)
    obj.add_lamp('fc_0047.SO2016A-019_0320.fits', 'CuHeAr', 65.69922083333333, -76.57505888888889)
else:
    sourcepath = '/user/simon/data/SOAR/work/multi-objects/'
    obj = ScienceObject('CVSO-34',
                        'bc_0099.cvso34_400M2_GG455.fits',
                        '2015-01-04T05:33:57.22',
                        81.41224583333333,
                        1.7139411111111111)
    obj.add_lamp('bc_0110.near_400M2_GG455.fits', 'NeAr', 81.45332916666668, 1.8253911111111112)
    obj.add_lamp('bc_0111.near_400M2_GG455.fits', 'NeAr', 81.45332916666668, 1.8253911111111112)

obj.print_all()

args = get_args()

p = Process(sourcepath,obj, args)
#targets = p.identify_spectra()
#for target in targets:
#    print(target)