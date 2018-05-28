#!/usr/bin/env python2
# -*- coding: utf8 -*-
"""Pipeline for GOODMAN spectra Extraction.

This program finds reduced images, i.e. trimmed, bias subtracted, flat fielded,
etc. that match the <pattern> in the source folder, then classify them in two
groups: Science or Lamps. For science images, finds the spectrum or spectra and
traces it doing some fit.
Simon Torres 2016-06-28

"""
# TODO (simon): Change all astropy.io.fits to astropy.CCDData.read
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .wavelength import WavelengthCalibration
from ..core import (classify_spectroscopic_data,
                    search_comp_group,
                    add_wcs_keys,
                    identify_targets,
                    trace_targets,
                    extraction,
                    save_extracted)

from ..core import (NoMatchFound,
                    NoTargetException,
                    ReferenceData)

import sys
import os
import textwrap
import argparse
import astropy.units as u
import glob
import logging
from ccdproc import CCDData

# import matplotlib
# matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import warnings

SHOW_PLOTS = False

warnings.filterwarnings('ignore')

__version__ = __import__('pipeline').__version__

def get_args(arguments=None):
    """Handles the argparse library and returns the arguments

    The list of arguments can be found with running `redspec -h`.

    Notes:
        The full list of arguments are not listed here as the may evolve in
        which case is impossible to keep this up to date.


    Returns:
        An object that contains all the variables parsed through the argument
        system

    """

    parser = argparse.ArgumentParser(
        description="Extracts goodman spectra and does automatic wavelength "
                    "calibration.\nPipeline Version: {:s}".format(__version__))

    parser.add_argument('--data-path',
                        action='store',
                        default=os.getcwd(),
                        type=str,
                        metavar='<Source Path>',
                        dest='source',
                        help='Path for location of raw data. Default <./>')

    parser.add_argument('--proc-path',
                        action='store',
                        default=os.getcwd(),
                        type=str,
                        metavar='<Destination Path>',
                        dest='destination',
                        help='Path for destination of processed data. Default '
                             '<./>')

    parser.add_argument('--search-pattern',
                        action='store',
                        default='cfzsto',
                        type=str,
                        metavar='<Search Pattern>',
                        dest='pattern',
                        help="Pattern for matching the goodman's reduced data.")

    parser.add_argument('--output-prefix',
                        action='store',
                        default='w',
                        metavar='<Out Prefix>',
                        dest='output_prefix',
                        help="Prefix to add to calibrated spectrum.")

    parser.add_argument('--extraction',
                        action='store',
                        default='fractional',
                        type=str,
                        metavar='<Extraction Type>',
                        dest='extraction_type',
                        choices=['fractional', 'optimal'],
                        help='Only fractional pixel extraction is implemented.')

    parser.add_argument('--reference-files',
                        action='store',
                        default='data/ref_comp/',
                        metavar='<Reference Dir>',
                        dest='reference_dir',
                        help="Directory of Reference files location")

    parser.add_argument('--debug',
                        action='store_true',
                        dest='debug_mode',
                        help="Debugging Mode")

    parser.add_argument('--max-targets',
                        action='store',
                        dest='max_n_targets',
                        metavar='<max targets>',
                        type=int,
                        default=3,
                        help="Maximum number of targets to be found in a "
                             "single image. Default 3")

    parser.add_argument('--save-plots',
                        action='store_true',
                        dest='save_plots',
                        help="Save all plots in a directory")

    parser.add_argument('--plot-results',
                        action='store_true',
                        dest='plot_results',
                        help="Show wavelength calibrated spectrum at the end.")

    args = parser.parse_args(args=arguments)

    return args


class MainApp(object):
    """Defines and initialize all important variables for processing the data

    The MainApp class controls the way the night is organized for further
    processing. It also sets the appropriate parameters that will allow for a
    smooth working in all the other modules.

    """
    def __init__(self):
        """Init method for MainApp class

        This method initializes the arguments for the class, if they are not
        provided it will get them.


        """
        self.log = logging.getLogger(__name__)
        self.args = None
        self.wavelength_solution_obj = None
        self.wavelength_calibration = None
        self.reference = None
        self._pipeline_version = __version__

    def __call__(self, args=None):
        """Call method for the MainApp class

        This method call the higher level functions in order to do the
        spectroscopic data reduction.

        Args:
            args (object): argparse.Namespace instance that contains all the
                arguments.

        """

        if args is None:
            self.args = get_args()
        else:
            self.args = args

        if not self._check_args():
            sys.exit()
      
        self.log.debug("Initializing reference data locator.")
        self.reference = ReferenceData(reference_dir=self.args.reference_dir)
        self.log.debug("Calling data classification procedure.")
        data_container = classify_spectroscopic_data(
            path=self.args.source,
            search_pattern=self.args.pattern)

        if data_container.is_empty:
            self.log.debug("Received empty data container.")
            sys.exit("Unable to find or classify data.")
        else:
            self.log.debug("Received non-empty data container.")

        self.log.debug("Calling _run method for MainApp")
        self._run(data_container=data_container,
                  extraction_type=self.args.extraction_type)

        self.log.info("END")

    def _run(self, data_container, extraction_type):
        assert data_container.is_empty is False
        assert any(extraction_type == option for option in ['fractional',
                                                            'optimal'])
        full_path = data_container.full_path

        for sub_container in [groups for groups in [
            data_container.spec_groups,
            data_container.object_groups]
                              if groups is not None]:
            for group in sub_container:
                self.wavelength_calibration = WavelengthCalibration(
                    args=self.args)

                object_group = group[group.obstype == 'OBJECT']
                obj_groupby = object_group.groupby(['object']).size(

                ).reset_index().rename(columns={0: 'count'})

                self.log.info("Processing Science Target: "
                              "{:s} with {:d} files."
                              "".format(obj_groupby.iloc[0]['object'],
                                        obj_groupby.iloc[0]['count']))
                # this has to be initialized here
                comp_group = None
                comp_ccd_list = []
                if 'COMP' in group.obstype.unique():
                    self.log.debug('Group has comparison lamps')
                    comp_group = group[group.obstype == 'COMP']
                    comp_group = self.reference.check_comp_group(comp_group)

                if comp_group is None:
                    self.log.debug('Group does not have comparison lamps')
                    if data_container.comp_groups is not None:
                        self.log.debug('There are comparison lamp group '
                                       'candidates')

                        try:

                            comp_group = search_comp_group(
                                object_group=object_group,
                                comp_groups=data_container.comp_groups,
                                reference_data=self.reference)

                            self.log.warning(
                                'This comparison lamp might not be optimal '
                                'if you are doing radial velocity studies')

                        except NoMatchFound:

                            self.log.error(
                                'It was not possible to find a comparison '
                                'group')
                    else:
                        self.log.warning('Data will be extracted but not '
                                         'calibrated')

                _combine = True

                if len(object_group.file.tolist()) > 1 and _combine:
                    self.log.debug("This can be combined")

                for spec_file in object_group.file.tolist():

                    self.log.info('Processing Science File: {:s}'.format(
                        spec_file))

                    file_path = os.path.join(full_path, spec_file)
                    ccd = CCDData.read(file_path, unit=u.adu)
                    ccd.header.set('GSP_PNAM', value=spec_file)
                    ccd = add_wcs_keys(ccd=ccd)

                    if comp_group is not None and comp_ccd_list == []:

                        for comp_file in comp_group.file.tolist():
                            comp_path = os.path.join(full_path, comp_file)
                            comp_ccd = CCDData.read(comp_path, unit=u.adu)
                            comp_ccd = add_wcs_keys(ccd=comp_ccd)
                            comp_ccd.header.set('GSP_PNAM', value=comp_file)
                            comp_ccd_list.append(comp_ccd)

                    else:

                        self.log.debug(
                            'Comp Group is None or comp list already exist')

                    # identify
                    self.log.debug("Calling procedure for target "
                                   "identification.")
                    target_list = identify_targets(ccd=ccd,
                                                   nfind=3,
                                                   plots=SHOW_PLOTS)

                    # trace
                    if len(target_list) > 0:
                        self.log.debug("Calling procedure for tracing target.")
                        trace_list = trace_targets(ccd=ccd,
                                                   target_list=target_list,
                                                   sampling_step=5,
                                                   pol_deg=2)
                    else:
                        self.log.error("The list of identified targets is "
                                       "empty.")
                        continue

                    # if len(trace_list) > 0:
                    extracted_target_and_lamps = []
                    for single_trace, single_profile in trace_list:
                        try:
                            # target extraction
                            extracted = extraction(
                                ccd=ccd,
                                target_trace=single_trace,
                                spatial_profile=single_profile,
                                extraction_name=extraction_type)
                            save_extracted(ccd=extracted,
                                           destination=self.args.destination)

                            # lamp extraction
                            all_lamps = []
                            if comp_ccd_list:
                                for comp_lamp in comp_ccd_list:
                                    extracted_lamp = extraction(
                                        ccd=comp_lamp,
                                        target_trace=single_trace,
                                        spatial_profile=single_profile,
                                        extraction_name=extraction_type)
                                    save_extracted(
                                        ccd=extracted_lamp,
                                        destination=self.args.destination)
                                    all_lamps.append(extracted_lamp)
                            extracted_target_and_lamps.append([extracted,
                                                               all_lamps])

                            if self.args.debug_mode:
                                plt.close('all')
                                fig, ax = plt.subplots(1, 1)

                                fig.canvas.set_window_title('Extracted Data')

                                manager = plt.get_current_fig_manager()

                                if plt.get_backend() == u'GTK3Agg':
                                    manager.window.maximize()
                                elif plt.get_backend() == u'Qt5Agg':
                                    manager.window.showMaximized()

                                for edata, comps in extracted_target_and_lamps:
                                    ax.plot(edata.data,
                                            label=edata.header['OBJECT'])
                                    if comps:
                                        for comp in comps:
                                            ax.plot(comp.data,
                                                    label=comp.header[
                                                         'OBJECT'])
                                ax.legend(loc='best')
                                if plt.isinteractive():
                                    plt.draw()
                                    plt.pause(1)
                                else:
                                    plt.show()

                        except NoTargetException:
                            self.log.error('No target was identified in file'
                                           ' {:s}'.format(spec_file))
                            continue
                    object_number = None
                    for sci_target, comp_list in extracted_target_and_lamps:
                        self.wavelength_solution_obj = \
                            self.wavelength_calibration(
                                ccd=sci_target,
                                comp_list=comp_list,
                                object_number=object_number)

        return True

    def _check_args(self):
        """Check arguments for redspec

        Returns:
            True if all the test succeed or False in case of any of them fails.

        """

        self.log.debug("Starting arguments check")

        try:
            ref_full_path = os.path.join(
                os.path.dirname(sys.modules['goodman.pipeline'].__file__),
                self.args.reference_dir)
        except KeyError as error:
            self.log.debug("KeyError {:s}".format(str(error)))
            ref_full_path = os.path.join(
                os.path.dirname(sys.modules['pipeline'].__file__),
                self.args.reference_dir)
        if not os.path.isdir(ref_full_path):
            self.log.critical("Reference files directory \"{:s}\" "
                              "doesn't exist.".format(ref_full_path))
            return False
        else:
            self.args.reference_dir = ref_full_path

        if not os.path.isabs(self.args.source):
            self.args.source = os.path.join(os.getcwd(), self.args.source)

        if os.path.isdir(self.args.source):
            if len(glob.glob(os.path.join(self.args.source,
                                          self.args.pattern + '*'))) == 0:
                self.log.critical("Source directory is empty or does not "
                                  "contain any files "
                                  "matching \"{:s}\" "
                                  "".format(self.args.pattern + "*"))
                return False
        else:
            self.log.critical(
                "Source Directory \"{:s}\" doesn't exist."
                "".format(self.args.source))
            return False

        if not os.path.isabs(self.args.destination):
            self.args.destination = os.path.join(os.getcwd(),
                                                 self.args.destination)

        if not os.path.isdir(self.args.destination):
            self.log.warning("Destination folder \"{:s}\" doesn't exist."
                             "".format(self.args.destination))
            try:
                os.path.os.makedirs(self.args.destination)
                self.log.info('Destination folder \"{:s}\" created.'
                              ''.format(self.args.destination))
            except OSError as err:
                self.log.critical('OSError: {:s}'.format(repr(err)))
                return False
        elif not os.access(self.args.destination, os.W_OK):
            self.log.critical('Destination folder \"{:s}\" does not have '
                              'writing permission.'
                              ''.format(self.args.destination))
            return False
        return True


if __name__ == '__main__':
    MAIN_APP = MainApp()
    try:
        MAIN_APP()
    except KeyboardInterrupt:
        sys.exit(0)
