# -*- coding: utf8 -*-
""" Line list per elements

This contains a dictionary of elements used for comparison lamps with their
main emission lines

The current elements present are:

    Hg: Mercury
    Ar: Argon
    Cu: Copper
    Ne: Neon
    He: Helium
    CuAr: Copper Argon

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import logging
import pandas
import os
import ccdproc
from ccdproc import CCDData
import astropy.units as u
import sys
from glob import glob
import re


# FORMAT = '%(levelname)s:%(filename)s:%(module)s: 	%(message)s'
# self.log.basicConfig(level=self.log.DEBUG, format=FORMAT)


class ReferenceData(object):
    """Contains spectroscopic reference lines values and filename to templates.

    This class stores:
        - file names for reference fits spectrum
        - file names for CSV tables with reference lines and relative
          intensities
        - line positions only for the elements used in SOAR comparison lamps
    """
    def __init__(self, reference_dir):
        """Init method for the ReferenceData class

        This methods uses ccdproc.ImageFileCollection on the reference_dir to
        capture all possible reference lamps. Also defines dictionaries
        containg line lists for several elements used in lamps.

        Args:
            reference_dir (str): full path to the reference data directory
        """
        self.log = logging.getLogger(__name__)
        self.reference_dir = reference_dir
        reference_collection = ccdproc.ImageFileCollection(self.reference_dir)
        self.ref_lamp_collection = reference_collection.summary.to_pandas()
        # print(self.ref_lamp_collection)
        self.lines_pixel = None
        self.lines_angstrom = None
        self._ccd = None

    def get_reference_lamp(self, header):
        """Finds a suitable template lamp from the catalog

        Args:
            header (object): FITS header of image we are looking a a reference
                lamp.

        Returns:
            full path to best matching reference lamp.

        """
        # print(header['CAM_TARG'], self.ref_lamp_collection['cam_targ'].to_string(index=False))
        # print(header['GRT_TARG'],
        #       self.ref_lamp_collection['grt_targ'].to_string(index=False))
        # print(header['wavmode'],
        #       self.ref_lamp_collection['wavmode'].to_string(index=False))
        # print(header['object'],
        #       self.ref_lamp_collection['object'].to_string(index=False))
        # print(header['GSP_FNAM'])

        filtered_collection = self.ref_lamp_collection[
            (self.ref_lamp_collection['object'] == header['object']) &
            (self.ref_lamp_collection['wavmode'] == header['wavmode'])]
        # print(filtered_collection)
        if filtered_collection.empty:
            raise NotImplementedError("It was not possible to find any lamps "
                                      "that match")
        elif len(filtered_collection) == 1:
            self.log.info("Lamp: {:s}".format(filtered_collection.file.to_string(
                                         index=False)))
            full_path = os.path.join(self.reference_dir,
                                     filtered_collection.file.to_string(
                                         index=False))
            self._ccd = CCDData.read(full_path, unit=u.adu)
            self._recover_lines()
            return self._ccd
        else:
            raise NotImplementedError

    def lamp_exists(self, object, grating, grt_targ, cam_targ):
        filtered_collection = self.ref_lamp_collection[
            (self.ref_lamp_collection['object'] == object) &
            (self.ref_lamp_collection['grating'] == grating) &
            (self.ref_lamp_collection['grt_targ'] == grt_targ) &
            (self.ref_lamp_collection['cam_targ'] == cam_targ)
        ]

        if filtered_collection.empty:
            return False
        elif len(filtered_collection) == 1:
            return True
        else:
            raise NotImplementedError
        
    def check_comp_group(self, comp_group):
        lamps = comp_group.groupby(['object',
                                    'grating',
                                    'grt_targ',
                                    'cam_targ']).size().reset_index(
                                    ).rename(columns={0: 'count'})

        for i in lamps.index:
            if self.lamp_exists(
                    object=lamps.iloc[i]['object'],
                    grating=lamps.iloc[i]['grating'],
                    grt_targ=lamps.iloc[i]['grt_targ'],
                    cam_targ=lamps.iloc[i]['cam_targ']):
                new_group = comp_group[
                    (comp_group['object'] == lamps.iloc[i]['object']) &
                    (comp_group['grating'] == lamps.iloc[i]['grating']) &
                    (comp_group['grt_targ'] == lamps.iloc[i]['grt_targ']) &
                    (comp_group['cam_targ'] == lamps.iloc[i]['cam_targ'])]
                # print(new_group.file)
                return new_group
            else:
                # print(lamps.iloc[i])
                # print(comp_group.file)
                self.log.warning("The target's comparison lamps do not have "
                                 "reference lamps.")
                self.log.debug("In this case a compatible lamp will be obtained"
                               "from all the lamps obtained in the data or"
                               "present in the files.")
        return None

    def _recover_lines(self):
        self.lines_pixel = []
        self.lines_angstrom = []
        pixel_keys = self._ccd.header['GSP_P*']
        for pixel_key in pixel_keys:
            if re.match(r'GSP_P\d{3}', pixel_key) is not None:
                angstrom_key = re.sub('GSP_P', 'GSP_A', pixel_key)
                assert pixel_key[-3:] == angstrom_key[-3:]
                assert angstrom_key in self._ccd.header
                if int(self._ccd.header[angstrom_key]) != 0:
                    self.lines_pixel.append(float(self._ccd.header[pixel_key]))
                    self.lines_angstrom.append(float(self._ccd.header[angstrom_key]))
                else:
                    self.log.debug("File: {:s}".format(self._ccd.header['GSP_FNAM']))
                    self.log.warning(
                        "Ignoring keywords: {:s}={:f}, {:s}={:f}".format(
                        pixel_key,
                        self._ccd.header[pixel_key],
                        angstrom_key,
                        self._ccd.header[angstrom_key]))

        # self.lines_pixel = np.asarray(self.lines_pixel, dtype=float)
        # self.lines_angstrom = np.asarray(self.lines_angstrom, dtype=float)

    def _validate_lines(self):
        """Calls all available validation methods

        Returns:
            True if none of the validation fails.
        """
        assert len(self.lines_pixel) == len(self.lines_angstrom)
        if not self._order_validation(self.lines_pixel):
            return False
        if not self._order_validation(self.lines_angstrom):
            return False
        # if not self._validate_line_existence():
        #     return False
        self._validate_line_existence()
        return True

    @staticmethod
    def _order_validation(lines_array):
        """Checks that the array of lines only increases."""
        previous = None
        for line_value in lines_array:
            # print(line_value)
            if previous is not None:
                try:
                    assert line_value > previous
                    previous = line_value
                except AssertionError:
                    print("Error: Line {:f} is not larger "
                          "than {:f}".format(line_value, previous))
                    return False
            else:
                previous = line_value
        return True

    def _load_nist_list(self, **kwargs):
        """Load all csv files from strong lines in nist."""
        nist_path = kwargs.get(
            'path',
            os.path.join(os.path.dirname(sys.modules['pipeline'].__file__),
                         'data/nist_list'))
        assert os.path.isdir(nist_path)
        nist_files = glob(os.path.join(nist_path, "*.txt"))
        for nist_file in nist_files:
            key = os.path.basename(nist_file)[22:-4]
            nist_data = pandas.read_csv(nist_file, names=['intensity',
                                                          'air_wavelength',
                                                          'spectrum',
                                                          'reference'])
            self.nist[key] = nist_data

    def _validate_line_existence(self):
        """Check if a line actually exists in any table of NIST
        Notes:
            It does not actually check NIST, it loads six csv tables
            from NIST's strong lines for the elements used in lamps.
            Ar,  Cu, Fe, He, Hg, Ne. It does not work perfect so far
            so the it is not checking existence actually but if it finds it
            it will get the value at "spectrum" column in NIST tables which
            correspond to the source of the line for instance Ne I, Ar II.
            """

        lamp_elements = []
        lamp_name = self._ccd.header['OBJECT']
        if len(lamp_name) % 2 == 0:
            for element_index in range(0, len(lamp_name), 2):
                element = lamp_name[element_index:element_index + 2]
                lamp_elements.append(element)

        if self.nist is None:
            self.nist = {}
            self._load_nist_list()

        self.spectrum = list(self.lines_angstrom)
        for i in range(len(self.lines_angstrom)):
            for element in lamp_elements:
                line_info = self.nist[element][
                    self.nist[element].air_wavelength == self.lines_angstrom[i]]
                if line_info.empty:
                    # print(self.lines_angstrom[i], 'no-info')
                    self.spectrum[i] = ''
                else:
                    self.spectrum[i] = line_info['spectrum'].to_string(index=False)
                    # print(self.lines_angstrom[i], line_info['spectrum'].to_string(index=False))
            # print(self.lines_angstrom[i], self.spectrum[i])

