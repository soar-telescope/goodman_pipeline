# translate from camera angle to wavmode
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import re
import pandas
import logging
import numpy as np
from ccdproc import CCDData
from astropy import units as u
from .core import get_central_wavelength

log = logging.getLogger('goodmanccd.wavmodetranslator')


class SpectroscopicMode(object):

    def __init__(self):
        """Init method for the Spectroscopic Mode

        This method defines a pandas.DataFrame instance that contains all the
        current standard wavelength modes for Goodman HTS.

        """
        columns = ['grating_freq', 'wavmode', 'camtarg', 'grttarg', 'ob_filter']
        spec_mode = [['400', 'm1', '11.6', '5.8', 'None'],
                     ['400', 'm2', '16.1', '7.5', 'GG455'],
                     ['600', 'UV', '15.25', '7.0', 'None'],
                     ['600', 'Blue', '17.0', '7.0', 'None'],
                     ['600', 'Mid', '20.0', '10.0', 'GG385'],
                     ['600', 'Red', '27.0', '12.0', 'GG495'],
                     ['930', 'm1', '20.6', '10.3', 'None'],
                     ['930', 'm2', '25.2', '12.6', 'None'],
                     ['930', 'm3', '29.9', '15.0', 'GG385'],
                     ['930', 'm4', '34.6', '18.3', 'GG495'],
                     ['930', 'm5', '39.4', '19.7', 'GG495'],
                     ['930', 'm6', '44.2', '22.1', 'OG570'],
                     ['1200', 'm0', '26.0', '16.3', 'None'],
                     ['1200', 'm1', '29.5', '16.3', 'None'],
                     ['1200', 'm2', '34.4', '18.7', 'None'],
                     ['1200', 'm3', '39.4', '20.2', 'None'],
                     ['1200', 'm4', '44.4', '22.2', 'GG455'],
                     ['1200', 'm5', '49.6', '24.8', 'GG455'],
                     ['1200', 'm6', '54.8', '27.4', 'GG495'],
                     ['1200', 'm7', '60.2', '30.1', 'OG570'],
                     ['1800', 'Custom', 'None', 'None', 'None'],
                     ['2100', 'Custom', 'None', 'None', 'None'],
                     ['2400', 'Custom', 'None', 'None', 'None']
                     ]
        self.modes_data_frame = pandas.DataFrame(spec_mode, columns=columns)
        # print(self.modes_data_frame)

    def __call__(self,
                 header=None,
                 grating=None,
                 camera_targ=None,
                 grating_targ=None,
                 blocking_filter=None):
        """Get spectroscopic mode

        This method can be called either parsing a header alone or the rest of
        values separated.
        Args:
            header (object): FITS header.
            grating (str): Grating as in the FITS header.
            camera_targ (str): Camera target angle as in the FITS header.
            grating_targ (str): Grating target angle as in the FITS header.
            blocking_filter (str): Order blocking filter as in the FITS header.

        Returns:
            string that defines the instrument wavelength mode.

        """

        if all(x is None for x in (
                grating, camera_targ, grating_targ, blocking_filter)) and \
                        header is not None:

            grating = str(re.sub('[A-Za-z_-]', '', header['grating']))
            camera_targ = str(header['cam_targ'])
            grating_targ = str(header['grt_targ'])
            blocking_filter = str(header['filter2'])
            # print(grating, camera_targ, grating_targ, blocking_filter)

            return self.get_mode(grating=grating,
                                 camera_targ=camera_targ,
                                 grating_targ=grating_targ,
                                 blocking_filter=blocking_filter)
        else:
            grating = re.sub('[A-Za-z_-]', '', grating)

            return self.get_mode(grating=grating,
                                 camera_targ=camera_targ,
                                 grating_targ=grating_targ,
                                 blocking_filter=blocking_filter)

    def get_mode(self, grating, camera_targ, grating_targ, blocking_filter):
        """Get the camera's optical configuration mode.

        This method is useful for data that does not have the WAVMODE keyword

        Args:
            grating (str): Grating frequency as string
            camera_targ (str): Camera target angle as in the header.
            grating_targ (str): Grating target angle as in the header.
            blocking_filter (str): Order blocking filter listed on the header.

        Returns:
            string that defines the wavelength mode used

        """

        # print(grating, camera_targ, grating_targ)
        if any(grat == grating for grat in ('1800', '2100', '2400')):
            central_wavelength = get_central_wavelength(grating=grating,
                                                        grt_ang=grating_targ,
                                                        cam_ang=camera_targ)
            return 'Custom_{:d}nm'.format(int(round(central_wavelength)))

        else:
            _mode = self.modes_data_frame[
                ((self.modes_data_frame['grating_freq'] == grating) &
                 (self.modes_data_frame['camtarg'] == camera_targ) &
                 (self.modes_data_frame['grttarg'] == grating_targ))]
            if _mode.empty:
                central_wavelength = get_central_wavelength(grating=grating,
                                                            grt_ang=grating_targ,
                                                            cam_ang=camera_targ)
                return 'Custom_{:d}nm'.format(int(round(central_wavelength)))
            else:
                # print('%s %s' %(grating, _mode.wavmode))
                # print(_mode['wavmode'].to_string)
                return _mode['wavmode'].to_string(index=False)


if __name__ == '__main__':
    pass




