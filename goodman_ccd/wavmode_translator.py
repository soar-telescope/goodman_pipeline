# translate from camera angle to wavmode

import pandas
from ccdproc import CCDData
from astropy import units as u
import numpy as np
import re
import logging

log = logging.getLogger('goodmanccd.wavmodetranslator')

# modes_dict = {'400': {'m1': {'cam_targ': '11.6', 'grt_targ': '5.8', 'filter2': ''},
#                  'm2': {'cam_targ': '16.1', 'grt_targ': '7.5', 'filter2': 'GG455'}
#                  },
#          '600': {'UV': {'cam_targ': '15.25', 'grt_targ': '7.0', 'filter2': ''},
#                  'Blue': {'cam_targ': '17.0', 'grt_targ': '7.0', 'filter2': ''},
#                  'Mid': {'cam_targ': '20.0', 'grt_targ': '10', 'filter2': 'GG385'},
#                  'Red': {'cam_targ': '27.0', 'grt_targ': '12.0', 'filter2': 'GG495'}
#                  },
#          '930': {'m1': {'cam_targ': '20.6', 'grt_targ': '10.3', 'filter2': ''},
#                  'm2': {'cam_targ': '', 'grt_targ': '', 'filter2': ''},
#                  'm3': {'cam_targ': '29.9', 'grt_targ': '15.0', 'filter2': 'GG385'},
#                  'm4': {'cam_targ': '', 'grt_targ': '', 'filter2': 'GG495'},
#                  'm5': {'cam_targ': '39.4', 'grt_targ': '19.7', 'filter2': 'GG495'},
#                  'm6': {'cam_targ': '44.2', 'grt_targ': '22.1', 'filter2': 'OG570'}
#                  },
#          '1200': {'m0': {'cam_targ': '26.0', 'grt_targ': '16.3', 'filter2': ''},
#                   'm1': {'cam_targ': '29.5', 'grt_targ': '16.3', 'filter2': ''},
#                   'm2': {'cam_targ': '34.4', 'grt_targ': '18.7', 'filter2': ''},
#                   'm3': {'cam_targ': '39.4', 'grt_targ': '20.2', 'filter2': ''},
#                   'm4': {'cam_targ': '44.4', 'grt_targ': '22.2', 'filter2': 'GG455'},
#                   'm5': {'cam_targ': '49.6', 'grt_targ': '24.8', 'filter2': 'GG455'},
#                   'm6': {'cam_targ': '54.8', 'grt_targ': '27.4', 'filter2': 'GG495'},
#                   'm7': {'cam_targ': '60.2', 'grt_targ': '30.1', 'filter2': 'OG570'}
#                   },
#          '1800': {'Custom': {'cam_targ': '', 'grt_targ': '', 'filter2': ''}
#                   },
#          '2100': {'Custom': {'cam_targ': '', 'grt_targ': '', 'filter2': ''}
#                   },
#          '2400': {'Custom': {'cam_targ': '', 'grt_targ': '', 'filter2': ''}
#                   }
#          }


class SpectroscopicMode(object):

    def __init__(self):
        columns = ['grating_freq', 'wavmode', 'camtarg', 'grttarg', 'ob_filter']
        spec_mode = [['400', 'm1', '11.6', '5.8', 'None'],
                     ['400', 'm2', '16.1', '7.5', 'GG455'],
                     ['600', 'UV', '15.25', '7.0', 'None'],
                     ['600', 'Blue', '17.0', '7.0', 'None'],
                     ['600', 'Mid', '20.0', '10.0', 'GG385'],
                     ['600', 'Red', '27.0', '12.0', 'GG495'],
                     ['930', 'm1', '20.6', '10.3', 'None'],
                     ['930', 'm2', 'None', 'None', 'None'],
                     ['930', 'm3', '29.9', '15.0', 'GG385'],
                     ['930', 'm4', 'None', 'None', 'GG495'],
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

    def __call__(self, header=None, grating=None, camera_targ=None, grating_targ=None, blocking_filter=None):
        """

        Args:
            header:
            grating:
            camera_targ:
            grating_targ:
            blocking_filter:

        Returns:

        """

        if all(x is None for x in (grating, camera_targ, grating_targ, blocking_filter)) and header is not None:
            grating = str(re.sub('[A-Za-z_-]', '', header['grating']))
            camera_targ = str(header['cam_targ'])
            grating_targ = str(header['grt_targ'])
            blocking_filter = str(header['filter2'])
            # print(grating, camera_targ, grating_targ, blocking_filter)
            return self.get_mode(grating, camera_targ, grating_targ, blocking_filter)
        else:
            grating = re.sub('[A-Za-z_-]', '', grating)
            return self.get_mode(grating, camera_targ, grating_targ, blocking_filter)

    def get_mode(self, grating, camera_targ, grating_targ, blocking_filter):
        """

        Args:
            grating:
            camera_targ:
            grating_targ:
            blocking_filter:

        Returns:

        """

        # print(grating, camera_targ, grating_targ)
        if any(grat == grating for grat in ('1800', '2100', '2400')):
            return 'Custom' + self.get_central_wavelength(grating, grating_targ, camera_targ)
        _mode = self.modes_data_frame[((self.modes_data_frame['grating_freq'] == grating)
                                                    & (self.modes_data_frame['camtarg'] == camera_targ)
                                                    & (self.modes_data_frame['grttarg'] == grating_targ))]
        # print('%s %s' %(grating, _mode.wavmode))
        # print(_mode['wavmode'].to_string)
        return _mode['wavmode'].to_string(index=False)

    @staticmethod
    def get_central_wavelength(grating, grt_ang, cam_ang):
        """

        Args:
            grating:
            grt_ang:
            cam_ang:

        Returns:

        """

        grating_frequency = float(grating)
        alpha = float(grt_ang)
        beta = float(cam_ang) - float(grt_ang)
        center_wavelength = (1e6 / grating_frequency) * (np.sin(alpha * np.pi / 180.) + np.sin(beta * np.pi / 180.))
        # log.error(center_wavelength)
        return '_{:d}nm'.format(int(round(center_wavelength)))


if __name__ == '__main__':
    spect_mode = SpectroscopicMode()
    # spect_mode(grating='lls_1200', camera_targ='49.6', grating_targ='24.8', blocking_filter='GG455')
    ccd = CCDData.read('/user/simon/data/soar/raw/spectroscopy_engineering_night/2017-02-07/0227_quartz_400M2_GG455.fits', unit=u.adu)
    spect_mode(header=ccd.header)
    # (spect_mode.grating_freq)
    # print(spect_mode.blocking_filter)




