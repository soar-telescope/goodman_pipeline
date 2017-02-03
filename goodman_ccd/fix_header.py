import os
import glob
from astropy.io import fits
from astropy import log

__author__ = 'David Sanmartim'
__date__ = '2016-07-15'
__version__ = "1.0"
__email__ = "dsanmartim@ctio.noao.edu"


def fix_header_and_shape(input_path, output_path, prefix, overwrite=False):
    """
    Remove/Update some  inconvenient parameters in the header of the Goodman FITS
    files. Some of these parameters contain non-printable ASCII characters. The ouptut
    files are created in the output_path. Also convert fits from 3D [1,X,Y] to 2D [X,Y].
    """
    for _file in sorted(glob.glob(os.path.join(input_path, '*.fits'))):

        ccddata, hdr = fits.getdata(_file, header=True, ignore_missing_end=True)
        # 3D to 2D
        ccddata = ccddata[0]
        # keywords to remove
        key_list_to_remove = ['PARAM0', 'PARAM61', 'PARAM62', 'PARAM63', 'NAXIS3', 'INSTRUME']

        # Keyword to be changed (3 --> 2)
        hdr['N_PARAM'] -= len(key_list_to_remove)
        hdr['NAXIS'] = 2

        # Specific keywords to be removed
        for key in key_list_to_remove:
            try:
                if (key in hdr) is True:
                    hdr.remove(keyword=key)
            except KeyError:
                pass

        # Removing duplicated keywords
        key_list = []
        for key in hdr.iterkeys():
            if key in key_list:
                hdr.remove(keyword=key)
            key_list.append(key)

        hdr.add_history('Header and Shape fixed.')
        fits.writeto(os.path.join(output_path, '') + prefix + os.path.basename(_file), ccddata, hdr,
                     clobber=overwrite)
        log.info('Keywords header of ' + os.path.basename(_file) + ' have been updated --> ' + prefix
                 + os.path.basename(_file))
    log.info('Done: all keywords header have been updated.')
    print('\n')
    return
