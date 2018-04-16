import glob
import os
import sys
from ccdproc import CCDData

search_pattern = '*fits'

if len(sys.argv) == 2:
    search_pattern = sys.argv[1]
elif len(sys.argv) > 2:
    sys.exit("Too many arguments, remember you just need a pattern. "
             "default is {:s}".format(search_pattern))

file_list = glob.glob(search_pattern)
for file in file_list:
    ccd = CCDData.read(file, unit='adu')
    basename = os.path.basename(file)

    if 'GSP_FNAM' in ccd.header:
        if ccd.header['GSP_FNAM'] == basename:
            print("GSP_FNAM is correct")
        else:
            ccd.header['GSP_FNAM'] = basename
    else:
        print('GSP_FNAM does not exists.')
        ccd.header.set('GSP_FNAM', value=basename, comment='Current file name')
    ccd.write(file, overwrite=True)