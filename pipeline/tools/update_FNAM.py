import glob
import os
import sys
from ccdproc import CCDData
import logging


class KeywordUpdate(object):

    def __init__(self, search_pattern='*fits'):
        self.log = logging.getLogger(__name__)
        self.search_pattern = search_pattern
        self.file_list = None
        self._file_name = None

    def __call__(self, keyword='GSP_FNAM', value='_file',
                 comment='Current file name'):
        if len(sys.argv) == 2:
            self.search_pattern = sys.argv[1]
        elif len(sys.argv) > 2:
            self.log.error("Too many arguments, remember you just need a "
                           "pattern. default is "
                           "{:s}".format(self.search_pattern))
            sys.exit()

        self.file_list = glob.glob(self.search_pattern)

        for self._file_name in self.file_list:
            ccd = CCDData.read(self._file_name, unit='adu')
            if value == '_file':
                value = os.path.basename(self._file_name)
            self._update_keyword(ccd=ccd,
                                 keyword=keyword,
                                 value=value,
                                 comment=comment)

    def _update_keyword(self, ccd, keyword, value, comment):
        if keyword in ccd.header:
            if ccd.header[keyword] == value:
                self.log.info("{:s} is correct".format(keyword))
            else:
                ccd.header.set(keyword,
                               value=value,
                               comment=comment)
        else:
            self.log.info('{:s} does not exists.'.format(keyword))
            ccd.header.set(keyword,
                           value=value,
                           comment=comment)

        ccd.write(self._file_name, overwrite=True)


if __name__ == '__main__':

    updater = KeywordUpdate()
    updater()



