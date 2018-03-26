from ccdproc import CCDData
import sys
import glob
import os
import re
import random
import sqlite3


class DatabaseManager(object):

    def __init__(self, db_name='spectlines.db', table_name='line_record'):
        self._conn = sqlite3.connect(db_name)
        self._cursor = self._conn.cursor()
        self.all_files = None
        self._table_fields = None
        self._table_name = table_name
        self.keywords = ['date',
                         'slit',
                         'obstype',
                         'object',
                         'exptime',
                         'obsra',
                         'obsdec',
                         'grating',
                         'cam_targ',
                         'grt_targ',
                         'filter',
                         'filter2',
                         'gain',
                         'rdnoise']

    def __call__(self, data_location=sys.argv[1], pattern='goodman*fits'):
        self.all_files = glob.glob(os.path.join(data_location, pattern))
        self._create_table()
        self._retrieve_information()

    def _retrieve_information(self):
        for lamp in self.all_files:
            ccd = CCDData.read(lamp, unit='adu')
            all_angstrom_key = ccd.header['GSP_A*']
            values = []
            for key in self.keywords:
                values.append(ccd.header[key])
            # print(self.keywords, sep=" ", end="\n", flush=True)
            # print(values, sep=" ", end="\n", flush=True)
            for ang_key in all_angstrom_key:
                ang_value = ccd.header[ang_key]
                pix_key = re.sub('GSP_A', 'GSP_P', ang_key)
                pix_value = ccd.header[pix_key]
                self._fill_table(file_name=os.path.basename(lamp),
                                 values=values,
                                 pix_keyword=pix_key,
                                 pix_value=pix_value,
                                 ang_keyword=ang_key,
                                 ang_value=ang_value)

    def _fill_table(self, file_name,
                    values,
                    pix_keyword,
                    pix_value,
                    ang_keyword,
                    ang_value):

        values_db = [file_name]
        values_db.extend(values)
        values_db.append(pix_keyword)
        values_db.append(pix_value)
        values_db.append(ang_keyword)
        values_db.append(ang_value)
        values_db.append('')
        values_string = "?"
        for i in range(len(values_db) - 1):
            values_string += ",?"
        # print(", ".join(values_db))
        query = '''INSERT into {:s} VALUES ({:s})'''.format(self._table_name,
                                                     values_string)
        self._cursor.execute(query, tuple(values_db))
        self._conn.commit()
        print("Adding record for file {:s}".format(file_name))

    def _create_table(self):
        self._table_fields = ['file_name text']

        sample_ccd = CCDData.read(random.choice(self.all_files), unit='adu')
        for key in self.keywords:

            field_type = type(sample_ccd.header[key])
            field_dbtype = ''
            if 'str' in repr(field_type):
                field_dbtype = 'text'
            elif 'float' in repr(field_type):
                field_dbtype = 'real'

            self._table_fields.append("{:s} {:s}".format(key, field_dbtype))
        self._table_fields.append("pix_key text")
        self._table_fields.append("pix_value real")
        self._table_fields.append("ang_key text")
        self._table_fields.append("ang_value real")
        self._table_fields.append("spectrum text")

        columns_name = ", ".join(self._table_fields)

        query = '''CREATE TABLE IF NOT EXISTS {:s} ({:s})'''.format(
            self._table_name,
            columns_name)
        print(query)
        self._cursor.execute(query)
        self._conn.commit()


if __name__ == '__main__':
    db_manager = DatabaseManager()
    db_manager()
