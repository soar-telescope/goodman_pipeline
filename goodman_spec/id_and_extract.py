from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os

import matplotlib
# matplotlib.use('Qt4Agg')
matplotlib.use('Qt4Agg')

import astropy.units as u
import logging
import numpy as np
import matplotlib.pyplot as plt
import re
import glob
import pandas

from astropy.modeling import (models, fitting, Model)
from ccdproc import CCDData
from ccdproc import ImageFileCollection
from goodman_ccd.core import NightDataContainer
from goodman_ccd.core import cosmicray_rejection
from goodman_ccd.core import identify_targets
from goodman_ccd.core import trace_targets
from goodman_ccd.core import get_extraction_zone
from goodman_ccd.core import remove_background
from goodman_ccd.core import ra_dec_to_deg
from goodman_ccd.core import classify_spectroscopic_data
from goodman_spec.wavelength import process_spectroscopy_data
from goodman_spec.redspec import get_args
from numpy import ma
from scipy import signal

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('goodmanccd')

# log.setLevel(level=logging.DEBUG)








# def process_spectroscopy_data(data_container, extraction_type='simple'):
#
#     assert data_container.is_empty is False
#     assert any(extraction_type == option for option in ['simple',
#                                                         'optimal'])
#
#     full_path = data_container.full_path
#
#     for spec_group in data_container.spec_groups:
#         comp_group = None
#         object_group = None
#         comp_ccd_list = []
#         obstypes = spec_group.obstype.unique()
#         if 'COMP' in obstypes:
#             comp_group = spec_group[spec_group.obstype == 'COMP']
#             # print(comp_group)
#         if 'OBJECT' in obstypes:
#             object_group = spec_group[spec_group.obstype == 'OBJECT']
#             # print(object_group)
#         for spec_file in object_group.file.tolist():
#             file_path = os.path.join(full_path, spec_file)
#             ccd = CCDData.read(file_path, unit=u.adu)
#             ccd.header['OFNAME'] = (spec_file, 'Original File Name')
#             if comp_group is not None:
#                 for comp_file in comp_group.file.tolist():
#                     comp_path = os.path.join(full_path, comp_file)
#                     comp_ccd = CCDData.read(comp_path, unit=u.adu)
#                     comp_ccd.header['OFNAME'] = (comp_file,
#                                                  'Original File Name')
#                     comp_ccd_list.append(comp_ccd)
#
#             extracted, comps = manage_extraction(ccd=ccd,
#                                                  extraction=extraction_type,
#                                                  comp_list=comp_ccd_list)
#             if True:
#                 manager = plt.get_current_fig_manager()
#                 if plt.get_backend() == u'GTK3Agg':
#                     manager.window.maximize()
#                 elif plt.get_backend() == u'Qt4Agg':
#                     manager.window.showMaximized()
#                 for edata in extracted:
#                     plt.plot(edata.data, label=edata.header['OBJECT'])
#                     if comps != []:
#                         for comp in comps:
#                             plt.plot(comp.data, label=comp.header['OBJECT'])
#                 plt.legend(loc='best')
#                 plt.show()
#
#     print('\nEND')


if __name__ == '__main__':

    pandas.set_option('display.expand_frame_repr', False)
    args = get_args()

    prefix = 'cfzsto'
    path = '/user/simon/data/soar/work/20161114_eng_3/RED3'

    data_cont = classify_spectroscopic_data(path=path, search_pattern=prefix)

    process_spectroscopy_data(data_container=data_cont, args=args)
