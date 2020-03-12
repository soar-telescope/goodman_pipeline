from __future__ import absolute_import

from astropy.io import fits
from unittest import TestCase, skip

from ccdproc import CCDData
from ...core import NightDataContainer
from ..image_processor import ImageProcessor
from ..goodman_ccd import get_args

import numpy as np


class ImageProcessorTest(TestCase):

    def setUp(self):
        arguments = ['--saturation_threshold', '1']
        args = get_args(arguments=arguments)
        data_container = NightDataContainer(path='/fake',
                                            instrument='Red',
                                            technique='Spectroscopy')
        self.image_processor = ImageProcessor(args=args,
                                              data_container=data_container)

        self.ccd = CCDData(data=np.ones((100, 100)),
                           meta=fits.Header(),
                           unit='adu')
        self.ccd.header.set('INSTCONF', value='Red')
        self.ccd.header.set('GAIN', value=1.48)
        self.ccd.header.set('RDNOISE', value=3.89)

        self.half_full_well = 69257

    def test___call__(self):
        self.image_processor()
