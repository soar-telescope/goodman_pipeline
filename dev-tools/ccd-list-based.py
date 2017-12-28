from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


import argparse
from goodman_ccd.image_processor import ImageProcessor


def get_args(arguments=None):
    parser = argparse.ArgumentParser(description='''Extracts pipeline spectra and
    does wavelength calibration.''')

    parser.add_argument('--bias-list',
                        action='store',
                        default='bias',
                        type=str,
                        metavar='bias list',
                        dest='bias_list',
                        help='Name of the bias list')

    parser.add_argument('--flat-list',
                        action='store',
                        default='flat',
                        type=str,
                        metavar='flat list',
                        dest='flat_list',
                        help='Name of the flat list')

    parser.add_argument('--object-list',
                        action='store',
                        default='object',
                        type=str,
                        metavar='object list',
                        dest='object_list',
                        help='Name of the object list')

    parser.add_argument('--instrument',
                        action='store',
                        default='blue',
                        choices=['red', 'blue'],
                        type=str,
                        metavar='instrument',
                        dest='instrument',
                        help='Name of the instrument or camera')

    parser.add_argument('--technique',
                        action='store',
                        default='spectroscopy',
                        choices=['imaging', 'spectroscopy'],
                        type=str,
                        metavar='technique',
                        dest='technique',
                        help='Name of the technique used')


    args = parser.parse_args(args=arguments)

    return args

if __name__ == '__main__':
    args = get_args()
    print(args)

