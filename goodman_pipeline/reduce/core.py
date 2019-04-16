import argparse
import os


def get_reduce_args(arguments=None):
    parser = argparse.ArgumentParser(
        description='Reduce data live')

    parser.add_argument('--folder',
                        action='store',
                        dest='folder',
                        default=os.getcwd(),
                        help='Data location')

    parser.add_argument('--save-to',
                        action='store',
                        metavar='<save_to>',
                        type=str,
                        default=os.path.join(os.getcwd(), 'RED'),
                        help="Path to reduced data.")

    parser.add_argument('--reset',
                        action='store_true',
                        dest='reset',
                        help='Reset database')

    parser.add_argument('--ignore-bias',
                        action='store_true',
                        dest='ignore_bias',
                        help="Ignore bias correction")

    parser.add_argument('--ignore-flats',
                        action='store_true',
                        dest='ignore_flats',
                        help="Ignore flat field correction")

    parser.add_argument('--saturation',
                        action='store',
                        default=1.,
                        dest='saturation_threshold',
                        metavar='<value>',
                        help="Maximum percent of pixels above saturation "
                             "threshold. Default 1 percent.")

    parser.add_argument('--dbname',
                        action='store',
                        dest='dbname',
                        default='livereductiondb',
                        help='Database to store live reduction data')

    args = parser.parse_args(args=arguments)

    return args
