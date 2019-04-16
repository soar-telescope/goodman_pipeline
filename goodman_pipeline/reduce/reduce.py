import argparse


__version__ = __import__('goodman_pipeline').__version__


def get_args(arguments=None):

    # log.loggin
    parser = argparse.ArgumentParser(
        description='Reduce data live')

    args = parser.parse_args(args=arguments)

    return args


class Reduce(object):

    def __init__(self, arguments=None):
        self.master = {
            'bias': None,
            'flat': None}
        self.args = None
        if self.__validate_args(arguments=arguments):
            self.args = get_args(arguments=arguments)
            self.__validate_args()

    def __validate_args(self, arguments=None):
        if arguments:
            return True
        else:
            print(self.args)
            return True

    def __update_args(self, new_args):
        pass

    def __call__(self, arguments):
        pass

