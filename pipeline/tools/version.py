"""
v1.0.0
    - First public release.
"""

import logging
import requests
import re

from .. import info


LATEST_URL = \
    'https://api.github.com/repos/soar-telescope/goodman/releases/latest'


def get_last(url=LATEST_URL):

    response = requests.get(url)
    tag_name = response.json()['tag_name']

    _version = re.findall(r'\d+', tag_name)

    _api = int(_version[0])
    _feature = int(_version[1])
    _bug = int(_version[2])

    return _api, _feature, _bug


def check_last(version):

    version = re.findall(r'\d+', version)

    api = int(version[0])
    feature = int(version[1])
    bug = int(version[2])

    last_api, last_feature, last_bug = get_last()

    if last_api > api or last_feature > feature or last_bug > bug:

        logging.info('A new version of the Goodman DRP is available: ')

        logging.info('\tcurrent version: v{}.{}.{}'.format(
            api, feature, bug))

        logging.info('\tlatest version: v{}.{}.{}'.format(
            last_api, last_feature, last_bug))

    else:
        logging.info('Goodman DRP {:s}'.format(info.__version__))
