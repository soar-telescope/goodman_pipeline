import logging
import requests
import re

latest_url='https://api.github.com/repos/soar-telescope/goodman/releases/latest'

api = 1
feature = 0
bug = 0

month = 0
year = 0

__str__ = "v{:d}.{:d}.{:d} - {:d}-{:d}".format(api, feature, bug, month, year)


def get_last(url=latest_url):

    response = requests.get(url)
    tag_name = response.json()['tag_name']

    _version = re.findall(r'\d+', tag_name)

    _api = int(_version[0])
    _feature = int(_version[1])
    _bug = int(_version[2])

    return _api, _feature, _bug


def check_last():

    last_api, last_feature, last_bug = get_last()

    if last_api > api or last_feature > feature or last_bug > bug:

        logging.warning('A new version of the Goodman DRP is available: ')

        logging.warning('\tcurrent version: v{}.{}.{}'.format(
            api, feature, bug))

        logging.warning('\tlatest version: v{}.{}.{}'.format(
            last_api, last_feature, last_bug))
