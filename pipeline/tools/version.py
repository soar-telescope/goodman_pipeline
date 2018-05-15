"""
v1.0.0
    - First public release.
"""

import logging
import requests
import re

logger = logging.getLogger(__name__)


LATEST_URL = \
    'https://api.github.com/repos/soar-telescope/goodman/releases/latest'


def get_last(url=LATEST_URL):
    """
    Returns the version of the last release on GitHub.

    Parameters
    ----------
        url (str, optinal) : the URL that is used to retrieve the information

    Returns
    -------
        api (int) : api version.
        feature (int) : feature version.
        bug (int) : bug version.
    """

    response = requests.get(url)

    if response.status_code == 200:

        tag_name = response.json()['tag_name']

        _version = re.findall(r'\d+', tag_name)

        _api = int(_version[0])
        _feature = int(_version[1])
        _bug = int(_version[2])

        return _api, _feature, _bug

    else:
        raise ConnectionRefusedError('Number of tests reached maximum for now.')


def check_last(version):

    version = re.findall(r'\d+', version)

    api = int(version[0])
    feature = int(version[1])
    bug = int(version[2])

    last_api, last_feature, last_bug = get_last()

    return last_api > api or last_feature > feature or last_bug > bug
